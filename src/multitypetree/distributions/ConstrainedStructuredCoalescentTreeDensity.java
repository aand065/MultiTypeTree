/*
 * Copyright (C) 2012 Tim Vaughan
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package multitypetree.distributions;

import beast.core.*;
import beast.core.parameter.IntegerParameter;

import java.util.*;

/**
 *
 * @author Bella Anderson
 */
@Description("Likelihood of ColouredTree under structured coalescent.")
public class ConstrainedStructuredCoalescentTreeDensity extends StructuredCoalescentTreeDensity {

      public Input<IntegerParameter> fromTypesInput = new Input<>(
              "fromTypes", "Location.");
      public Input<IntegerParameter> toTypesInput = new Input<>(
               "toTypes", "Location.");
        public Input<Double> migrationTimeInput = new Input<>(
               "migrationTime", "Time of migration.", 0.25);

    protected List<Integer> fromTypes;
    protected List<Integer> toTypes;
    protected double migrationTime;


    @Override
    public void initAndValidate() {
        super.initAndValidate();
        migrationTime = migrationTimeInput.get();

        fromTypes = new ArrayList<>();
        toTypes = new ArrayList<>();
        fromTypes.addAll(Arrays.asList(fromTypesInput.get().getValues()));
        toTypes.addAll(Arrays.asList(toTypesInput.get().getValues()));
    }

    @Override
    public double calculateLogP() {

        // Check validity of tree if required:
        if (checkValidity && !mtTree.isValid())
            return Double.NEGATIVE_INFINITY;

        // Ensure sequence of events is up-to-date:
        updateEventSequence();

        // Start from the tips of the tree, working up.
        logP = 0;

        // Note that the first event is always a sample. We begin at the first
        // _interval_ and the event following that interval.
        for (int eventIdx = 1; eventIdx<eventList.size(); eventIdx++) {

            SCEvent event = eventList.get(eventIdx);
            Integer[] lineageCount = lineageCountList.get(eventIdx);
            double delta_t = event.time - eventList.get(eventIdx - 1).time;
            boolean isValidMigration = true;

               for(int i =0; i< fromTypes.size(); i++){
                 Integer currentFromType = fromTypes.get(i);
                 Integer currentToType = toTypes.get(i);

                   if(currentFromType== event.type && currentToType==event.destType && migrationTime<event.time){
                       isValidMigration=false;
                       break;
                   }

               }

            // Interval contribution:
            if( !isValidMigration ){
                logP=Double.NEGATIVE_INFINITY;
            }else if (delta_t>0) {
                double lambda = 0.0;
                for (int c = 0; c<lineageCount.length; c++) {
                    int k = lineageCount[c];
                    double Nc = migrationModel.getPopSize(c);
                    lambda += k*(k-1)/(2.0*Nc);

                    for (int cp = 0; cp<lineageCount.length; cp++) {
                        if (cp==c)
                            continue;

                        double m = migrationModel.getRate(c, cp);
                        lambda += k*m;
                    }
                }
                logP += -delta_t*lambda;
            }

            // Event contribution:
            switch (event.kind) {
                case COALESCE:
                    double N = migrationModel.getPopSize(event.type);
                    logP += Math.log(1.0/N);
                    break;

                case MIGRATE:
                    double m = 0;
                    if(isValidMigration ){
                        m = migrationModel.getRate(event.type, event.destType);
                    }
                    logP += Math.log(m);
                    break;

                case SAMPLE:
                    // Do nothing here: only effect of sampling event is
                    // to change the lineage counts in subsequent intervals.
                    break;
            }
        }

        return logP;
    }
}