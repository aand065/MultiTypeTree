/*
 * Copyright (C) 2013 Tim Vaughan <tgvaughan@gmail.com>
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
package beast.evolution.operators;

import beast.core.Description;
import beast.evolution.tree.MultiTypeNode;
import beast.evolution.tree.MultiTypeTree;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.HashSet;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Special move specific to two taxon trees for debugging only.")
public class SpecialTypeBirthDeath extends MultiTypeTreeOperator {
    
    @Override
    public void initAndValidate() { };

    @Override
    public double proposal() {
        mtTree = multiTypeTreeInput.get();
                
        if (mtTree.getLeafNodeCount() != 2)
            throw new IllegalArgumentException("SpecialTypeBirthDeath only valid for 2 taxon trees.");
        
        // Cannot produce or operate on trees without any changes
        if (mtTree.getTotalNumberOfChanges()==0)
            return Double.NEGATIVE_INFINITY;
        
        //String startTree = mtTree.toString();
        
        MultiTypeNode root = (MultiTypeNode)mtTree.getRoot();
        MultiTypeNode left = (MultiTypeNode)root.getLeft();
        MultiTypeNode right = (MultiTypeNode)root.getRight();
        
        // Select event:
        int event = Randomizer.nextInt(mtTree.getTotalNumberOfChanges());
        
        MultiTypeNode node, sister;
        int changeIdx;
        if (event<left.getChangeCount()) {
            node = left;
            sister = right;
            changeIdx = event;
        } else {
            node = right;
            sister = left;
            changeIdx = event - left.getChangeCount();
        }
        
        // Implement birth or death
        double logHR;
        boolean birth = Randomizer.nextBoolean();
        if (birth)
            logHR = birthProposal(node, changeIdx, sister);
        else
            logHR = deathProposal(node, changeIdx, sister);
        
//        if (!birth && (logHR != Double.NEGATIVE_INFINITY) && (Randomizer.nextInt(10)==4)) {
//            PrintStream outFile = null;
//            try {
//                outFile = new PrintStream("death_distrib.txt");
//                outFile.println("# " + startTree);
//                outFile.println("event changeType HR");
//                for (int i=0; i<10000; i++) {
//                    mtTree.restore();
//                    root = (MultiTypeNode)mtTree.getRoot();
//                    left = (MultiTypeNode)root.getLeft();
//                    right = (MultiTypeNode)root.getRight();
//                    
//                    //event = Randomizer.nextInt(mtTree.getTotalNumberOfChanges());
//                    if (event<left.getChangeCount()) {
//                        node = left;
//                        sister = right;
//                        changeIdx = event;
//                    } else {
//                        node = right;
//                        sister = left;
//                        changeIdx = event - left.getChangeCount();
//                    }
//                    deathProposal(node, changeIdx, sister);
//                    
//                    outFile.format("%d %d %g\n", event, node.getChangeType(changeIdx), Math.exp(logHR));
//                }
//                outFile.close();
//                
//                System.exit(0);
//                
//            } catch (FileNotFoundException ex) {
//                Logger.getLogger(SpecialTypeBirthDeath.class.getName()).log(Level.SEVERE, null, ex);
//            } finally {
//                outFile.close();
//            }
//        }
        
        return logHR;
    }

    private double birthProposal(MultiTypeNode node, int changeIdx, MultiTypeNode sister) {
        double logHR = 0.0;

        double tmin, tmax;

        Set<Integer> illegalTypesBirth = new HashSet<Integer>();
        Set<Integer> illegalTypesDeath = new HashSet<Integer>();
        
        if (changeIdx==0)
            illegalTypesDeath.add(node.getNodeType());
        else
            illegalTypesDeath.add(node.getChangeType(changeIdx-1));
        
        illegalTypesBirth.add(node.getChangeType(changeIdx));
        
        tmin = node.getChangeTime(changeIdx);
        
        if (changeIdx+1>=node.getChangeCount()) {
            if (sister.getChangeCount()>0) {
                if (sister.getChangeCount()>1) {
                    illegalTypesDeath.add(sister.getChangeType(sister.getChangeCount()-2));
                    illegalTypesBirth.add(sister.getChangeType(sister.getChangeCount()-2));
                } else {
                    illegalTypesDeath.add(sister.getNodeType());
                    illegalTypesBirth.add(sister.getNodeType());
                }
            } else {
                return Double.NEGATIVE_INFINITY;
            }
            tmax = node.getParent().getHeight();
        } else {
            illegalTypesDeath.add(node.getChangeType(changeIdx+1));
            illegalTypesBirth.add(node.getChangeType(changeIdx+1));
            tmax = node.getChangeTime(changeIdx+1);
        }

        // Backward move probability:
        int Cdeath = mtTree.getNTypes() - illegalTypesDeath.size();
        logHR += Math.log(1.0/(Cdeath*(mtTree.getTotalNumberOfChanges()+1)));
        
        // Forward move probability:
        int Cbirth = mtTree.getNTypes() - illegalTypesBirth.size();
        logHR -= Math.log(1.0/(Cbirth*mtTree.getTotalNumberOfChanges()*(tmax-tmin)));
        
        // Select new change time and colour:
        double tnew = tmin + Randomizer.nextDouble()*(tmax-tmin);
        int n = Randomizer.nextInt(Cbirth);
        int changeType;
        for (changeType = 0; changeType<mtTree.getNTypes(); changeType++) {
            if (illegalTypesBirth.contains(changeType))
                continue;
            if (n==0)
                break;
            n -= 1;
        }
        node.insertChange(changeIdx+1, changeType, tnew);
        
        // Propagate colour changes:
        if (changeIdx+2>=node.getChangeCount()) {
            ((MultiTypeNode)node.getParent()).setNodeType(changeType);
            sister.setChangeType(sister.getChangeCount()-1, changeType);
        }
        
        return logHR;
    }

    private double deathProposal(MultiTypeNode node, int changeIdx, MultiTypeNode sister) {

        // Zero probability of proposing death of absent node.
        if (changeIdx+1>=node.getChangeCount())
            return Double.NEGATIVE_INFINITY;
        
        double logHR = 0.0;

        double tmin, tmax;

        Set<Integer> illegalTypesBirth = new HashSet<Integer>();
        Set<Integer> illegalTypesDeath = new HashSet<Integer>();
        
        if (changeIdx==0)
            illegalTypesDeath.add(node.getNodeType());
        else
            illegalTypesDeath.add(node.getChangeType(changeIdx-1));
        
        illegalTypesBirth.add(node.getChangeType(changeIdx));
        
        tmin = node.getChangeTime(changeIdx);
        
        if (changeIdx+2>=node.getChangeCount()) {
            if (sister.getChangeCount()>0) {
                if (sister.getChangeCount()>1) {
                    illegalTypesDeath.add(sister.getChangeType(sister.getChangeCount()-2));
                    illegalTypesBirth.add(sister.getChangeType(sister.getChangeCount()-2));
                } else {
                    illegalTypesDeath.add(sister.getNodeType());
                    illegalTypesBirth.add(sister.getNodeType());
                }
            } else {
                return Double.NEGATIVE_INFINITY;
            }
            tmax = node.getParent().getHeight();
        } else {
            illegalTypesDeath.add(node.getChangeType(changeIdx+2));
            illegalTypesBirth.add(node.getChangeType(changeIdx+2));
            tmax = node.getChangeTime(changeIdx+2);
        }

        // Backward move probability:
        int Cbirth = mtTree.getNTypes() - illegalTypesBirth.size();
        logHR += Math.log(1.0/(Cbirth*(mtTree.getTotalNumberOfChanges()-1)*(tmax-tmin)));
        
        // Forward move probability:
        int Cdeath = mtTree.getNTypes() - illegalTypesDeath.size();
        logHR -= Math.log(1.0/(Cdeath*mtTree.getTotalNumberOfChanges()));
        
        // Remove change:
        node.removeChange(changeIdx+1);
        
        // Select new colour:
        int n = Randomizer.nextInt(Cdeath);
        int changeType;
        for (changeType = 0; changeType<mtTree.getNTypes(); changeType++) {
            if (illegalTypesDeath.contains(changeType))
                continue;
            if (n==0)
                break;
            n -= 1;
        }
        
        // Update colouring
        node.setChangeType(changeIdx, changeType);
        if (changeIdx+1>=node.getChangeCount()) {
            ((MultiTypeNode)node.getParent()).setNodeType(changeType);
            sister.setChangeType(sister.getChangeCount()-1, changeType);
        }
        
        return logHR;
    }
    
}