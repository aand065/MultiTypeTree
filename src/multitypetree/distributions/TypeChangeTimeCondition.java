package multitypetree.distributions;

import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.MultiTypeNode;
import beast.evolution.tree.Node;

import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Conditions AGAINST specific type changes in given region.")
public class TypeChangeTimeCondition extends MultiTypeTreeDistribution {

    public Input<RealParameter> h1Input = new Input<>(
            "h1",
            "Smaller height in height range.",
            new RealParameter(new Double[] {0.0}));

    public Input<RealParameter> h2Input = new Input<>(
            "h2",
            "Larger height in height range.",
            new RealParameter(new Double[] {Double.POSITIVE_INFINITY}));

    public Input<IntegerParameter> fromTypesInput = new Input<>(
            "fromTypes",
            "Specify type changes FROM these types.");

    public Input<String> fromTypeNamesInput = new Input<>(
            "fromTypeNames",
            "Specify comma-delimited list of FROM type names.",
            Input.Validate.XOR, fromTypesInput);

    public Input<IntegerParameter> toTypesInput = new Input<>(
            "toTypes",
            "Specify type changes TO these types.",
            Input.Validate.REQUIRED);

    public Input<String> toTypeNamesInput = new Input<>(
            "toTypeNames",
            "Specify comma-delimited list of FROM type names.",
            Input.Validate.XOR, toTypesInput);

    protected Set<Integer> fromTypes, toTypes;
    protected boolean needsUpdate;

    protected RealParameter h1, h2;

    @Override
    public void initAndValidate() throws Exception {
        super.initAndValidate();

        h1 = h1Input.get();
        h2 = h2Input.get();

        fromTypes = new HashSet<>();
        toTypes = new HashSet<>();

        needsUpdate = true;
    }

    @Override
    public double calculateLogP() throws Exception {
        update();

        logP = 0.0;

        for (Node node : mtTree.getNodesAsArray()) {
            if (node.isRoot())
                continue;

            if (node.getHeight()>h2.getValue()
                || node.getParent().getHeight()<h1.getValue())
                continue;

            MultiTypeNode mtNode = (MultiTypeNode)node;

            int lastType = mtNode.getNodeType();
            for (int i=0; i<mtNode.getChangeCount(); i++) {
                if (mtNode.getChangeTime(i)>h1.getValue()
                        && mtNode.getChangeTime(i)<h2.getValue()
                        && toTypes.contains(lastType)
                        && fromTypes.contains(mtNode.getChangeType(i))) {
                    logP = Double.NEGATIVE_INFINITY;
                    return logP;
                }

                lastType = mtNode.getChangeType(i);
            }
       }

        return logP;
    }

    private void update() {
        if (!needsUpdate)
            return;

        fromTypes.clear();
        if (fromTypesInput.get() != null) {
            fromTypes.addAll(Arrays.asList(fromTypesInput.get().getValues()));
        } else {
            for (String typeName : fromTypeNamesInput.get().split(",")) {
                int type = mtTree.getTypeFromString(typeName.trim());
                if (type<0)
                    throw new IllegalArgumentException("Type name '" + typeName
                            + "' specified in input to TypeChangeTimeCondition not found.");

                fromTypes.add(type);
            }
        }

        toTypes.clear();
        if (toTypesInput.get() != null) {
            toTypes.addAll(Arrays.asList(toTypesInput.get().getValues()));
        } else {
            for (String typeName : toTypeNamesInput.get().split(",")) {
                int type = mtTree.getTypeFromString(typeName.trim());
                if (type<0)
                    throw new IllegalArgumentException("Type name '" + typeName
                            + "' specified in input to TypeChangeTimeCondition not found.");

                toTypes.add(type);
            }
        }

        needsUpdate = false;
    }

    @Override
    protected boolean requiresRecalculation() {
        if ((fromTypesInput.get() != null && fromTypesInput.get().somethingIsDirty())
            || (toTypesInput.get() != null && toTypesInput.get().somethingIsDirty()))
            needsUpdate = true;

        return true;
    }
}
