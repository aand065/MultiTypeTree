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
import beast.core.Input;
import beast.evolution.tree.Node;
import beast.util.Randomizer;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Operator implementing a very similar move to ColouredUniform, "
        + " with the difference that this operator draws new node positions"
        + " uniformly between the parent and oldest child heights and"
        + " recolours the modified branches.  Additionally, this operator"
        + " acts on the root node.")
public class ColouredUniformRecolour extends UniformizationRecolourOperator {
    
    public Input<Boolean> rootOnlyInput = new Input<Boolean>("rootOnly",
            "Always select root node for height adjustment.", false);
    
    public Input<Double> rootScaleFactorInput = new Input<Double>("rootScaleFactor",
            "Scale factor used in root height proposals. (Default 1.5)", 1.5);
    
    @Override
    public void initAndValidate() { }
    
    @Override
    public double proposal() {
        
        cTree = colouredTreeInput.get();
        tree = cTree.getUncolouredTree();
        
        // Select internal node to adjust:
        Node node;
        if (rootOnlyInput.get())
            node = tree.getRoot();
        else
            node = tree.getNode(Randomizer.nextInt(tree.getInternalNodeCount()));
                
        // Generate relevant proposal:
        if (node.isRoot())
            return rootProposal(node);
        else
            return nonRootProposal(node);
    }
    
    /**
     * Root node proposal.
     * @param root
     * @return log of HR
     */
    private double rootProposal(Node root) {
        
        double logHR = 0.0;
        
        // Record probability of current colouring:
        logHR += getBranchColourProb(root.getLeft())
                + getBranchColourProb(root.getRight());
        
        // Select new root height:
        double u = Randomizer.nextDouble();
        double f = u*rootScaleFactorInput.get() + (1-u)/rootScaleFactorInput.get();
        double oldestChildHeight = Math.max(
                root.getLeft().getHeight(),
                root.getRight().getHeight());
        root.setHeight(oldestChildHeight + f*(root.getHeight()-oldestChildHeight));
        logHR -= f;
        
        // Select new root node colour:
        setNodeColour(root, Randomizer.nextInt(cTree.getNColours()));
        
        // Recolour branches below root:
        try {
            logHR -= recolourBranch(root.getLeft())
                    + recolourBranch(root.getRight());
        } catch (RecolouringException ex) {
            if (cTree.discardWhenMaxExceeded()) {
                ex.discardMsg();
                return Double.NEGATIVE_INFINITY;
            } else
                ex.throwRuntime();
        }
        
        return logHR;
    }
    
    /**
     * Non-root internal node proposal.
     * @param node
     * @return log of HR of move
     */
    private double nonRootProposal(Node node) {
        
        double logHR = 0.0;
        
        // Record probability of current colouring:
        logHR += getBranchColourProb(node)
                + getBranchColourProb(node.getLeft())
                + getBranchColourProb(node.getRight());
        
        // Select new node height:        
        double upperBound = node.getParent().getHeight();
        double lowerBound = Math.max(
                node.getLeft().getHeight(),
                node.getRight().getHeight());
        node.setHeight(lowerBound+(upperBound-lowerBound)*Randomizer.nextDouble());
        
        // Select new node colour:
        setNodeColour(node, Randomizer.nextInt(cTree.getNColours()));
        
        // Recolour branches connected to node:
        try {
            logHR -= recolourBranch(node)
                    + recolourBranch(node.getLeft())
                    + recolourBranch(node.getRight());
        } catch (RecolouringException ex) {
            if (cTree.discardWhenMaxExceeded()) {
                ex.discardMsg();
                return Double.NEGATIVE_INFINITY;
            } else
                ex.throwRuntime();
        }
        
        return logHR;        
    }
    
}