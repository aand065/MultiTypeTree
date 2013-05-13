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
import beast.evolution.tree.Node;
import beast.util.Randomizer;
import java.util.HashSet;
import java.util.Set;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Migration birth/death operator from Ewing et al., 2004.")
public class TypeBirthDeath extends MultiTypeTreeOperator {

    private Set<Integer> illegalTypes;
    
    @Override
    public void initAndValidate() {
        illegalTypes = new HashSet<Integer>();
    }
    
    @Override
    public double proposal() {
        
        mtTree = multiTypeTreeInput.get();
        
        // Immediate reject if <3 types in model
        if (mtTree.getNTypes()<3)
            return Double.NEGATIVE_INFINITY;
        
        // Select event at random:
        int event = Randomizer.nextInt(mtTree.getTotalNumberOfChanges()
                + mtTree.getInternalNodeCount()-1);
        
        MultiTypeNode node = null;
        int changeIdx = -1;
        if (event < mtTree.getInternalNodeCount()-1)
            node = (MultiTypeNode)mtTree.getNode(event + mtTree.getLeafNodeCount());           
        
        else {
            event -= mtTree.getInternalNodeCount()-1;
            for (Node thisNode : mtTree.getNodesAsArray()) {
                if (thisNode.isRoot())
                    continue;
                
                if (event < ((MultiTypeNode)thisNode).getChangeCount()) {
                    node = (MultiTypeNode)thisNode;
                    changeIdx = event;
                    break;
                }
                
                event -= ((MultiTypeNode)thisNode).getChangeCount();
            }
        }

        double logHR;
        boolean birth = Randomizer.nextBoolean();
        if (birth)
            logHR = birthMove(node, changeIdx);
        else
            logHR = deathMove(node, changeIdx);
        
        return logHR;
    }
    
    private double birthMove(MultiTypeNode node, int changeIdx) {
        double logHR = 0.0;
        
        // Construct the set of illegal change types for the reverse move:
        try {
            getIllegalTypes(changeIdx, node);
        } catch (Exception ex) {
            // Subtree contains leaf node:
            return Double.NEGATIVE_INFINITY;
        }
        
        // Record number of legal change types in reverse move for HR:
        int Cdeath = mtTree.getNTypes() - illegalTypes.size();
        if (Cdeath == 0)
            // Reverse move is impossible
            return Double.NEGATIVE_INFINITY;
        
        // Reverse move HR contribution:
        logHR += Math.log(1.0/(Cdeath
                *(mtTree.getTotalNumberOfChanges()+1 + mtTree.getInternalNodeCount()-1)));
        
        // Determine time boundaries for new event:
        double tmin = changeIdx<0
                ? node.getHeight()
                : node.getChangeTime(changeIdx);
        double tmax = changeIdx+1>=node.getChangeCount()
                ? node.getParent().getHeight()
                : node.getChangeTime(changeIdx+1);
        
        double tnew = tmin + (tmax-tmin)*Randomizer.nextDouble();
        
        // Insert new change with dummy type:
        node.insertChange(changeIdx+1, 0, tnew);
        
        // Construct the set of illegal change types for the forward move:
        try {
            getIllegalTypes(changeIdx+1, node);
        } catch (Exception ex) {
            return Double.NEGATIVE_INFINITY;
        }
        
        // Record number of legal change types in forward move for HR:
        int Cbirth = mtTree.getNTypes() - illegalTypes.size();
        
        // No legal moves
        if (Cbirth == 0)
            return Double.NEGATIVE_INFINITY;
       
        // Select change type:
        int changeType = selectLegalChangeType();
        
        // Propagate changes across subtree:
        retypeSubtree(changeIdx+1, node, changeType);
        
        // Forward move HR contribution:
        logHR -= Math.log(1.0/(Cbirth
                *(mtTree.getTotalNumberOfChanges()-1+mtTree.getInternalNodeCount()-1)
                *(tmax-tmin)));
        
        return logHR;
    }
    
    private double deathMove(MultiTypeNode node, int changeIdx) {
        double logHR = 0.0;
        
        // Reject if event above (node,changeIdx) is not a migration
        if (changeIdx+1>=node.getChangeCount())
            return Double.NEGATIVE_INFINITY;
        
        // Construct set of illegal change types for reverse move:
        try {
            getIllegalTypes(changeIdx+1, node);
        } catch (Exception ex) {
            return Double.NEGATIVE_INFINITY;
        }
        
        // Record number of legal change types for reverse move HR
        int Cbirth = mtTree.getNTypes() - illegalTypes.size();
        if (Cbirth==0)
            // Reverse move impossible
            return Double.NEGATIVE_INFINITY;
        
        double tmin = changeIdx<0
                ? node.getHeight()
                : node.getChangeTime(changeIdx);
        double tmax = changeIdx+2>=node.getChangeCount()
                ? node.getParent().getHeight()
                : node.getChangeTime(changeIdx+2);
        
        // Reverse move HR contribution
        logHR += Math.log(1.0/(Cbirth
                *(mtTree.getTotalNumberOfChanges()-1+mtTree.getInternalNodeCount()-1)
                *(tmax-tmin)));
        
        // Remove dying type change WITHOUT MODIFYING NODE COLOURS:
        // (Must do this here, otherwise the illegal types calculation will
        // include the type of the dying move.  This means that changeIdx+1
        // now refers to the change _above_ the removed change.)
        node.removeChange(changeIdx+1);
        
        // Construct set of illegal change types for forward move:
        try {
            getIllegalTypes(changeIdx, node);
        } catch (Exception ex) {
            // Subtree contains leaf
            return Double.NEGATIVE_INFINITY;
        }
        
        // Record number of legal change types in forward move for HR
        int Cdeath = mtTree.getNTypes() - illegalTypes.size();
        
        // Reject if no legal change types:
        if (Cdeath == 0)
            return Double.NEGATIVE_INFINITY;
        
        // Select legal change type:
        int changeType = selectLegalChangeType();

        // Implement subtree type changes
        retypeSubtree(changeIdx, node, changeType);
        
        // Forward move HR contribution
        logHR -= Math.log(1.0/(Cdeath
                *(mtTree.getTotalNumberOfChanges()+1+mtTree.getInternalNodeCount()-1)));        
        
        return logHR;
    }
    
    /**
     * Populates illegalTypes set with types of subtree containing node and no
     * intervening migration events.
     * 
     * @param changeIdx Position of event r.
     * @param node Node on which event r sits.
     * @throws Exception when subtree contains a leaf node.
     */
    public void getIllegalTypes(int changeIdx, MultiTypeNode node) throws Exception {
        illegalTypes.clear();
        
        if (changeIdx<0) {
            MultiTypeNode startNode = findDecendentNodeWithMigration((MultiTypeNode)node.getLeft());
            
            if (startNode.getChangeCount()==0)
                // Termiated at leaf: move would be irreversable
                throw new Exception("Leaf in sub-tree.");
            
            if (startNode.getChangeCount()==1)
                illegalTypes.add(startNode.getNodeType());
            else
                illegalTypes.add(startNode.getChangeType(startNode.getChangeCount()-2));
            
            getIllegalTypesRecurse(illegalTypes,
                    (MultiTypeNode)startNode.getParent(), startNode);
            
        } else {
            if (changeIdx==0) {
                illegalTypes.add(node.getNodeType());
            } else {
                illegalTypes.add(node.getChangeType(changeIdx-1));
            }
            
            if (changeIdx+1<node.getChangeCount()) {
                illegalTypes.add(node.getChangeType(changeIdx+1));
            } else {
                getIllegalTypesRecurse(illegalTypes,
                        (MultiTypeNode)node.getParent(), node);
            }
        }
    }
    
    /**
     * Recursive method used by getIllegalTypes.
     * 
     * @param node
     * @param prevNode 
     */
    private void getIllegalTypesRecurse(Set<Integer> illegalTypes,
            MultiTypeNode node, MultiTypeNode prevNode) throws Exception {
        
        if (node.isLeaf())
            throw new Exception("Leaf in sub-tree.");
        
        if (prevNode == node.getParent()) {
            
            MultiTypeNode left = (MultiTypeNode)node.getLeft();
            MultiTypeNode right = (MultiTypeNode)node.getRight();
            
            if (left.getChangeCount()>0) {
                if (left.getChangeCount()>1)
                    illegalTypes.add(left.getChangeType(left.getChangeCount()-2));
                else
                    illegalTypes.add(left.getNodeType());
                
            } else
                getIllegalTypesRecurse(illegalTypes, left, node);
            
            if (right.getChangeCount()>0)
                if (right.getChangeCount()>1)
                    illegalTypes.add(right.getChangeType(right.getChangeCount()-2));
                else
                    illegalTypes.add(right.getNodeType());
            else
                getIllegalTypesRecurse(illegalTypes, right, node);
            
        } else {
            
            MultiTypeNode sister = (MultiTypeNode)getOtherChild(node, prevNode);
            
            if (!node.isRoot()) {
                if (node.getChangeCount()>0)
                    illegalTypes.add(node.getChangeType(0));
                else
                    getIllegalTypesRecurse(illegalTypes, (MultiTypeNode)node.getParent(), node);
            }

            if (sister.getChangeCount()>0) {
                if (sister.getChangeCount()>1)
                    illegalTypes.add(sister.getChangeType(sister.getChangeCount()-2));
                else
                    illegalTypes.add(sister.getNodeType());
            } else
                getIllegalTypesRecurse(illegalTypes, sister, node);
        }
        
    }
    
    /**
     * Retype subtree including edge from changeIdx to changeIdx+1 on node
     * with new type changeType.
     * 
     * @param changeIdx
     * @param node
     * @param changeType 
     */
    private void retypeSubtree(int changeIdx, MultiTypeNode node, int changeType) {
        if (changeIdx<0) {
            MultiTypeNode startNode = findDecendentNodeWithMigration((MultiTypeNode)node.getLeft());
            startNode.setChangeType(startNode.getChangeCount()-1, changeType);
            retypeSubtreeRecurse(changeType, (MultiTypeNode)startNode.getParent(), startNode);
        } else {
            node.setChangeType(changeIdx, changeType);
            if (changeIdx+1>=node.getChangeCount())
                retypeSubtreeRecurse(changeType, (MultiTypeNode)node.getParent(), node);
        }
    }
    
    /**
     * Assigns a new type to all branches within the subtree attached to node.
     * 
     * @param type
     * @param node
     * @param prevNode 
     */
    private void retypeSubtreeRecurse(int type, MultiTypeNode node, MultiTypeNode prevNode) {
        if (node.isLeaf())
            throw new IllegalArgumentException("Leaf attached to subtree: cannot retype!");
        
        node.setNodeType(type);
        
        if (prevNode == node.getParent()) {
                        
            MultiTypeNode left = (MultiTypeNode)node.getLeft();
            MultiTypeNode right = (MultiTypeNode)node.getRight();
            
            if (left.getChangeCount()>0)
                left.setChangeType(left.getChangeCount()-1, type);
            else
                retypeSubtreeRecurse(type, left, node);
            
            if (right.getChangeCount()>0)
                right.setChangeType(right.getChangeCount()-1, type);
            else
                retypeSubtreeRecurse(type, right, node);
            
        } else {
            
            MultiTypeNode sister = (MultiTypeNode)getOtherChild(node, prevNode);
            
            if (!node.isRoot() && node.getChangeCount()==0)
                retypeSubtreeRecurse(type, (MultiTypeNode)node.getParent(), node);

            if (sister.getChangeCount()>0)
                sister.setChangeType(sister.getChangeCount()-1, type);
            else
                retypeSubtreeRecurse(type, sister, node);
            
        }
            
    }
    
    /**
     * Does what it says.
     * 
     * @param node
     * @return leaf node or node with a migration event on it.
     */
    private MultiTypeNode findDecendentNodeWithMigration(MultiTypeNode node) {
        if (node.isLeaf() || node.getChangeCount()>0)
            return node;
        else
            return findDecendentNodeWithMigration((MultiTypeNode)node.getLeft());
    }
    
    
    /**
     * Choose a change type uniformly at random from the types not found
     * in the this.illegalTypes set.
     * 
     * @return chosen type
     */
    private int selectLegalChangeType() {
        int n = Randomizer.nextInt(mtTree.getNTypes()- illegalTypes.size());
        int changeType;
        for (changeType=0; changeType<mtTree.getNTypes(); changeType++) {
            if (illegalTypes.contains(changeType))
                continue;
            if (n==0)
                break;
            n -= 1;
        }
        
        return changeType;
    }
    
}