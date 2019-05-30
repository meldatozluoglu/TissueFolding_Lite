#ifndef Node_H
#define Node_H

#include <stdio.h>
#include <iostream>
#include <vector>
#include <array>
#include <memory>

class ShapeBase;
/**
 *  The node class
 *  */
class Node{
private:

public:
    Node(int id, int dim, std::array<double,3> pos, int tissuePos, int tissueType);	///<Constructer of the node
	~Node();
    std::array<bool,3>      FixedPos;                       ///< The boolean array stating if the node's position is fixed in any direction, format: [x y z], true = fixed
    int                     Id;								///< The unique identification number of the node
    size_t                  nDim;							///< The number of dimensions of the node, (2 or 3)
    std::array<double,3>    Position;						///< The position array of the node.
    std::array<double,3> 	RKPosition;						///< The position array for position during a Runge-Kutta step array of the node.
    std::array<bool,3>      externalViscositySetInFixing;	///< The boolean array stating if the external viscosity of any axis has been set in node fixing options. The node fixing is carried out before the physical parameter settings in most cases. The boolean check is carried out not to overwrite the existing set viscosity in normal viscosity assignment.
    std::array<double,3> 	externalViscosity;				///< External viscosity of the node, defined by its placement within the tissue. This can be defined as an external adhesion, ECM remodelling, or any other form of viscosity.
    std::array<double,3> 	initialExternalViscosity;		///< External viscosity of the node before any modifications by perturbation (for example ECM stiffness/ viscosity change)
    std::array<double,3> 	maximumExternalViscosity;       ///< Maximum external viscosity the node can reach as a result of modifications by perturbation (for example ECM stiffness/ viscosity change)
    std::array<double,3> 	minimumExternalViscosity;       ///< Minimum external viscosity the node can reach as a result of modifications by perturbation (for example ECM stiffness/ viscosity change)
    std::array<double,3> 	ECMViscosityChangePerHour;		///< The change in ECM viscosity per one hour. The double array of size (1,3) stores the viscosity change in x, y, and z directions, respectively.
    double                  displacement;					///< the displacement of the node from previous time step;
    int                     tissuePlacement;				///< The tissue placement is 0 for basal nodes, 1 for apical nodes, and 2 for middle range
    int                     tissueType;		 				///< The tissue type is 0 for columnar layer, 1 for peripodial membrane, and 2 for linker zone
    bool                    atCircumference;				///< Boolean defining if the node is at the circumference of the columnar layer of the tissue.
    double                  mass;							///< The mass of the node, calculated via the elements that use the node as a vertex
    double                  viscositySurface;				///< The surface of the node, calculated for application of external viscosity via surface. It is positive for a surface that is to feel external viscosity, zero otherwise.
    double                  zProjectedArea;         		///< The surface of the node, as projected in Z, calculated from apical or pasal surfces of elements, lateral surfaces are not included.
    std::vector<int>        immediateNeigs;                 ///< The list of Id's for immediate neighbours of the node, i.e. the nodes that are shared by the elements that utilise the owner of this node.
    std::vector<int>        connectedElementIds;			///< The list of Id's for elements that are utilising this node.
    std::vector<double>     connectedElementWeights;		///< The list of weights (normalised mass) for elements that are utilising this node, order is linked to Node#connectedElementIds.
    bool                    hasLateralElementOwner;			///< The boolean stating if any lateral element uses this node
    bool                    atSymmetricityBorder;			///< The boolean stating if the node is at the border of symmetricity
    bool                    insideEllipseBand;				///< The boolean stating if the node is inside a marker ellipse
    int                     coveringEllipseBandId;			///< The id of the marker ellipse that the node is covered by. If hte node is not covered by any marker ellipse, the value is -1.
    bool                    allOwnersECMMimicing;			///< The boolean stating all the elements making use of the node are ECM mimicking elements.
    bool                    isMaster[3];                    ///< The boolean array stating if the corresponding degree of freedom of the node is a master of any other degree of freedom
    int                     slaveTo[3];                     ///< The int array storing the master node ids if the corresponding degree of freedom of the node is a slave to the respective degree of freedom of another node.
    std::vector<int>        collapsedWith;                  ///< The vector of int storing the ids of the nodes this node is collapsed and bound with.
    int                     adheredTo;                      ///< The id of the node this node is adhered to.
    bool                    onFoldInitiation;               ///< The boolean stating if this node is on a fold initiation region, automatically detected via curvature.
    bool                    checkOwnersforEllipseAsignment; ///< When the ellipse ids are assigned to nodes by owner elemetns, you can end up with elements that have all their nodes engulfed in an ellipse, but the element itself is not assigned into an ellipse. This flag will check for that.
    bool                    positionUpdateOngoing;          ///< The boolean stating if this node has its in the process of updating its position during staged collapse of nodes upon adhesion.
    int                     positionUpdateCounter;          ///< The counter for steps the nodes position update upon collapse with another node has been ongoing.
    bool                    haveCommonOwner(Node* nodeSlave);   ///< The function to check if the node pointed by the input pointer (nodeSlave) has common owner with this node (returns boolean).
    int                     getCommonOwnerId(Node* nodeSlave);  ///< The function to return the ShapeBase#Id the owner element this node shares with the node pointed by the input pointer (nodeSlave).
    void                    setExternalViscosity(double ApicalVisc,double BasalVisc, bool extendExternalViscosityToInnerTissue);    ///< The function to set the viscosity of the node.
    bool                    checkIfNeighbour(int IdToCheck); 	///< The function to check if the node with input Id (IdToCheck) is an immediate neighbour of the node
    bool                    checkIfNodeHasPacking();			///< The function to check if the node is eligible for packing.
    std::array<double,3>    getCurrentPosition();				///< return the current position of the node
    void                    displayConnectedElementIds();		///< This function will print out a list of connected element Id's
    void                    displayConnectedElementWeights();	///< This function will print out the weights of the connected elements, in the order of  Id s given in connectedElementIds
    void                    displayPosition();                  ///< The helper function to display the position of the node during debugging or in catasthropic failure.
    void                    addToImmediateNeigs(int newNodeId);	///< This function adds the input node Id (int) to the list of neighbours of the node (Node#immediateNeigs)
    void                    addToConnectedElements(int newElementId, double volumePerNode);	///< This function adds the input newElementId (int) to the list of elements connected by this node, updating the mass, and weights of mass per connected element in the process.
    void                    removeFromConnectedElements(int newElementId, double volumePerNode);///< This function removes the input newElementId (int) from the list of elements connected by this node, updating the mass, and weights of mass per connected element in the process.
    bool                    isMyNeig(int nodeId);               ///< The function to check if the node with the input node id is a neighbour of this node (returns boolean).
    bool                    isNeigWithMyCollapsedNodes(int NodeId, const std::vector <std::unique_ptr<Node>>& Nodes);   ///< The function to check if the node with the input Node#Id (nodeId) is collapsed with a neighbour of this node (returns boolean).
    void                    getNewCollapseListAndAveragePos(std::vector<int> &newCollapseList, std::array<double,3> avrPos, std::array<bool,3> fix, const std::vector<std::unique_ptr<Node>> &Nodes, int masterNodeId); ///< This function appends the list of node Ids this node is collapsed with to the input array.
    void                    clearDuplicatesFromCollapseList();  ///< This function clears the Node#collapsedWith array from duplicates.
    void                    collapseOnNode(const std::vector<std::unique_ptr<Node>>& Nodes, int masterNodeId); ///< This function collapses this node on the node with the Node#Id with input masterNodeId, at one time step.
    void                    collapseOnNodeInStages( std::vector<int> &newCollapseList, std::array<double,3> avrPos, std::array<bool,3> fix, const std::vector<std::unique_ptr<Node>>& Nodes); ///< This function collapses this node on the node with the Node#Id with input masterNodeId, at multiple time steps, perturbing the position of the node smoothly to ensure stability.
    void                    updatePositionTowardsPoint(std::array<double,3> avrPos,std::array<bool,3> fix); ///< This function moves this node towards the coordinates given in input array avrPos, excluding movement of fixed position.
    bool                    isECMChangeAppliedToNode(bool changeApicalECM, bool changeBasalECM, std::vector<int> &ECMChangeEllipseBandIds, int numberOfECMChangeEllipseBands); ///< This function checks if the node falls within the range of the marker ellipse band ids for the ECM property perturbations.
    int                     getId();                            ///< The function returns the Id (Node#Id) of the node
};
#endif
