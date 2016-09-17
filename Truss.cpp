#include "Truss.h"
#include <cmath>

Truss::Truss (int nn, int nm, vector<double> nc, vector<int> ncon, vector<int> ndof, vector<double> nl, vector<double> mcs, vector<double> mym){
	numNodes = nn;
	numMembers = nm;
	
	Matrix Node_Coor(nn, 3); //Node coordinates in 3-d space (x,y,z)
	Matrix Node_Conn(nn, 2);
	Matrix Node_Freedom(nn, 3);
	Matrix Node_Loads(nn, 3); //Have to figure out how to place the load (at the three force member)
		
	Matrix Member_CS(nm, 1); //Cross sections of all members
	Matrix Member_YM(nm, 1); //Young's modulous of all memebers
	
	//Initialized to zeros matrices
	Matrix ForceMatrix(nm, 1);
	Matrix DisplacementMatrix(nn, 3);
	Matrix ReactionMatrix(nn, 3);
	
	int count =0;
	for (int i=0; i<nn; i++){
		for (int j=0; j<3; j++){
		Node_Coor(i, j) = nc.at(count);
		Node_Freedom (i,j) = ndof.at(count); //1->fixed, 0->free
		Node_Loads(i,j) = nl.at(count);
		if(j<2)
			Node_Conn(i,j) = ncon.at(count);
		count ++; 
		}
	}
	
	count =0;
	for(i=0; i<nm; i++){
		for (j=0; j<2; j++){
			Member_CS(i,j) = mcs.at(count);
			Member_YM(i,j) = mym.at(count);
			count ++;
		}
	}
}

void Truss::PopulateMatrices(){
	Matrix stiffnessMatrix (3*this.numNodes, 3* this.numNodes);
	Matrix displacementMatrix = this.Node_Freedom;
	for(int i=0; i<this.numNodes; i++)
		for(int j=0; j<3; j++)
			displacementMatrix (i,j) = this.Node_Freedom(i,j) == 1 ? 0 : 1;
	
}
//Joint::Joint (double x_coor, double y_coor){
//	m_x = x_coor;
//	m_y = y_coor;
//	m_z = x_coor;
//	m_appliedForceX = 0;
//	m_appliedForceY = 0;
//	m_appliedForceZ = 0;
//}
//
//double Joint::getX() const{
//	return m_x;
//}
//
//double Joint::getY() const{
//	return m_y;
//}
//
//double Joint::getZ() const{
//	return m_z;
//}
//
//double Joint::numMembers() const{
//	return members.size();
//}
//
//double Joint::appliedX() const{
//	return appliedForceX;
//}
//
//double Joint::appliedY() const{
//	return appliedForceY;
//}
//
//double Joint::appliedz() const{
//	return appliedForceZ;
//}
//void Joint::setX(double x_coor){
//	m_x = x_coor;
//}
//
//void Joint::setY(double y_coor){
//	m_y = y_coor;
//}
//
//void Joint::setZ(double z_coor){
//	
//}
//
//void Joint::setAppliedForce (double x_applied, double y_applied, double z_applied){
//	m_appliedForceX = x_applied;
//	m_appliedForceY = y_applied;
//	m_appliedForceZ = z_applied;
//}
//
////Member Class
//
//int Member::m_nextID = 0; //Start member numbering at one
//
//Member::Member(Joint* left, Joint* right){
//	leftParent = left;
//	rightParent = right;
//	m_id = m_nextID;
//	m_nextID ++;
//	force = 0;
//}
//
//double Member::length() const{
//	deltax = leftParent->getX() - rightParent->getX();
//	deltay = leftParent->getY() - rightParent->getY();
//	
//	return (sqrt(deltax*deltax + deltay*deltay)); //Return XY distance between the two joints
//}
//
//int Member::id() const{
//	return m_id;
//}
//
//Joint* Member::leftJoint() const{
//	return leftParent;
//}
//
//Joint* Member::rightJoint() const{
//	return rightParent;
//}
