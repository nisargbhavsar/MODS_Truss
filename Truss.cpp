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
		//for (j=0; j<2; j++){
		//	Member_CS(i,j) = mcs.at(count);
		//	Member_YM(i,j) = mym.at(count);
		//	count ++;
		//}
		Member_CS(i, 0) = mcs.at(count);
		Member_YM(i, 0) = mym.at(count);
		count++;
	}
}

void Truss::PopulateMatrices(){
	Matrix stiffnessMatrix (3*this.numNodes, 3* this.numNodes);
	Matrix displacementMatrix = this.Node_Freedom; //Node nmber by row, xyz by coloumn
	for(int i=0; i<this.numNodes; i++)
		for(int j=0; j<3; j++)
			displacementMatrix (i,j) = this.Node_Freedom(i,j) == 1 ? 0 : 1;
	
	//
	Matrix Tj (this.numMembers, 3); //Stiffness vector
	//Loop through the number of members
	for(i =0; i<this.numMembers; i++)
	{
		int firstNode = this.Node_Conn(i, 0);
		int secondNode = this.Node_Conn(i,1);
		deltax = this.Node_Coor(firstNode, 0) - this.Node_Coor(secondNode,0);
		deltay = this.Node_Coor(firstNode, 1) - this.Node_Coor(secondNode,1);
		deltaz = this.Node_Coor(firstNode, 2) - this.Node_Coor(secondNode,2);
		
		Matrix delta (3,1);
		delta.setElement(0,0, deltax);
		delta.setElement(1,0, deltay);
		delta.setElement(2,0, deltaz);
		
		double length = norm(delta);
		
		Matrix unitVector = delta;
		unitVector = (1/length) * unitVector;
		
		Matrix m_stiffness = unitVector * unitVector.transpose();
		G = this.Member_CS(i, 0) * this.Member_YM(i, 0) / length; //Stiffness constant (Cross-Sectional Area * Young's Modulous / Length)
		
		Tj.setElement(i,0, delta.getElement(0,0) * G);
		Tj.setElement(i,1, delta.getElement(1,0) * G);
		Tj.setElement(i,2, delta.getElement(2,0) * G);
		
		vector<int> e;
		for(int temp = 3* firstNode -3; temp < 3*firstNode-1; temp++)
			e.push_back(temp);
		for(temp = 3*secondNode-3; temp < 3*secondNode-1; temp++)
			e.push_back(temp);
		Matrix s_n (6,6);
		for (temp =0; temp<6; temp++){
			for(int temp2=0; temp2<6; temp2++){
				double element = m_stiffness.getElement(temp%3, temp2%3);
				if((temp<3 && temp2<3) || (temp>2 && temp2 >2))
					{}//double element = m_stiffness.getElement(temp%3, temp2%3);
				if((temp<3 && temp2>2) || (temp >2 && temp2<3))
					element *= -1; //double element = -1 * m_stiffness.getElement(temp%3, temp2%3);
			s_n.setElement(temp, temp2, element);
			}
		}
		s_n *= G;
		for(temp =0; temp<e.size(); temp++)
		{
			double temp1 = stiffnessMatrix.getElement(e.at(temp), e.at(temp))
			stiffnessMatrix.setElement(e.at(temp), e.at(temp), (temp1 + s_n.getElement(temp, temp)));
			s_n.setElement(0,0) = m_stiffness;
			s_n.setElement(0,1) = -1*m_stiffness;
			s_n.setElement(1,0) = -1*m_stiffness;
			s_n.setElement(1,1) = m_stiffness;
			temp += G *  
		}
   S(e,e)=S(e,e)+G*[s -s;-s s];
               % add this members stiffness to stiffness matrix
	}
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
