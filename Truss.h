#include <vector>
#include "matrix.h"

//class Member;
//class Joint;

class Truss {
	public:
		
		Truss(int num_node, int num_member, vector <double> node_coor, vector<int> node_conn, vector<int> node_dof, vector<double> node_l, vector<double> member_cs, vector<double> member_ym);
		void populateMatrices();
		
		 
//		double getMass();
//		Matrix getForceMatrix();
//		double getFitness();
		int numNodes, numMembers;
		
		
	private:
//		std::vector<Joint*> moveableJoints;
//		std::vector<Joint*> fixedJoints; 
//		std::vector<Member*> members;
		Matrix Node_Coor; //Node coordinates in 3-d space (x,y,z)
		Matrix Node_Conn;
		Matrix Node_Freedom;
		Matrix Node_Loads; //Have to figure out how to place the load (at the three force member)
		
		Matrix Member_CS; //Cross sections of all members
		Matrix Member_YM; //Young's modulous of all memebers
		
		Matrix ForceMatrix;
		Matrix DisplacementMatrix;
		Matrix ReactionMatrix;
}
//class Joint {
//public:
//  Joint(double x, double y, double z); //(x,y, z) -> Points in 3d space
//  
//  double getX()const; 
//  double getY()const;
//  double getZ()const;
//  int numMembers()const;
//  double appliedX()const;
//  double appliedY()const;
//  double appliedZ()const;
//  
//  void setX(double x);
//  void setY(double y);
//  void setZ(double z);
//  void setAppliedForce(double x, double y, double z);
//  
//
//  std::vector<Member*> members;
//
//private:
//  double m_x;
//  double m_y;
//  double m_z;
//  double m_appliedForceX;  
//  double m_appliedForceY;
//  double m_appliedForceZ;
//
//};
//
//class Member {
//private:
//  Joint* leftParent;
//  Joint* rightParent;
//  static int m_nextId;
//  int m_id;
//public:
//  double force;
//  Member(Joint* left, Joint* right);
//  double length()const;
//  int id()const;
//  Joint* leftJoint()const;
//  Joint* rightJoint()const;
//};
