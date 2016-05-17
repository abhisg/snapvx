#ifndef ADMM_H
#define ADMM_H

#include "Solver/Node.hpp"
#include "Solver/Edge.hpp"
#include "Solver/NodeVar.hpp"
#include "Solver/ProximalMap.hpp"
#include "CVXcanon/include/Eigen/Sparse"
#include "CVXcanon/include/Eigen/Core"
#include <unordered_map>
#include <string>
#include <utility>
#include <unordered_map>
#include <thread>

class ADMM
{
	protected:
		/*std::vector<LinOp* > node_objectives;
		std::vector<LinOp* > edge_objectives;
		std::vector<std::vector< LinOp *> > node_constraints;
		std::vector<std::vector< LinOp *> > edge_constraints;*/

		//list of nodes and edges for the problem
		std::vector<Node *> node_list;
		std::vector<Edge *> edge_list;
		
		//update variables
		std::unordered_map<int,Node_Var> node_x_vals;
		std::unordered_map<int,Eigen::MatrixXd> edge_z_vals;
		std::unordered_map<int,Eigen::MatrixXd> edge_u_vals;
		int x_var_size;
		int z_var_size;

		std::map<std::string,double > solver_options;

		//supplementary variables
		double size_x;
		double size_z;
		double rho;
		double lambda;
		double eps_abs;
		double eps_rel;

		//solver functions
		//void ADMM_node(Node *);
		//void ADMM_edge(Edge *,double&,double&,double&,double&,double&);
	public:
		ADMM();
		//void LoadNodes(std::vector<LinOp* > &,std::vector<std::vector< LinOp *> > &);
		//void LoadEdges(std::vector<LinOp* > &,std::vector<std::vector< LinOp *> > &);
		void LoadNodeProximal(std::string, std::vector<int> &,std::vector<std::string>&,
			std::vector<std::vector<int> >  &,std::vector<int>  &,std::vector<int> &,
			std::vector<std::map<std::string,Eigen::MatrixXd> >  &);
		void LoadEdgeProximal(std::string,std::vector<int>,std::vector<int>&,std::vector<int>&,std::vector<int>&,int);
		void Solve(double,double,double,double);
		void PrintSolution();
		
};

#endif		

