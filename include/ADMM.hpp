#ifndef ADMM_H
#define ADMM_H

#include "Solver/Node.hpp"
#include "Solver/Edge.hpp"
#include "Solver/NodeVar.hpp"
#include "CVXcanon/include/Eigen/Sparse"
#include "CVXcanon/include/Eigen/Core"
#include <unordered_map>
#include <string>
#include <utility>
#include <unordered_map>

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
		ProximalOperator edge_prox;
		ProximalOperator node_prox;
		int prox_edge_arg;
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
		void LoadNodes(std::vector<LinOp* > &,std::vector<std::vector< LinOp *> > &);
		void LoadEdges(std::vector<LinOp* > &,std::vector<std::vector< LinOp *> > &);
		void LoadNodesProximal(ProximalOperator,std::vector<std::vector<int> > &,std::vector<std::vector<std::string> > &,
					std::vector<std::vector<std::vector<int> > > &,std::vector<std::vector<int> > &,std::vector<std::vector<std::vector<double> > > &);
		void LoadEdgesProximal(ProximalOperator,std::vector<std::vector<std::pair<int,int> > > &,std::vector<std::vector<std::pair<int,int> > > &,int);
		void LoadNodeProximal(ProximalOperator, std::vector<int> &,std::vector<std::string>&,std::vector<std::vector<int> >  &,std::vector<int>  &,std::vector<std::map<std::string,Eigen::MatrixXd> >  &);
		void LoadEdgeProximal(ProximalOperator,std::vector<int>,std::vector<int>&,std::vector<int>&,std::vector<int>&,int);
		void Solve(double,double,double,double);
		void PrintSolution();

		//other helper functions
		/*std::vector<double> numpyToVector(double *array,int n){
			std::vector<double> v;
			v.assign(array,array+n);
			return v;
		}*/
		
};

#endif		
//std::unordered_map<int,Solution> Solve(Sense sense, std::unordered_map<int,LinOp* objective> >,std::unordered_map<int,std::vector< LinOp* > constraints >);

