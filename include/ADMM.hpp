#ifndef ADMM_H
#define ADMM_H

#include "CVXcanon/src/CVXcanon.hpp"
#include <unordered_map>
#include <string>
#include <utility>
#include <unordered_map>

typedef enum ProximalOperator
{
	SQUARE,
	LASSO,
	NONE
}ProximalOperator;

typedef struct Node 
{
	LinOp* node_objective;
	std::vector<LinOp *> node_constraints;
	std::vector<std::vector<int> > neighbour_var_idx;
	std::vector<int> x_var_idx;
	std::vector<Eigen::MatrixXd> args;
}Node;

typedef struct Edge
{
	LinOp* edge_objective;
	std::vector<LinOp *> edge_constraints;
	std::vector<std::vector<int> > edge_var_idx;
	std::vector<std::vector<int> > node_var_idx;
}Edge;

typedef struct Node_Var
{
	Eigen::MatrixXd value;
	std::string name;
	int nodeId;
}Node_Var;

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
		/*std::vector<Node_Var> node_x_vals;
		std::vector<Eigen::MatrixXd> edge_z_vals;
		std::vector<Eigen::MatrixXd> edge_u_vals;*/
		ProximalOperator edge_prox;
		ProximalOperator node_prox;
		int prox_edge_arg;
		int x_var_size;
		int z_var_size;

		std::map<std::string,double > solver_options;

		//for convergence
		/*Eigen::MatrixXd primal_res;
		Eigen::MatrixXd dual_res;
		Eigen::MatrixXd xnorm;
		Eigen::MatrixXd znorm;
		Eigen::MatrixXd unorm;*/
		double primal_res;
		double dual_res;
		double xnorm;
		double znorm;
		double unorm;
		double size_x;
		double size_z;

		//solver functions
		void ADMM_node(Node *);
		void ADMM_edge(Edge *);
		//void copyLinop(LinOp* const&,LinOp* const&);
	public:
		ADMM();
		/*void LoadNodes(std::vector<LinOp* > &,std::vector<std::vector< LinOp *> > &);
		void LoadEdges(std::vector<LinOp* > &,std::vector<std::vector< LinOp *> > &);*/
		void LoadNodes(std::vector<LinOp* > &,std::vector<std::vector< LinOp *> > &);
		void LoadEdges(std::vector<LinOp* > &,std::vector<std::vector< LinOp *> > &);
		void LoadNodesProximal(ProximalOperator,std::vector<std::vector<int> > &,std::vector<std::vector<std::string> > &,
					std::vector<std::vector<std::vector<int> > > &,std::vector<std::vector<int> > &,std::vector<std::vector<std::vector<double> > > &);
		void LoadEdgesProximal(ProximalOperator,std::vector<std::vector<std::vector<int> > > &,std::vector<std::vector<std::vector<int> > > &,int);
		void LoadNodeProximal(ProximalOperator, std::vector<int> &,std::vector<std::string>&,std::vector<std::vector<int> >  &,std::vector<int>  &,std::vector<std::vector<double> >  &);
		void LoadEdgeProximal(ProximalOperator,std::vector<std::vector<int> >  &,std::vector<std::vector<int> >  &,int);
		void Solve();
		void PrintSolution();
		/*std::vector<std::unordered_map<int, Eigen::MatrixXd> > get_node_x_vals();
		std::vector<std::unordered_map<int, Eigen::MatrixXd> > get_edge_z_vals();
		std::vector<std::unordered_map<int, Eigen::MatrixXd> > get_edge_u_vals();*/
};

#endif		
//std::unordered_map<int,Solution> Solve(Sense sense, std::unordered_map<int,LinOp* objective> >,std::unordered_map<int,std::vector< LinOp* > constraints >);

