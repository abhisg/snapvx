#include "ADMM.h"
#include <iostream>

ADMM::ADMM()
{
	node_list.clear();
	edge_list.clear();
	node_x_vals.clear();
	edge_z_vals.clear();
	edge_u_vals.clear();
	edge_prox = NONE;
	node_prox = NONE;
	prox_node_arg = 0;
	prox_edge_arg = 0;
	solver_options = std::map<std::string,double>();
}
	
//TODO : improve the shitty memory management
void ADMM::LoadNodes(std::vector<LinOp* > &objectives,std::vector<std::vector< LinOp *> > &constraints)
{
for ( int i = 0 ; i < objectives.size(); ++i){
		Node *newnode = new Node;
		newnode->node_objective = objectives[i];
		newnode->node_constraints = constraints[i];	//sloppy but we do not reclaim the memory before the end of the python function
		newnode->z_var_idx = std::vector<std::vector<int> >();
		newnode->u_var_idx = std::vector<std::vector<int> >();
		newnode->x_var_idx = std::vector<std::vector<int> >();
		node_list.push_back(newnode);
	}	
}

void ADMM::LoadEdges(std::vector<LinOp* > &objectives,std::vector<std::vector< LinOp *> > &constraints)
{
	for ( int i = 0 ; i < objectives.size(); ++i){
		Node *newedge = new Edge;
		newedge->edge_objective = objectives[i];
		newedge->edge_constraints = constraints[i];	//sloppy but we do not reclaim the memory before the end of the python function
		newedge->zij_var_idx = std::vector<int>();
		newedge->zji_var_idx = std::vector<int>();
		newedge->uij_var_idx = std::vector<int>();
		newedge->uji_var_idx = std::vector<int>();
		newedge->xi_var_idx = std::vector<int>();
		newedge->xj_var_idx = std::vector<int>();
		edge_list.push_back(newedge);
	}	
}

void LoadNodesProximal(ProximalOperator prox,std::vector<std::vector<int> > x_var_idx,std::vector<std::vector<std::vector<int> > > neighbour_var_idx,arg=0)
{
	node_prox = prox;
	prox_node_arg = arg;
	for ( int i = 0 ; i < x_var_idx.size(); ++i){
		Node *newnode = new Node;
		newnode->node_objective = NULL;
		newnode->node_constraints = std::vector<LinOp *>();
		newnode->neighbour_var_idx = neighbour_var_idx[i];
		newnode->x_var_idx = x_var_idx[i];
		node_list.push_back(newnode);
	}
}

void LoadEdgesProximal(ProximalOperator prox,std::vector<std::vector<std::pair<int,int> > > edge_var_idx,std::vector<std::vector<std::pair<int,int> > > node_var_idx,arg=0)
{
	edge_prox = prox;
	prox_edge_arg = arg;
	for ( int i = 0 ; i < edge_var_idx.size(); ++i ){
		Edge *newedge = new Edge;
		newedge->edge_objective = NULL;
		newedge->edge_constraints = std::vector<LinOp *>();
		newedge->edge_var_idx = edge_var_idx[i];
		newedge->node_var_idx = node_var_idx[i];
		edge_list.push_back(newedge);
	}
}

void ADMM::Solve()
{
	//TODO : full implementation with convergence
	for(int iter = 0 ; iter < 1; ++iter){
		for ( int i = 0 ; i < node_list.size(); ++i ){
			ADMM_x(i);
		}
		for ( int i = 0 ; i < edge_list.size(); ++i ){
			ADMM_z(i);
			ADMM_u(i);
		}
	}
}

void ADMM::ADMM_x(int i){
	//TODO : change the rho of the objective
	//*node_x_vals[i] = solve(MINIMIZE,node_objectives[i],node_constraints[i],solver_options);
	if ( node_list[i]->node_objective != NULL ){
		Solution soln = solve(MINIMIZE,node_objectives[i],node_constraints[i],solver_options);
	}
	else{
		//solve using prox operator update
	}
}

void ADMM::ADMM_z(int i){
	//TODO : change the rho of the objective
	//*edge_z_vals[i] = solve(MINIMIZE,edge_objectives[i],edge_constraints[i],solver_options);
	if ( edge_list[i]->edge_objective != NULL ){
		Solution soln = solve(MINIMIZE,edge_objectives[i],edge_constraints[i],solver_options);
	}
}

void ADMM::ADMM_u(int i){
	std::cout << "ADMM u called\n";
	//TODO
}

std::vector<std::map<int, Eigen::MatrixXd> > ADMM::get_node_x_vals()
{
	std::vector<std::map<int, Eigen::MatrixXd> > node_x;
	for ( int i = 0 ; i < node_x_vals.size(); ++i ){
		node_x.push_back(node_x_vals[i]->primal_values);
	}
	return node_x;
}

std::vector<std::map<int, Eigen::MatrixXd> > ADMM::get_edge_z_vals()
{
	std::vector<std::map<int, Eigen::MatrixXd> > edge_z;
	for ( int i = 0 ; i < edge_z_vals.size(); ++i ){
		edge_z.push_back(edge_z_vals[i]->primal_values);
	}
	return edge_z;
}

std::vector<std::map<int, Eigen::MatrixXd> > ADMM::get_edge_u_vals()
{
	std::vector<std::map<int, Eigen::MatrixXd> > edge_u;
	for ( int i = 0 ; i < edge_u_vals.size(); ++i ){
		edge_u.push_back(edge_u_vals[i]->primal_values);
	}
	return edge_u;
}
	
