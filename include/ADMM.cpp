#include "ADMM.h"

void ADMM::LoadNodes(std::vector<LinOp* > &objectives,std::vector<std::vector< LinOp *> > &constraints)
{
	for(int i = 0 ; i < objectives.size(); ++i ){
		LinOp *newobj = new LinOp;
		newobj->type = objectives[i]->type;
		newobj->size = objectives[i]->size;
		newobj->sparse = objectives[i]->sparse;
		newobj->dense_data = objectives[i]->dense_data;
		newobj->slice = objectives[i]->slice;
		node_objectives.push_back(newobj);
		node_x_vals.push_back(new Solution);
	}
	for ( int i = 0 ; i < constraints.size(); ++i ){
		std::vector<LinOp *> newcons;
		for ( int j = 0 ; j < constraints[i].size(); ++j){
			LinOp *cons = new LinOp;
			cons->type = constraints[i][j]->type;
			cons->size = constraints[i][j]->size;
			cons->sparse = constraints[i][j]->sparse;
			cons->dense_data = constraints[i][j]->dense_data;
			cons->slice = constraints[i][j]->slice;
			newcons.push_back(cons);
		}
		node_constraints.push_back(newcons);
	}
}

void ADMM::LoadEdges(std::vector<LinOp* > &objectives,std::vector<std::vector< LinOp *> > &constraints)
{
	for(int i = 0 ; i < objectives.size(); ++i ){
		LinOp *newobj = new LinOp;
		newobj->type = objectives[i]->type;
		newobj->size = objectives[i]->size;
		newobj->sparse = objectives[i]->sparse;
		newobj->dense_data = objectives[i]->dense_data;
		newobj->slice = objectives[i]->slice;
		edge_objectives.push_back(newobj);
		edge_z_vals.push_back(new Solution);
		edge_u_vals.push_back(new Solution);
	}
	for ( int i = 0 ; i < constraints.size(); ++i ){
		std::vector<LinOp *> newcons;
		for ( int j = 0 ; j < constraints[i].size(); ++j){
			LinOp *cons = new LinOp;
			cons->type = constraints[i][j]->type;
			cons->size = constraints[i][j]->size;
			cons->sparse = constraints[i][j]->sparse;
			cons->dense_data = constraints[i][j]->dense_data;
			cons->slice = constraints[i][j]->slice;
			newcons.push_back(cons);
		}
		edge_constraints.push_back(newcons);
	}
}

void ADMM::Solve()
{
	//TODO : full implementation with convergence
	for(int iter = 0 ; iter < 10; ++iter){
		for ( int i = 0 ; i < node_objectives.size(); ++i ){
			ADMM_x(i);
		}
		for ( int i = 0 ; i < edge_objectives.size(); ++i ){
			ADMM_z(i);
			ADMM_u(i);
		}
	}
}

void ADMM::ADMM_x(int i){
	//TODO : change the rho of the objective
	*node_x_vals[i] = solve(MINIMIZE,node_objectives[i],node_constraints[i],solver_options);
}

void ADMM::ADMM_z(int i){
	//TODO : change the rho of the objective
	*edge_z_vals[i] = solve(MINIMIZE,edge_objectives[i],edge_constraints[i],solver_options);
}

void ADMM::ADMM_u(int i){
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
	
