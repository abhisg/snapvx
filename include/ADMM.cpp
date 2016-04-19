#include "ADMM.h"
#include <iostream>

ADMM::ADMM()
{
	node_objectives.clear();
	node_constraints.clear();
	edge_objectives.clear();
	edge_constraints.clear();
}
	
//TODO : improve the shitty memory management
void ADMM::LoadNodes(std::vector<LinOp* > &objectives,std::vector<std::vector< LinOp *> > &constraints)
{
	for(int i = 0 ; i < objectives.size(); ++i ){
		node_objectives.push_back(objectives[i]);
	}
	for ( int i = 0 ; i < constraints.size(); ++i ){
		std::vector<LinOp*> cons_vec;
		std::cout << "Constraint " << i << " " <<  constraints[i].size() << "\n";
		for ( int j = 0 ; j < constraints[i].size(); j++ ){
			cons_vec.push_back(constraints[i][j]);
		}
		node_constraints.push_back(cons_vec);
	}
}

void ADMM::LoadEdges(std::vector<LinOp* > &objectives,std::vector<std::vector< LinOp *> > &constraints)
{
	for(int i = 0 ; i < objectives.size(); ++i ){
		edge_objectives.push_back(objectives[i]);
	}
	for ( int i = 0 ; i < constraints.size(); ++i ){
		std::vector<LinOp*> cons_vec;
		std::cout << "Constraint " << i << " " <<  constraints[i].size() << "\n";
		for ( int j = 0 ; j < constraints[i].size(); j++ ){
			cons_vec.push_back(constraints[i][j]);
		}
		edge_constraints.push_back(cons_vec);
	}
}


/*void ADMM::LoadEdges(std::vector<LinOp* > &objectives,std::vector<std::vector< LinOp *> > &constraints)
{
	for(int i = 0 ; i < objectives.size(); ++i ){
		//Solution *usol = new Solution;
		//Solution *zsol = new Solution;
		edge_objectives.push_back(objectives[i]);
		//edge_z_vals.push_back(zsol);
		//edge_u_vals.push_back(usol);
	}
	for ( int i = 0 ; i < constraints.size(); ++i ){
		/*std::vector<LinOp *> newcons;
		for ( int j = 0 ; j < constraints[i].size(); ++j){
			newcons.push_back(constraints[i][j]);
		}
		edge_constraints.push_back(newcons);
		edge_constraints.push_back(constraints[i]);
	}
}*/

void ADMM::Solve()
{
	//TODO : full implementation with convergence
	for(int iter = 0 ; iter < 1; ++iter){
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
	//*node_x_vals[i] = solve(MINIMIZE,node_objectives[i],node_constraints[i],solver_options);
	Solution soln = solve(MINIMIZE,node_objectives[i],node_constraints[i],solver_options);
}

void ADMM::ADMM_z(int i){
	//TODO : change the rho of the objective
	//*edge_z_vals[i] = solve(MINIMIZE,edge_objectives[i],edge_constraints[i],solver_options);
	Solution soln = solve(MINIMIZE,edge_objectives[i],edge_constraints[i],solver_options);
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
	
