#include "ADMM.hpp"
#include <iostream>

ADMM::ADMM()
{
	node_list.clear();
	edge_list.clear();
	node_x_vals.clear();
	edge_z_vals.clear();
	edge_u_vals.clear();
	x_var_size = 0;
	z_var_size = 0;
	solver_options = std::map<std::string,double>();
	size_x = 0;
	size_z = 0;
	rho = 1;
	lambda = 1;
	eps_abs = 0.001;
	eps_rel = 0.001;
	ProximalMap::LoadOperators();
}
	
//TODO : improve the shitty memory management
/*void ADMM::LoadNodes(std::vector<LinOp* > &objectives,std::vector<std::vector< LinOp *> > &constraints)
{
	for ( int i = 0 ; i < objectives.size(); ++i){
		Node *newnode = new Node;
		newnode->node_objective = objectives[i];
		newnode->node_constraints = constraints[i];	//sloppy but we do not reclaim the memory before the end of the python function
		newnode->neighbour_var_idx = std::vector<std::vector<int> >();
		newnode->x_var_idx = std::vector<int>();
		node_list.push_back(newnode);
	}	
}

void ADMM::LoadEdges(std::vector<LinOp* > &objectives,std::vector<std::vector< LinOp *> > &constraints)
{
	for ( int i = 0 ; i < objectives.size(); ++i){
		Edge *newedge = new Edge;
		newedge->edge_objective = objectives[i];
		newedge->edge_constraints = constraints[i];	//sloppy but we do not reclaim the memory before the end of the python function
		edge_list.push_back(newedge);
	}	
}*/

void ADMM::LoadNodeProximal(std::string prox,std::vector<int>  &x_var_idx,
				std::vector<std::string> &x_var_names,
				std::vector<std::vector<int> >  &neighbour_var_idx,
				std::vector<int>  &sizes_i,
				std::vector<int>  &sizes_j,
				std::vector<std::map<std::string,Eigen::MatrixXd > >  &args)
{
	//need zij and uij for all neighbours of i
	x_var_size += x_var_idx.size();
	Node *newnode = ProximalMap::getNodeInstance(prox);
	if ( newnode != NULL){
		newnode->node_objective = NULL;
		newnode->node_constraints = std::vector<LinOp *>();
		newnode->neighbour_var_idx = neighbour_var_idx;
		newnode->x_var_idx = x_var_idx;
		newnode->args = std::vector<std::map<std::string,Eigen::MatrixXd> >();
		for ( int j = 0 ; j < x_var_idx.size(); ++j){
			node_x_vals[x_var_idx[j]] = {Eigen::MatrixXd::Constant(sizes_i[j],sizes_j[j],0),x_var_names[j],0};
			size_x += sizes_i[j]*sizes_j[j];
			size_z += (sizes_i[j]*sizes_j[j])*neighbour_var_idx[j].size();
			z_var_size += neighbour_var_idx[j].size();
			for ( int k = 0 ; k < neighbour_var_idx[j].size(); ++k){
				edge_u_vals[neighbour_var_idx[j][k]] = Eigen::MatrixXd::Constant(sizes_i[j],sizes_j[j],0);
				edge_z_vals[neighbour_var_idx[j][k]] = Eigen::MatrixXd::Constant(sizes_i[j],sizes_j[j],0);
			}
			newnode->args.push_back(args[j]);
		}
	}
	node_list.push_back(newnode);
}

//for bulk loading of nodes

void ADMM::LoadEdgeProximal(std::string prox,std::vector<int> edge_var_idx_left,std::vector<int> &edge_var_idx_right,std::vector<int> &node_var_idx_left,std::vector<int> &node_var_idx_right,int arg=0)
{
	//need xi and uij for all zij
	Edge *newedge = ProximalMap::getEdgeInstance(prox);
	if ( newedge != NULL){
		newedge->edge_objective = NULL;
		newedge->edge_constraints = std::vector<LinOp *>();
		newedge->edge_var_idx_left = edge_var_idx_left;
		newedge->edge_var_idx_right = edge_var_idx_right;
		newedge->node_var_idx_left = node_var_idx_left;
		newedge->node_var_idx_right = node_var_idx_right;
	}
	edge_list.push_back(newedge);
}


void ADMM::Solve(double rho, double eps_abs, double eps_rel, double lambda)
{
	this->rho = rho;
	this->eps_abs = eps_abs;
	this->eps_rel = eps_rel;
	this->lambda = lambda;
	for(int iter = 0 ; iter <= 250; ++iter){
		double primal_res = 0;
		double dual_res = 0;
		double xnorm = 0;
		double unorm = 0;
		double znorm = 0;
		//std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
		#if defined(_OPENMP)
			#pragma omp parallel for schedule(dynamic, 32)
		#endif
		for ( int i = 0 ; i < node_list.size(); ++i ){
			if ( node_list[i] ){
				node_list[i]->ADMM_node(node_x_vals,edge_z_vals,edge_u_vals,rho);
			}
		}
		//std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
		//std::cout << "Time difference in x = " << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() <<std::endl;
		#if defined(_OPENMP)
			#pragma omp parallel for schedule(dynamic, 64) reduction(+:primal_res,dual_res,xnorm,unorm,znorm)
		#endif
		for ( int i = 0 ; i < edge_list.size(); ++i ){
			std::vector<double> norms(5,0);
			if ( edge_list[i] ){
				edge_list[i]->SetLambda(lambda);
			 	norms = edge_list[i]->ADMM_edge_solver(node_x_vals,edge_z_vals,edge_u_vals,rho);
			}
			primal_res += norms[0];
			dual_res += norms[1];
			znorm += norms[2];
			xnorm += norms[3];
			unorm += norms[4];
		}
		//std::chrono::steady_clock::time_point t3 = std::chrono::steady_clock::now();
		//std::cout << "Time difference in z and u = " << std::chrono::duration_cast<std::chrono::microseconds>(t3 - t2).count() <<std::endl;
		double e_pri = sqrt(size_x) * eps_abs + eps_rel * sqrt(std::max(xnorm,znorm)) + 0.0001;
		double e_dual = sqrt(size_z) * eps_abs + eps_rel * sqrt(unorm) * rho + 0.0001;
		std::cout << "Size " << sqrt(size_z) << " e_abs " << eps_abs << " rho " << rho << " " << " unorm " << sqrt(unorm) << "\n";
		std::cout << "r: " << sqrt(primal_res) << " e_pri: " << e_pri << " s: " << rho * sqrt(dual_res) << " e_dual: " << e_dual << "\n";
		if ( sqrt(primal_res) <= e_pri && rho * sqrt(dual_res) <= e_dual ){
			break;
		}
	}
}

void ADMM::PrintSolution()
{
	for(int i = 0 ; i < x_var_size; ++i){
		std::cout <<"Node ID " << node_x_vals[i].nodeId << "\n" << node_x_vals[i].name << " " << node_x_vals[i].value.transpose() << "\n";
	}
}


	
