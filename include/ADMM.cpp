#include "ADMM.hpp"
#include <iostream>
#include <thread>

ADMM::ADMM()
{
	node_list.clear();
	edge_list.clear();
	node_x_vals.clear();
	edge_z_vals.clear();
	edge_u_vals.clear();
	edge_prox = NONE;
	node_prox = NONE;
	x_var_size = 0;
	z_var_size = 0;
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
}

void ADMM::LoadNodeProximal(ProximalOperator prox,std::vector<int>  &x_var_idx,
				std::vector<std::string> &x_var_names,
				std::vector<std::vector<int> >  &neighbour_var_idx,
				std::vector<int>  &sizes,
				std::vector<std::vector<double> >  &args)
{
	//need zij and uij for all neighbours of i
	node_prox = prox;

	
	Node *newnode = new Node;
	newnode->node_objective = NULL;
	newnode->node_constraints = std::vector<LinOp *>();
	newnode->neighbour_var_idx = neighbour_var_idx;
	newnode->x_var_idx = x_var_idx;
	newnode->args = std::vector<Eigen::MatrixXd>();
	for ( int j = 0 ; j < x_var_idx.size(); ++j){
		node_x_vals[x_var_idx[j]] = {Eigen::MatrixXd::Constant(sizes[j],1,0),x_var_names[j],0};
		size_x += sizes[j];
		for ( int k = 0 ; k < neighbour_var_idx[j].size(); ++k){
			edge_u_vals[neighbour_var_idx[j][k]] = Eigen::MatrixXd::Constant(sizes[j],1,0);
			edge_z_vals[neighbour_var_idx[j][k]] = Eigen::MatrixXd::Constant(sizes[j],1,0);
			size_z += sizes[j];
		}
		Eigen::MatrixXd argmat = Eigen::MatrixXd(sizes[j],1);
		for ( int k = 0 ; k < args[j].size(); ++k ){
			argmat(k,0) = args[j][k];
		}
		newnode->args.push_back(argmat);
	}
	node_list.push_back(newnode);
}

//for bulk loading of nodes
void ADMM::LoadNodesProximal(ProximalOperator prox,std::vector<std::vector<int> > &x_var_idx,
				std::vector<std::vector<std::string> > &x_var_names,
				std::vector<std::vector<std::vector<int> > > &neighbour_var_idx,
				std::vector<std::vector<int> > &sizes,
				std::vector<std::vector<std::vector<double> > > &args)
{
	//need zij and uij for all neighbours of i
	node_prox = prox;

	//calculate the size required for the solution arrays
	for ( int i = 0 ; i < x_var_idx.size(); ++i ){
		x_var_size += x_var_idx[i].size();
		for ( int j = 0 ; j < x_var_idx[i].size(); ++j ){
			z_var_size += neighbour_var_idx[i][j].size();
		}
	}
	//declare the arrays
	//node_x_vals.resize(x_var_size);// = new Node_Var[x_var_size];
	//edge_z_vals.resize(z_var_size);// = new Eigen::MatrixXd[z_var_size];
	//edge_u_vals.resize(z_var_size);// = new Eigen::MatrixXd[z_var_size];
	
	for ( int i = 0 ; i < x_var_idx.size(); ++i){
		Node *newnode = new Node;
		newnode->node_objective = NULL;
		newnode->node_constraints = std::vector<LinOp *>();
		newnode->neighbour_var_idx = neighbour_var_idx[i];
		newnode->x_var_idx = x_var_idx[i];
		newnode->args = std::vector<Eigen::MatrixXd>();
		for ( int j = 0 ; j < x_var_idx[i].size(); ++j){
			node_x_vals[x_var_idx[i][j]] = {Eigen::MatrixXd::Constant(sizes[i][j],1,0),x_var_names[i][j],i};
			size_x += sizes[i][j];
			for ( int k = 0 ; k < neighbour_var_idx[i][j].size(); ++k){
				edge_u_vals[neighbour_var_idx[i][j][k]] = Eigen::MatrixXd::Constant(sizes[i][j],1,0);
				edge_z_vals[neighbour_var_idx[i][j][k]] = Eigen::MatrixXd::Constant(sizes[i][j],1,0);
				size_z += sizes[i][j];
			}
			Eigen::MatrixXd argmat = Eigen::MatrixXd(sizes[i][j],1);
			for ( int k = 0 ; k < args[i][j].size(); ++k ){
				argmat(k,0) = args[i][j][k];
			}
			newnode->args.push_back(argmat);
		}
		node_list.push_back(newnode);
	}
	/*primal_res = Eigen::MatrixXd::Constant(z_var_size,1,0);
	dual_res = Eigen::MatrixXd::Constant(z_var_size,1,0);
	xnorm = Eigen::MatrixXd::Constant(z_var_size,1,0);
	unorm = Eigen::MatrixXd::Constant(z_var_size,1,0);
	znorm =  Eigen::MatrixXd::Constant(z_var_size,1,0);*/
}

void ADMM::LoadEdgeProximal(ProximalOperator prox,std::vector<std::pair<int,int> > &edge_var_idx,std::vector<std::pair<int,int> > &node_var_idx,int arg=0)
{
	//need xi and uij for all zij
	edge_prox = prox;
	prox_edge_arg = arg;
	Edge *newedge = new Edge;
	newedge->edge_objective = NULL;
	newedge->edge_constraints = std::vector<LinOp *>();
	newedge->edge_var_idx = edge_var_idx;
	newedge->node_var_idx = node_var_idx;
	edge_list.push_back(newedge);
}

//for bulk loading of edges
void ADMM::LoadEdgesProximal(ProximalOperator prox,std::vector<std::vector<std::pair<int,int> > > &edge_var_idx,std::vector<std::vector<std::pair<int,int> > > &node_var_idx,int arg=0)
{
	//need xi and uij for all zij
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
	for(int iter = 0 ; iter <= 1000; ++iter){
		double primal_res = 0;
		double dual_res = 0;
		double xnorm = 0;
		double unorm = 0;
		double znorm = 0;
		std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
		//#if defined(_OPENMP)
			#pragma omp parallel for schedule(dynamic, 32)
		//#endif
		for ( int i = 0 ; i < node_list.size(); ++i ){
			ADMM_node(node_list[i]);
		}
		std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
		std::cout << "Time difference in x = " << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() <<std::endl;
		//#if defined(_OPENMP)
			#pragma omp parallel for schedule(dynamic, 64) reduction(+:primal_res,dual_res,xnorm,unorm,znorm)
		//#endif
		for ( int i = 0 ; i < edge_list.size(); ++i ){
			ADMM_edge(edge_list[i],primal_res,dual_res,xnorm,unorm,znorm);
		}
		std::chrono::steady_clock::time_point t3 = std::chrono::steady_clock::now();
		std::cout << "Time difference in z and u = " << std::chrono::duration_cast<std::chrono::microseconds>(t3 - t2).count() <<std::endl;
		double e_pri = sqrt(size_x) * 0.01 + 0.01 * sqrt(std::max(xnorm,znorm)) + 0.0001;
		double e_dual = sqrt(size_z) * 0.01 + 0.01 * sqrt(unorm) + 0.0001;
		std::cout << sqrt(primal_res) << " " << e_pri << " " << sqrt(dual_res) << " " << e_dual << "\n";
		if ( sqrt(primal_res) <= e_pri && sqrt(dual_res) <= e_dual ){
			break;
		}
	}
}

void ADMM::PrintSolution()
{
	for(int i = 0 ; i < x_var_size; ++ i){
		std::cout <<"Node ID " << node_x_vals[i].nodeId << "\n" << node_x_vals[i].name << " " << node_x_vals[i].value.transpose() << "\n";
	}
}
	
void ADMM::ADMM_node(Node *node){
	//TODO : change the rho of the objective
	//TODO : initialise the increment variable,add parameters in the square objective
	//*node_x_vals[i] = solve(MINIMIZE,node_objectives[i],node_constraints[i],solver_options);
	if ( node->node_objective != NULL ){
		Solution soln = solve(MINIMIZE,node->node_objective,node->node_constraints,solver_options);
	}
	else{
		//solve using prox operator update
		switch(node_prox){
			case SQUARE: 	
					for ( int j = 0 ; j < node->x_var_idx.size(); ++j ){
						Eigen::MatrixXd increment(2*node->args[j]);
						for ( int k = 0 ; k < node->neighbour_var_idx[j].size(); ++k ){
							increment += edge_z_vals[node->neighbour_var_idx[j][k]] - edge_u_vals[node->neighbour_var_idx[j][k]];
						}
						node_x_vals[node->x_var_idx[j]].value = increment*1.0/(2+node->neighbour_var_idx[j].size());
						//std::cout<<"x " << node_list[i]->x_var_idx[j] <<" "<<node_x_vals[node_list[i]->x_var_idx[j]].value<<"\n";
					}
					
					break;
			case LASSO:
					break;
			default:std::cout<<"not implemented yet";exit(-1);
		}
			
	}
}

void ADMM::ADMM_edge(Edge *edge,double &primal_res,double &dual_res,double &xnorm,double &unorm,double &znorm){
	//TODO : change the rho of the objective
	if ( edge->edge_objective != NULL ){
		Solution soln = solve(MINIMIZE,edge->edge_objective,edge->edge_constraints,solver_options);
		//yank out variables,update the eigen matrices
	}
	else{
		switch(edge_prox){
			case SQUARE: 	std::cout << "yada yada\n";
					break;
			case LASSO:	//std::cout << "yada yada z\n";
					for ( int j = 0 ; j < edge->edge_var_idx.size(); ++j ){
						if ( edge->edge_var_idx[j].first != 0 || edge->edge_var_idx[j].second != 0 ){
							Eigen::MatrixXd z_ij = edge_z_vals[edge->edge_var_idx[j].first];
							Eigen::MatrixXd z_ji = edge_z_vals[edge->edge_var_idx[j].second];
							Eigen::MatrixXd u_ij = edge_u_vals[edge->edge_var_idx[j].first];
							Eigen::MatrixXd u_ji = edge_u_vals[edge->edge_var_idx[j].second];
							Eigen::MatrixXd x_i = node_x_vals[edge->node_var_idx[j].first].value;
							Eigen::MatrixXd x_j = node_x_vals[edge->node_var_idx[j].second].value;
							double theta = std::max(1-1.0/(x_i+u_ij - x_j - u_ji).norm(),0.5);
							//std::cout << " Old values " << edge_list[i]->edge_var_idx[j].first << " " << old_z_val_first << " " << edge_list[i]->edge_var_idx[j].second << " " << old_z_val_second << "\n";
							Eigen::MatrixXd sum_i = (x_i + u_ij),sum_j = (x_j + u_ji);
							edge_z_vals[edge->edge_var_idx[j].first] = theta * sum_i + (1-theta) * sum_j;
							edge_z_vals[edge->edge_var_idx[j].second] = ( 1- theta ) * sum_i + theta * sum_j;
							edge_u_vals[edge->edge_var_idx[j].first] += x_i - edge_z_vals[edge->edge_var_idx[j].first];
							edge_u_vals[edge->edge_var_idx[j].second] += x_j - edge_z_vals[edge->edge_var_idx[j].second];
							
							//#pragma omp critical
							{	
								primal_res += (x_i - z_ij).squaredNorm() + (x_j - z_ji).squaredNorm();
								dual_res += (edge_z_vals[edge->edge_var_idx[j].first] - z_ij ).squaredNorm()+
										(edge_z_vals[edge->edge_var_idx[j].second] - z_ji).squaredNorm();
								znorm +=  edge_z_vals[edge->edge_var_idx[j].first].squaredNorm() + edge_z_vals[edge->edge_var_idx[j].second].squaredNorm();
								xnorm += x_i.squaredNorm() +  x_j.squaredNorm();	
								unorm += edge_u_vals[edge->edge_var_idx[j].first].squaredNorm() + edge_u_vals[edge->edge_var_idx[j].second].squaredNorm();
							}
						}
					}
					break;
			default:	std::cout<<"not implemented yet\n";
					exit(-1);
		}
	}
}


	
