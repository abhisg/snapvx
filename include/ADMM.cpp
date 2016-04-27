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

void ADMM::LoadNodesProximal(ProximalOperator prox,std::vector<std::vector<int> > &x_var_idx,
				std::vector<std::vector<std::string> > &x_var_names,
				std::vector<std::vector<std::vector<int> > > &neighbour_var_idx,
				std::vector<std::vector<int> > &sizes,
				std::vector<std::vector<std::vector<double> > > &args)
{
	//need zij and uij for all neighbours of i
	node_prox = prox;
	for ( int i = 0 ; i < x_var_idx.size(); ++i){
		Node *newnode = new Node;
		newnode->node_objective = NULL;
		newnode->node_constraints = std::vector<LinOp *>();
		newnode->neighbour_var_idx = neighbour_var_idx[i];
		newnode->x_var_idx = x_var_idx[i];
		newnode->args = std::vector<Eigen::MatrixXd>();
		for ( int j = 0 ; j < x_var_idx[i].size(); ++j){
			node_x_vals[x_var_idx[i][j]] = {Eigen::MatrixXd(sizes[i][j],1),x_var_names[i][j],i};
			size_x += sizes[i][j];
			for ( int k = 0 ; k < neighbour_var_idx[i][j].size(); ++k){
				edge_u_vals[neighbour_var_idx[i][j][k]] = Eigen::MatrixXd(sizes[i][j],1);
				edge_z_vals[neighbour_var_idx[i][j][k]] = Eigen::MatrixXd(sizes[i][j],1);
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
	primal_res = Eigen::MatrixXd(edge_z_vals.size(),1);
	dual_res = Eigen::MatrixXd(edge_z_vals.size(),1);
	xnorm = Eigen::MatrixXd(edge_z_vals.size(),1);
	unorm = Eigen::MatrixXd(edge_z_vals.size(),1);
	znorm =  Eigen::MatrixXd(edge_z_vals.size(),1);
}

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
	//TODO : full implementation with convergence
	for(int iter = 0 ; iter <= 1000; ++iter){
		double e_pri = sqrt(size_x) * 0.01 + 0.01 * sqrt(std::max(xnorm.sum(),znorm.sum())) + 0.0001;
		double e_dual = sqrt(size_z) * 0.01 + 0.01 * sqrt(unorm.sum()) + 0.0001;
		for ( int i = 0 ; i < node_list.size(); ++i ){
			ADMM_x(i);
		}
		for ( int i = 0 ; i < edge_list.size(); ++i ){
			ADMM_z(i);
		}
		for ( int i = 0 ; i < edge_list.size(); ++i ){
			ADMM_u(i);
		}
		std::cout << sqrt(primal_res.sum()) << " " << e_pri << " " << sqrt(dual_res.sum()) << " " << e_dual << "\n";
		if ( sqrt(primal_res.sum()) <= e_pri && sqrt(dual_res.sum()) <= e_dual ){
			break;
		}
	}
}

void ADMM::PrintSolution()
{
	for(int i = 0 ; i < node_x_vals.size(); ++ i){
		std::cout <<"Node ID " << node_x_vals[i].nodeId << "\n" << node_x_vals[i].name << " " << node_x_vals[i].value.transpose() << "\n";
	}
}
	
void ADMM::ADMM_x(int i){
	//TODO : change the rho of the objective
	//TODO : initialise the increment variable,add parameters in the square objective
	//*node_x_vals[i] = solve(MINIMIZE,node_objectives[i],node_constraints[i],solver_options);
	if ( node_list[i]->node_objective != NULL ){
		Solution soln = solve(MINIMIZE,node_list[i]->node_objective,node_list[i]->node_constraints,solver_options);
	}
	else{
		//solve using prox operator update
		switch(node_prox){
			case SQUARE: 	
					for ( int j = 0 ; j < node_list[i]->x_var_idx.size(); ++j ){
						Eigen::MatrixXd increment(2*node_list[i]->args[j]);
						for ( int k = 0 ; k < node_list[i]->neighbour_var_idx[j].size(); ++k ){
							increment += edge_z_vals[node_list[i]->neighbour_var_idx[j][k]] - edge_u_vals[node_list[i]->neighbour_var_idx[j][k]];
						}
						node_x_vals[node_list[i]->x_var_idx[j]].value = increment*1.0/(2+node_list[i]->neighbour_var_idx[j].size());
						//std::cout<<"x " << node_list[i]->x_var_idx[j]<<" "<<node_x_vals[node_list[i]->x_var_idx[j]]<<"\n";
					}
					
					break;
			case LASSO:
					break;
			default:std::cout<<"not implemented yet";exit(-1);
		}
			
	}
}

void ADMM::ADMM_z(int i){
	//TODO : change the rho of the objective
	//*edge_z_vals[i] = solve(MINIMIZE,edge_objectives[i],edge_constraints[i],solver_options);
	if ( edge_list[i]->edge_objective != NULL ){
		Solution soln = solve(MINIMIZE,edge_list[i]->edge_objective,edge_list[i]->edge_constraints,solver_options);
		//yank out variables,update the eigen matrices
	}
	else{
		switch(edge_prox){
			case SQUARE: 	std::cout << "yada yada\n";
					break;
			case LASSO:	//std::cout << "yada yada z\n";
					for ( int j = 0 ; j < edge_list[i]->edge_var_idx.size(); ++j ){
						if ( edge_list[i]->edge_var_idx[j].first != 0 || edge_list[i]->edge_var_idx[j].second != 0 ){
							double theta = std::max(1-1.0/(node_x_vals[edge_list[i]->node_var_idx[j].first].value+
										edge_u_vals[edge_list[i]->edge_var_idx[j].first]-
										node_x_vals[edge_list[i]->node_var_idx[j].second].value-
										edge_u_vals[edge_list[i]->edge_var_idx[j].second]).norm(),0.5);
							Eigen::MatrixXd old_z_val_first = edge_z_vals[edge_list[i]->edge_var_idx[j].first];
							Eigen::MatrixXd old_z_val_second = edge_z_vals[edge_list[i]->edge_var_idx[j].second];
							edge_z_vals[edge_list[i]->edge_var_idx[j].first] = theta * (node_x_vals[edge_list[i]->node_var_idx[j].first].value + 
														edge_u_vals[edge_list[i]->edge_var_idx[j].first])+
												(1-theta) * (node_x_vals[edge_list[i]->node_var_idx[j].second].value+
										edge_u_vals[edge_list[i]->edge_var_idx[j].second]);
							edge_z_vals[edge_list[i]->edge_var_idx[j].second] = theta * (node_x_vals[edge_list[i]->node_var_idx[j].second].value + 
														edge_u_vals[edge_list[i]->edge_var_idx[j].second])+
												(1-theta) * (node_x_vals[edge_list[i]->node_var_idx[j].first].value+
										edge_u_vals[edge_list[i]->edge_var_idx[j].first]);
							primal_res(edge_list[i]->edge_var_idx[j].first,0) = (node_x_vals[edge_list[i]->node_var_idx[j].first].value - \
													edge_z_vals[edge_list[i]->edge_var_idx[j].first]).squaredNorm();
							primal_res(edge_list[i]->edge_var_idx[j].second,0) = (node_x_vals[edge_list[i]->node_var_idx[j].second].value - \
													edge_z_vals[edge_list[i]->edge_var_idx[j].second]).squaredNorm();
							dual_res(edge_list[i]->edge_var_idx[j].first,0) = (edge_z_vals[edge_list[i]->edge_var_idx[j].first] - old_z_val_first ).squaredNorm();
							dual_res(edge_list[i]->edge_var_idx[j].second,0) = (edge_z_vals[edge_list[i]->edge_var_idx[j].second] - old_z_val_second).squaredNorm();
							znorm(edge_list[i]->edge_var_idx[j].first,0) = edge_z_vals[edge_list[i]->edge_var_idx[j].first].squaredNorm();
							znorm(edge_list[i]->edge_var_idx[j].second,0) = edge_z_vals[edge_list[i]->edge_var_idx[j].second].squaredNorm();
							xnorm(edge_list[i]->edge_var_idx[j].first,0) = node_x_vals[edge_list[i]->node_var_idx[j].first].value.squaredNorm();
							xnorm(edge_list[i]->edge_var_idx[j].second,0) = node_x_vals[edge_list[i]->node_var_idx[j].second].value.squaredNorm();
							
						}
					}
					break;
			default:	std::cout<<"not implemented yet\n";
					exit(-1);
		}
	}
}

void ADMM::ADMM_u(int i){
	for ( int j = 0 ; j < edge_list[i]->edge_var_idx.size(); ++j ){
		if ( edge_list[i]->edge_var_idx[j].first != 0 || edge_list[i]->edge_var_idx[j].second != 0 ){
			edge_u_vals[edge_list[i]->edge_var_idx[j].first] += node_x_vals[edge_list[i]->node_var_idx[j].first].value - edge_z_vals[edge_list[i]->edge_var_idx[j].first];
			edge_u_vals[edge_list[i]->edge_var_idx[j].second] += node_x_vals[edge_list[i]->node_var_idx[j].second].value - edge_z_vals[edge_list[i]->edge_var_idx[j].second];
			unorm(edge_list[i]->edge_var_idx[j].first,0) = edge_u_vals[edge_list[i]->edge_var_idx[j].first].squaredNorm();
			unorm(edge_list[i]->edge_var_idx[j].second,0) = edge_u_vals[edge_list[i]->edge_var_idx[j].second].squaredNorm();
		}

	}
		
}

	
