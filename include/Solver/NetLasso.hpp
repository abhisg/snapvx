#ifndef NETLASSO_H_
#define NETLASSO_H_
#include "Edge.hpp"

class NetLasso : public Edge{
protected:
	double lambda;
public:
	NetLasso():Edge(){lambda=1.0;}

	void SetLambda(double lambda){
		this->lambda = lambda;
	}

	virtual std::vector<double> ADMM_edge(std::unordered_map<int,Node_Var> &node_x_vals,
				std::unordered_map<int,Eigen::MatrixXd> &edge_z_vals,
				std::unordered_map<int,Eigen::MatrixXd>	&edge_u_vals,
				double &rho){
		//std::vector<double> norms(5,0);
		SaveState(node_x_vals,edge_z_vals,edge_u_vals);
		for ( int j = 0 ; j < this->edge_var_idx_left.size(); ++j ){
			if ( this->edge_var_idx_left[j] != 0 || this->edge_var_idx_right[j] != 0 ){
				Eigen::MatrixXd z_ij = edge_z_vals[this->edge_var_idx_left[j]];
				Eigen::MatrixXd z_ji = edge_z_vals[this->edge_var_idx_right[j]];
				Eigen::MatrixXd u_ij = edge_u_vals[this->edge_var_idx_left[j]];
				Eigen::MatrixXd u_ji = edge_u_vals[this->edge_var_idx_right[j]];
				Eigen::MatrixXd x_i = node_x_vals[this->node_var_idx_left[j]].value;
				Eigen::MatrixXd x_j = node_x_vals[this->node_var_idx_right[j]].value;

				//actual solver
				double theta = std::max(1-lambda/(rho*(x_i+u_ij - x_j - u_ji).norm()),0.5);
				Eigen::MatrixXd sum_i = (x_i + u_ij),sum_j = (x_j + u_ji);
				edge_z_vals[this->edge_var_idx_left[j]] = theta * sum_i + (1-theta) * sum_j;
				edge_z_vals[this->edge_var_idx_right[j]] = ( 1- theta ) * sum_i + theta * sum_j;
				edge_u_vals[this->edge_var_idx_left[j]] += x_i - edge_z_vals[this->edge_var_idx_left[j]];
				edge_u_vals[this->edge_var_idx_right[j]] += x_j - edge_z_vals[this->edge_var_idx_right[j]];

				/*norms[0] += (x_i - z_ij).squaredNorm() + (x_j - z_ji).squaredNorm();
				norms[1] += ((edge_z_vals[this->edge_var_idx_left[j]] - z_ij ).squaredNorm()+
						(edge_z_vals[this->edge_var_idx_right[j]] - z_ji).squaredNorm());
				norms[2] +=  edge_z_vals[this->edge_var_idx_left[j]].squaredNorm() + edge_z_vals[this->edge_var_idx_right[j]].squaredNorm();
				norms[3] += x_i.squaredNorm() +  x_j.squaredNorm();	
				norms[4] += (edge_u_vals[this->edge_var_idx_left[j]].squaredNorm() + 
					edge_u_vals[this->edge_var_idx_right[j]].squaredNorm());*/
			}
		}
		return CalculateNorms(node_x_vals,edge_z_vals,edge_u_vals);
	}
};
#endif
