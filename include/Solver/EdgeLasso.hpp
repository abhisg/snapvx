#ifndef EDGELASSO_H_
#define EDGELASSO_H_
#include "Edge.hpp"

class EdgeLasso : public Edge{
protected:
	double lambda;
public:
	EdgeLasso():Edge(){lambda=1.0;}

	void SetLambda(double lambda){
		this->lambda = lambda;
	}

	virtual void ADMM_edge(std::unordered_map<int,Node_Var> &node_x_vals,
				std::unordered_map<int,Eigen::MatrixXd> &edge_z_vals,
				std::unordered_map<int,Eigen::MatrixXd>	&edge_u_vals,
				double &rho){
		//std::vector<double> norms(5,0);
		//SaveState(node_x_vals,edge_z_vals,edge_u_vals);
		int p = int p = (int)((-1+sqrt(1+4*this->x_var_idx.size()))/2.0);
		for ( int j = 0; j < p; ++j){
			edge_z_vals[this->edge_var_idx_left[j]] = node_x_vals[this->node_var_idx_left[j]].value + 
														edge_u_vals[this->edge_var_idx_left[j]];
		}
		for ( int j = p; j < this->edge_var_idx_left.size(); ++j){
			if ( j % p == 0 ){
				edge_z_vals[this->edge_var_idx_left[j]] = node_x_vals[this->node_var_idx_left[j]].value + 
														edge_u_vals[this->edge_var_idx_left[j]];
			}
			else{
				for ( int k = 0; k < edge_z_vals[this->edge_var_idx_left[j]].cols(); ++k){
					Eigen::MatrixXd sum = edge_u_vals[this->edge_var_idx_left[j]].col(k) + 
										node_x_vals[this->node_var_idx_left[j]].value.col(k);
					if ( sum.norm() > 1.0/lambda ){
						edge_z_vals[this->edge_var_idx_left[j]].col(k) = sum * ( 1- 1.0/(lambda * sum.norm()));
					}
					else{
						edge_z_vals[this->edge_var_idx_left[j]].col(k) = 0;
					}
				}
			}
		}
		for ( int j = 0; j < this->edge_var_idx_left.size(); ++j){
			edge_u_vals[this->edge_var_idx_left[j]]] += edge_z_vals[this->edge_var_idx_left[j]]] +
														node_x_vals[this->node_var_idx_left[j]].value;
		}
	}
};
#endif
