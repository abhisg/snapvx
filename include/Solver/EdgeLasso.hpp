#ifndef EDGELASSO_H_
#define EDGELASSO_H_
#include "Edge.hpp"

class EdgeLasso : public Edge{
public:
	EdgeLasso():Edge(){lambda=1.0;}

	virtual void ADMM_edge(std::unordered_map<int,Node_Var> &node_x_vals,
				std::unordered_map<int,Eigen::MatrixXd> &edge_z_vals,
				std::unordered_map<int,Eigen::MatrixXd>	&edge_u_vals,
				double &rho){
		//std::vector<double> norms(5,0);
		//SaveState(node_x_vals,edge_z_vals,edge_u_vals);
		int p = (int)((-1+sqrt(1+4*(this->edge_var_idx_left.size()-1)))/2.0);
		edge_z_vals[this->edge_var_idx_left[0]] = node_x_vals[this->node_var_idx_left[0]].value + 
														edge_u_vals[this->edge_var_idx_left[0]];
		std::cout << "z " << this->edge_var_idx_left[0] << " " << edge_z_vals[this->edge_var_idx_left[0]] << "\n";
		for ( int j = 1; j < p + 1; ++j){
			//edge_z_vals[this->edge_var_idx_left[j]] = node_x_vals[this->node_var_idx_left[j]].value + 
			//											edge_u_vals[this->edge_var_idx_left[j]];
			//std::cout << "z " << this->edge_var_idx_left[j] << " " << edge_z_vals[this->edge_var_idx_left[j]] << "\n";
			Eigen::MatrixXd sum = edge_u_vals[this->edge_var_idx_left[j]].col(0) + 
										node_x_vals[this->node_var_idx_left[j]].value.col(0);
			double theta = std::max(1-lambda/(sum.norm()*rho),0.0);
			//std::cout << "sum " << edge_u_vals[this->edge_var_idx_left[i]].col(l) << " " << node_x_vals[this->node_var_idx_left[i]].value.col(l) << " " << "eta " << 1.0/lambda << "\n";
			edge_z_vals[this->edge_var_idx_left[j]].col(0) = theta * sum ;
		}
		int i = p + 1;
		for ( int j = 1; j <= p; ++j ){
			for ( int k = 1; k <= p; ++k){
				if ( j == k ){
					edge_z_vals[this->edge_var_idx_left[i]] = node_x_vals[this->node_var_idx_left[i]].value + 
														edge_u_vals[this->edge_var_idx_left[i]];
				}
				else{
					for ( int l = 0; l < edge_z_vals[this->edge_var_idx_left[i]].cols(); ++l){
						Eigen::MatrixXd sum = edge_u_vals[this->edge_var_idx_left[i]].col(l) + 
										node_x_vals[this->node_var_idx_left[i]].value.col(l);
						double theta = std::max(1-lambda/(sum.norm()*rho),0.0);
						std::cout << "sum " << edge_u_vals[this->edge_var_idx_left[i]].col(l) << " " << node_x_vals[this->node_var_idx_left[i]].value.col(l) << " " << "eta " << 1.0/lambda << "\n";
						edge_z_vals[this->edge_var_idx_left[i]].col(l) = theta * sum ;
					}
				}
				std::cout << "z " << this->edge_var_idx_left[i] << " " << edge_z_vals[this->edge_var_idx_left[i]] << "\n";
				i++;
			}
		}
		for ( int j = 0; j < this->edge_var_idx_left.size(); ++j){
			edge_u_vals[this->edge_var_idx_left[j]] += node_x_vals[this->node_var_idx_left[j]].value -
												edge_z_vals[this->edge_var_idx_left[j]];
			std::cout << "u " << this->edge_var_idx_left[j] << " " << edge_u_vals[this->edge_var_idx_left[j]] << "\n";
		}
	}
}; 
#endif
