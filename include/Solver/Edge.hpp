#ifndef EDGE_H_
#define EDGE_H_
#include "../CVXcanon/src/CVXcanon.hpp"
#include "NodeVar.hpp"
#include <vector>
#include <unordered_map>

class Edge
{
	public:
		LinOp* edge_objective;
		std::vector<LinOp *> edge_constraints;	

		std::vector<int> edge_var_idx_left;
		std::unordered_map<int,Eigen::MatrixXd> z_ij;

		std::vector<int> edge_var_idx_right;
		std::unordered_map<int,Eigen::MatrixXd> z_ji;

		std::vector<int> node_var_idx_left;

		std::vector<int> node_var_idx_right;

		virtual void ADMM_edge(std::unordered_map<int,Node_Var> &,
							std::unordered_map<int,Eigen::MatrixXd> &,
							std::unordered_map<int,Eigen::MatrixXd>	&,
							double &) = 0;

		virtual std::vector<double> ADMM_edge_solver(std::unordered_map<int,Node_Var> &node_x_vals,
							std::unordered_map<int,Eigen::MatrixXd> &edge_z_vals,
							std::unordered_map<int,Eigen::MatrixXd>	&edge_u_vals,
							double &rho)
		{
			for ( int j = 0 ; j < this->edge_var_idx_left.size(); ++j ){
				//if ( this->edge_var_idx_left[j] != 0 || this->edge_var_idx_right[j] != 0 ){
					z_ij[this->edge_var_idx_left[j]] = edge_z_vals[this->edge_var_idx_left[j]];
					//z_ji[this->edge_var_idx_right[j]] = edge_z_vals[this->edge_var_idx_right[j]];
				//}
			}
			ADMM_edge(node_x_vals,edge_z_vals,edge_u_vals,rho);
			std::vector<double> norms(5,0);
			for ( int j = 0 ; j < this->edge_var_idx_left.size(); ++j ){
				if ( this->edge_var_idx_left[j] != 0 || this->edge_var_idx_right[j] != 0 ){
					norms[0] += (node_x_vals[this->node_var_idx_left[j]].value - edge_z_vals[this->edge_var_idx_left[j]]).squaredNorm();
							 //+ (node_x_vals[this->node_var_idx_right[j]].value - edge_z_vals[this->edge_var_idx_right[j]]).squaredNorm();
					norms[1] += (edge_z_vals[this->edge_var_idx_left[j]] - z_ij[this->edge_var_idx_left[j]] ).squaredNorm();
						//+(edge_z_vals[this->edge_var_idx_right[j]] - z_ji[this->edge_var_idx_right[j]]).squaredNorm());
					norms[2] +=  edge_z_vals[this->edge_var_idx_left[j]].squaredNorm();
						//+ edge_z_vals[this->edge_var_idx_right[j]].squaredNorm();
					norms[3] += node_x_vals[this->node_var_idx_left[j]].value.squaredNorm();
							 //+node_x_vals[this->node_var_idx_right[j]].value.squaredNorm();	
					norms[4] += edge_u_vals[this->edge_var_idx_left[j]].squaredNorm();
						//+edge_u_vals[this->edge_var_idx_right[j]].squaredNorm());
				}
			}
			return norms;
		}

};
#endif
