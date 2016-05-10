#ifndef SQUARE_H_
#define SQUARE_H_
#include "Node.hpp"

class Square : public Node
{
public:
	Square():Node(){}

	virtual void ADMM_node(std::unordered_map<int,Node_Var> &node_x_vals,
					std::unordered_map<int,Eigen::MatrixXd> &edge_z_vals,
					std::unordered_map<int,Eigen::MatrixXd>	&edge_u_vals,
					double &rho){
		for ( int j = 0 ; j < this->x_var_idx.size(); ++j ){
			Eigen::MatrixXd increment(2*this->args[j]["a"]);
			for ( int k = 0 ; k < this->neighbour_var_idx[j].size(); ++k ){
				increment += rho*(edge_z_vals[this->neighbour_var_idx[j][k]] - edge_u_vals[this->neighbour_var_idx[j][k]]);
			}
			node_x_vals[this->x_var_idx[j]].value = increment/(2+rho*this->neighbour_var_idx[j].size());
		}
	}
};
#endif