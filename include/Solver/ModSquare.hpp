#include "Node.hpp"

class Square : public Node
{
public:
	Square():Node(){}

	void ADMM_node(std::unordered_map<int,Node_Var> &node_x_vals,
					std::unordered_map<int,Eigen::MatrixXd> &edge_z_vals,
					std::unordered_map<int,Eigen::MatrixXd>	&edge_u_vals,
					double &rho){
		for ( int j = 0 ; j < this->x_var_idx.size(); ++j ){
			Eigen::MatrixXd increment(this->args[j]["rhs"]);
			for ( int k = 0 ; k < this->neighbour_var_idx[j].size(); ++k ){
				increment += rho/2*(edge_z_vals[this->neighbour_var_idx[j][k]] - edge_u_vals[this->neighbour_var_idx[j][k]]);
			}
			node_x_vals[this->x_var_idx[j]].value = this->args[j]["lhs"]*increment;
			//std::cout<<"x " << node->x_var_idx[j] <<" "<<node_x_vals[node->x_var_idx[j]].value<<"\n";
		}
	}
}