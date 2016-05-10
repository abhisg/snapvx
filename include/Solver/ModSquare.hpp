#include "Node.h"

class Square : public Node
{
public:
	Square():Node(){}

	void ADMM_node(Node * node, std::unordered_map<int,Node_Var> &node_x_vals,
					std::unordered_map<int,Eigen::MatrixXd> &edge_z_vals,
					std::unordered_map<int,Eigen::MatrixXd>	&edge_u_vals){
		for ( int j = 0 ; j < node->x_var_idx.size(); ++j ){
			Eigen::MatrixXd increment(node->args[j]["rhs"]);
			for ( int k = 0 ; k < node->neighbour_var_idx[j].size(); ++k ){
				increment += rho/2*(edge_z_vals[node->neighbour_var_idx[j][k]] - edge_u_vals[node->neighbour_var_idx[j][k]]);
			}
			node_x_vals[node->x_var_idx[j]].value = node->args[j]["lhs"]*increment;
			//std::cout<<"x " << node->x_var_idx[j] <<" "<<node_x_vals[node->x_var_idx[j]].value<<"\n";
		}
	}
}