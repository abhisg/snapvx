#include "Node.h"

class Square : public Node
{
public:
	Square():Node(){}

	void ADMM_node(Node * node, std::unordered_map<int,Node_Var> &node_x_vals,
					std::unordered_map<int,Eigen::MatrixXd> &edge_z_vals,
					std::unordered_map<int,Eigen::MatrixXd>	&edge_u_vals){
		for ( int j = 0 ; j < node->x_var_idx.size(); ++j ){
			Eigen::MatrixXd increment(2*node->args[j]["a"]);
			for ( int k = 0 ; k < node->neighbour_var_idx[j].size(); ++k ){
				increment += rho*(edge_z_vals[node->neighbour_var_idx[j][k]] - edge_u_vals[node->neighbour_var_idx[j][k]]);
			}
			node_x_vals[node->x_var_idx[j]].value = increment/(2+rho*node->neighbour_var_idx[j].size());
		}
	}
}