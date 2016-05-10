#ifndef EDGE_H_
#define EDGE_H_
#include "../CVXcanon/src/CVXcanon.hpp"
#include "NodeVar.hpp"
#include <vector>

class Edge
{
	public:
		LinOp* edge_objective;
		std::vector<LinOp *> edge_constraints;	
		std::vector<int> edge_var_idx_left;
		std::vector<int> edge_var_idx_right;
		std::vector<int> node_var_idx_left;
		std::vector<int> node_var_idx_right;

		virtual std::vector<double> ADMM_edge(std::unordered_map<int,Node_Var> &,
							std::unordered_map<int,Eigen::MatrixXd> &,
							std::unordered_map<int,Eigen::MatrixXd>	&,
							double &) = 0;
};
#endif
