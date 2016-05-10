#ifndef NODE_H_
#define NODE_H_
#include "CVXcanon/src/CVXcanon.hpp"
#include "NodeVar.h"
#include <unordered_map>
#include <vector>


class Node 
{
	public:
		LinOp* node_objective;
		std::vector<LinOp *> node_constraints;
		std::vector<std::vector<int> > neighbour_var_idx;
		std::vector<int> x_var_idx;
		std::vector<std::map<std::string,Eigen::MatrixXd> > args;

		virtual void ADMM_node(std::unordered_map<int,Node_Var> &,
						std::unordered_map<int,Eigen::MatrixXd> &,
						std::unordered_map<int,Eigen::MatrixXd>	&,
						double) = 0;
};
#endif