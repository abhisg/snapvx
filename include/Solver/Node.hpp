#ifndef NODE_H_
#define NODE_H_
#include "../CVXcanon/src/CVXcanon.hpp"
#include "../CVXcanon/include/Eigen/Eigenvalues"
#include "NodeVar.hpp"
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

		Node()
		{
			node_objective = NULL;
			node_constraints = std::vector<LinOp *>();
			x_var_idx = std::vector<int>();
			args = std::vector<std::map<std::string,Eigen::MatrixXd> >();
		}

		virtual void ADMM_node(std::unordered_map<int,Node_Var> &,
						std::unordered_map<int,Eigen::MatrixXd> &,
						std::unordered_map<int,Eigen::MatrixXd>	&,
						double &) = 0;

		virtual void LoadNodeProximal(std::vector<int>  &x_var_idx,
									std::vector<std::vector<int> >  &neighbour_var_idx,
									std::vector<int>  &sizes_i,
									std::vector<int>  &sizes_j,
									std::vector<std::map<std::string,Eigen::MatrixXd > >  &args)
		{
			this->node_objective = NULL;
			this->node_constraints = std::vector<LinOp *>();
			this->neighbour_var_idx = neighbour_var_idx;
			this->x_var_idx = x_var_idx;
			this->args = std::vector<std::map<std::string,Eigen::MatrixXd> >();
			for ( int j = 0 ; j < static_cast<int>(x_var_idx.size()); ++j){
				this->args.push_back(args[j]);
			}
			std::cout<<"Inside parent func\n";
		}
};
#endif
