#ifndef ADMM_H
#define ADMM_H

#include "CVXcanon/src/CVXcanon.hpp"
#include <map>
#include <string>
#include <utility>

class ADMM
{
	protected:
		std::vector<LinOp* > node_objectives;
		std::vector<LinOp* > edge_objectives;
		std::vector<std::vector< LinOp *> > node_constraints;
		std::vector<std::vector< LinOp *> > edge_constraints;
		std::vector<Solution *> node_x_vals;
		std::vector<Solution *> edge_z_vals;
		std::vector<Solution *> edge_u_vals;
		std::map<std::string,double > solver_options;
		void ADMM_x(int);
		void ADMM_z(int);
		void ADMM_u(int);
		//void copyLinop(LinOp* const&,LinOp* const&);
	public:
		ADMM();
		void LoadNodes(std::vector<LinOp* > &,std::vector<std::vector< LinOp *> > &);
		void LoadEdges(std::vector<LinOp* > &,std::vector<std::vector< LinOp *> > &);
		void Solve();
		std::vector<std::map<int, Eigen::MatrixXd> > get_node_x_vals();
		std::vector<std::map<int, Eigen::MatrixXd> > get_edge_z_vals();
		std::vector<std::map<int, Eigen::MatrixXd> > get_edge_u_vals();
};

#endif		
//std::map<int,Solution> Solve(Sense sense, std::map<int,LinOp* objective> >,std::map<int,std::vector< LinOp* > constraints >);

