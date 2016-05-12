#ifndef NETLAPLACE_H_
#define NETLAPLACE_H_
#include "Node.hpp"

class NetLaplace : public Node
{
public:
	NetLaplace():Node(){}

	virtual void ADMM_node(std::unordered_map<int,Node_Var> &node_x_vals,
					std::unordered_map<int,Eigen::MatrixXd> &edge_z_vals,
					std::unordered_map<int,Eigen::MatrixXd>	&edge_u_vals,
					double &rho){
		//construct theta,Z and U
		int p = (int)((-1+sqrt(1+4*this->x_var_idx.size()))/2.0);
		int d = 0;
		for ( int j = 0 ; j < p; ++j ){
			d += node_x_vals[this->x_var_idx[j]].value.rows();
		}
		Eigen::MatrixXd Theta(d+1,d+1);
		Eigen::MatrixXd Z(d+1,d+1);
		Eigen::MatrixXd U(d+1,d+1);
		Eigen::MatrixXd A(d+1,d+1);
		int i = 0;
		for ( int j = 0; j < p; ++j ){
			int m = edge_z_vals[this->neighbour_var_idx[j][0]].rows();
			Z.block(i,0,1,m) = edge_z_vals[this->neighbour_var_idx[j][0]].transpose();
			Z.block(0,i,m,1) = edge_z_vals[this->neighbour_var_idx[j][0]];
			U.block(i,0,1,m) = edge_u_vals[this->neighbour_var_idx[j][0]].transpose();
			U.block(0,i,m,1) = edge_u_vals[this->neighbour_var_idx[j][0]];
			A.block(i,0,1,m) = this->args[j]["A"].transpose();
			A.block(0,i,m,1) = this->args[j]["A"];
			i += m;
		}
		i=0,k=0;
		for ( int j = p; j < this->x_var_idx.size(); ++j){
			int m = edge_z_vals[this->neighbour_var_idx[j][0]].rows();
			int n = edge_z_vals[this->neighbour_var_idx[j][0]].cols();
			Z.block(i,k,m,n) = edge_z_vals[this->neighbour_var_idx[j][0]];
			U.block(i,k,m,n) = edge_u_vals[this->neighbour_var_idx[j][0]];
			A.block(i,k,m,n) = this->args[j]["A"];
			i += m;
			k += n;
		}
		Eigen::EigenSolver decomp((Z-U+Z.transpose()-U.transpose())/2-A);
		MatrixXcd lambda = decomp.eigenvalues();
		for (int i = 0; i < lambda.row(); ++i ){
			lambda(i,0) += sqrt(lambda(i,0) * lambda(i,0) + 4*rho);
		}
		Eigen::MatrixXd D = lambda.asDiagonal().real();
		Eigen::MatrixXd Q = es.eigenvectors().real();
		Theta = 1/(2*rho)*Q*D*Q.transpose();
		i = 0;
		for ( int j = 0; j < p; ++j){
			int m = node_x_vals[this->x_var_idx[j]].rows();
			node_x_vals[this->x_var_idx[j]] = Theta.block(i,0,m,1);
			i += m;
		}
		i=0,k=0;
		for ( int j = p; j < this->x_var_idx.size(); ++j){
			int m = node_x_vals[this->x_var_idx[j]].rows();
			int n = node_x_vals[this->x_var_idx[j]].cols();
			node_x_vals[this->x_var_idx[j]] = Theta.block(i,k,m,n);
			i += m;
			k += n;
		}
	}
};
#endif