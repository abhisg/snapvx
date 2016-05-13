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
		int p = (int)((-1+sqrt(1+4*(this->x_var_idx.size()-1)))/2.0);
		int d = 0;
		for ( int j = 0 ; j < p + 1; ++j ){
			d += node_x_vals[this->x_var_idx[j]].value.rows();
		}
		Eigen::MatrixXd Theta(d,d);
		Eigen::MatrixXd Z(d,d);
		Eigen::MatrixXd U(d,d);
		Eigen::MatrixXd A(d,d);
		Z.block(0,0,1,1) = edge_z_vals[this->neighbour_var_idx[0][0]];
		U.block(0,0,1,1) = edge_u_vals[this->neighbour_var_idx[0][0]];
		int i = 1,k = 1;
		for ( int j = 1; j < p + 1; ++j ){
			int m = edge_z_vals[this->neighbour_var_idx[j][0]].rows();
			Z.block(0,i,1,m) = edge_z_vals[this->neighbour_var_idx[j][0]].transpose();
			Z.block(i,0,m,1) = edge_z_vals[this->neighbour_var_idx[j][0]];
			U.block(0,i,1,m) = edge_u_vals[this->neighbour_var_idx[j][0]].transpose();
			U.block(i,0,m,1) = edge_u_vals[this->neighbour_var_idx[j][0]];
			A.block(0,i,1,m) = this->args[j]["A"].transpose();
			A.block(i,0,m,1) = this->args[j]["A"];
			i += m;
		}
		i=p+1;
		std::vector<int> lr(p,1);
		for ( int j = 1;j <= p ; ++j ){
			int lc = 1;
			for ( int k = 1; k <= p ; ++k ){
				int m = edge_z_vals[this->neighbour_var_idx[i][0]].rows();
				int n = edge_z_vals[this->neighbour_var_idx[i][0]].cols();
				Z.block(lr[k-1],lc,m,n) = edge_z_vals[this->neighbour_var_idx[i][0]];
				U.block(lr[k-1],lc,m,n) = edge_u_vals[this->neighbour_var_idx[i][0]];
				A.block(lr[k-1],lc,m,n) = this->args[i]["A"];
				lr[k-1] += m;
				lc += n;
				i++;
			}
		}
		std::cout << "A" << A << "\n";
		Eigen::EigenSolver<Eigen::MatrixXd> decomp(rho*((Z-U+Z.transpose()-U.transpose())/2)-A);
		Eigen::MatrixXcd lambda = decomp.eigenvalues();
		for (int i = 0; i < lambda.rows(); ++i ){
			lambda(i,0) += sqrt(lambda(i,0) * lambda(i,0) + 4*rho);
		}
		Eigen::MatrixXd D = lambda.real().asDiagonal();
		Eigen::MatrixXd Q = decomp.eigenvectors().real();
		Theta = 1/(2*rho)*Q*D*Q.transpose();
		//std::cout << "Theta " << Theta << "\n";
		node_x_vals[this->x_var_idx[0]].value = Theta.block(0,0,1,1);
		i = 1;
		for ( int j = 1; j < p + 1; ++j){
			int m = node_x_vals[this->x_var_idx[j]].value.rows();
			node_x_vals[this->x_var_idx[j]].value = Theta.block(i,0,m,1);
			std::cout << "x " << this->x_var_idx[j] << " " << m << " " << node_x_vals[this->x_var_idx[j]].value << "\n";
			i += m;
		}
		i=p+1;
		lr = std::vector<int>(p,1);
		for ( int j = 1;j <= p ; ++j ){
			int lc = 1;
			for ( int k = 1; k <= p ; ++k ){
				int m = edge_z_vals[this->neighbour_var_idx[i][0]].rows();
				int n = edge_z_vals[this->neighbour_var_idx[i][0]].cols();
				node_x_vals[this->x_var_idx[i]].value = Theta.block(lr[k-1],lc,m,n);
				lr[k-1] += m;
				std::cout << "x " << this->x_var_idx[i] << " " << m << " " << n << " " << lr[k-1] << " " << lc << " " << node_x_vals[this->x_var_idx[i]].value << "\n";
				lc += n;
				i++;
			}
		}
	}
};
#endif