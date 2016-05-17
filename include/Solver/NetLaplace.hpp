#ifndef NETLAPLACE_H_
#define NETLAPLACE_H_
#include "Node.hpp"

class NetLaplace : public Node
{
protected:
	int p; //number of unknowns
	int d;
	Eigen::MatrixXd Theta;
	Eigen::MatrixXd Z;
	Eigen::MatrixXd U;
	Eigen::MatrixXd A;

public:

	NetLaplace():Node()
	{
		p = 0;
		d = 0;
		Theta = Eigen::MatrixXd(1,1);
		Z = Eigen::MatrixXd(1,1);
		U = Eigen::MatrixXd(1,1);
		A = Eigen::MatrixXd(1,1);
	}

	virtual void LoadNodeProximal(std::vector<int>  &x_var_idx,
									std::vector<std::vector<int> >  &neighbour_var_idx,
									std::vector<int>  &sizes_i,
									std::vector<int>  &sizes_j,
									std::vector<std::map<std::string,Eigen::MatrixXd > >  &args)
	{
		Node::LoadNodeProximal(x_var_idx,neighbour_var_idx,sizes_i,sizes_j,args);
		std::cout<<"Inside child func\n";
		p = static_cast<int>((-1+sqrt(1+4*(this->x_var_idx.size()-1)))/2.0);
		for ( int j = 0 ; j < p + 1; ++j ){
			d += sizes_i[j];
		}
		Theta = Eigen::MatrixXd(d,d);
		Z = Eigen::MatrixXd(d,d);
		U = Eigen::MatrixXd(d,d);
		A = Eigen::MatrixXd(d,d);

		A.block(0,0,1,1) = this->args[0]["A"];
		int i = 1;
		for ( int j = 1; j < p + 1; ++j ){
			int m = sizes_i[j];
			A.block(0,i,1,m) = this->args[j]["A"].transpose();
			A.block(i,0,m,1) = this->args[j]["A"];
			i += m;
		}
		i=p+1;
		std::vector<int> lr(p,1);
		for ( int j = 1;j <= p ; ++j ){
			int lc = 1;
			for ( int k = 1; k <= p ; ++k ){
				int m = sizes_i[i];
				int n = sizes_j[i];
				A.block(lr[k-1],lc,m,n) = this->args[i]["A"];
				lr[k-1] += m;
				lc += n;
				i++;
			}
		}
		std::cout << "A: " << A << "\n";
	}

	virtual void ADMM_node(std::unordered_map<int,Node_Var> &node_x_vals,
					std::unordered_map<int,Eigen::MatrixXd> &edge_z_vals,
					std::unordered_map<int,Eigen::MatrixXd>	&edge_u_vals,
					double &rho){
		//construct theta,Z and U
		/*int p = (int)((-1+sqrt(1+4*(this->x_var_idx.size()-1)))/2.0);
		int d = 0;
		for ( int j = 0 ; j < p + 1; ++j ){
			d += node_x_vals[this->x_var_idx[j]].value.rows();
		}
		Eigen::MatrixXd Theta(d,d);
		Eigen::MatrixXd Z(d,d);
		Eigen::MatrixXd U(d,d);
		Eigen::MatrixXd A(d,d);*/
		Z.block(0,0,1,1) = edge_z_vals[this->neighbour_var_idx[0][0]];
		U.block(0,0,1,1) = edge_u_vals[this->neighbour_var_idx[0][0]];
		//A.block(0,0,1,1) = this->args[0]["A"];
		int i = 1;
		for ( int j = 1; j < p + 1; ++j ){
			int m = edge_z_vals[this->neighbour_var_idx[j][0]].rows();
			Z.block(0,i,1,m) = edge_z_vals[this->neighbour_var_idx[j][0]].transpose();
			Z.block(i,0,m,1) = edge_z_vals[this->neighbour_var_idx[j][0]];
			U.block(0,i,1,m) = edge_u_vals[this->neighbour_var_idx[j][0]].transpose();
			U.block(i,0,m,1) = edge_u_vals[this->neighbour_var_idx[j][0]];
			/*A.block(0,i,1,m) = this->args[j]["A"].transpose();
			A.block(i,0,m,1) = this->args[j]["A"];*/
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
				//A.block(lr[k-1],lc,m,n) = this->args[i]["A"];
				lr[k-1] += m;
				lc += n;
				i++;
			}
		}
		//std::cout << "A: " << A << "\n";
		//std::cout << "matrix: " << rho*((Z-U+Z.transpose()-U.transpose())/2)-A << "\n";
		Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> decomp(rho*((Z-U+Z.transpose()-U.transpose())/2)-A);
		Eigen::VectorXd lambda = decomp.eigenvalues();
		for (int i = 0; i < lambda.rows(); ++i ){
			lambda(i,0) = (lambda(i,0)+sqrt(lambda(i,0) * lambda(i,0) + 4*rho))/(2*rho);
		}
		Eigen::MatrixXd D = lambda.asDiagonal();
		Eigen::MatrixXd Q = decomp.eigenvectors();
		Theta = Q*D*Q.transpose();
		std::cout << "Theta " << Theta << "\n";
		node_x_vals[this->x_var_idx[0]].value = Theta.block(0,0,1,1);
		i = 1;
		for ( int j = 1; j < p + 1; ++j){
			int m = node_x_vals[this->x_var_idx[j]].value.rows();
			node_x_vals[this->x_var_idx[j]].value = Theta.block(i,0,m,1);
			//std::cout << "x " << this->x_var_idx[j] << " " << m << " " << node_x_vals[this->x_var_idx[j]].value << "\n";
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
				//std::cout << "x " << this->x_var_idx[i] << " " << m << " " << n << " " << lr[k-1] << " " << lc << " " << node_x_vals[this->x_var_idx[i]].value << "\n";
				lc += n;
				i++;
			}
		}
	}
};
#endif