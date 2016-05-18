#ifndef NODEVAR_H_
#define NODEVAR_H_
#include "../CVXcanon/src/CVXcanon.hpp"
#include <string>

typedef struct Node_Var
{
	Eigen::MatrixXd value;
	std::string name;
	int nodeId;
}Node_Var;
#endif