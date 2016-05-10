#include "CVXcanon/src/CVXcanon.hpp"
#include <string>

typedef struct Node_Var
{
	Eigen::MatrixXd value;
	std::string name;
	int nodeId;
}Node_Var;