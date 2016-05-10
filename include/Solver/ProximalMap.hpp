#ifndef PROXIMAL_MAP_H_
#define PROXIMAL_MAP_H_

#include <map>
#include <string>
#include "Node.hpp"
#include "Edge.hpp"
#include "Square.hpp"
#include "ModSquare.hpp"
#include "NetLasso.hpp"

template<typename T> Node * createNodeInstance() { return new T; }
template<typename T> Edge * createEdgeInstance() { return new T; }

//include all the solver functions here
class ProximalMap{
protected:

	//static map for node objectives
	static std::map<std::string, Node*(*)()> node_map;
	static std::map<std::string, Edge*(*)()> edge_map;

public:
	static void LoadOperators(){
		node_map["SQUARE"] = &createNodeInstance<Square>;
		node_map["MOD_SQUARE"] = &createNodeInstance<ModSquare>;
		edge_map["NETLASSO"] = &createEdgeInstance<NetLasso>;
	}

	static Node *getNodeInstance(std::string op){
		return node_map[op]();
	}

	static Edge *getEdgeInstance(std::string op){
		return edge_map[op]();
	}
};
#endif

