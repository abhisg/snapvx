#ifndef PROXIMAL_MAP_H_
#define PROXIMAL_MAP_H_

#include <map>
#include <string>
#include "Node.hpp"
#include "Edge.hpp"

//node solvers
#include "Square.hpp"
#include "ModSquare.hpp"
#include "NetLaplace.hpp"

//edge solvers
#include "NetLasso.hpp"
#include "EdgeLasso.hpp"

template<typename T> Node * createNodeInstance() { return new T; }
template<typename T> Edge * createEdgeInstance() { return new T; }
static std::map<std::string, Node*(*)()> node_map;
static std::map<std::string, Edge*(*)()> edge_map;

//include all the solver functions here
class ProximalMap{

public:
	//static map for node objectives

	static void LoadOperators(){
		node_map["SQUARE"] = &createNodeInstance<Square>;
		node_map["MOD_SQUARE"] = &createNodeInstance<ModSquare>;
		node_map["NETLAPLACE"] = &createNodeInstance<NetLaplace>;
		edge_map["NETLASSO"] = &createEdgeInstance<NetLasso>;
		edge_map["EDGELASSO"] = &createEdgeInstance<EdgeLasso>;
	}

	static Node *getNodeInstance(std::string op){
		if ( node_map.find(op) != node_map.end() ){
			return node_map[op]();
		}
		else{
			return NULL;
		}
	}

	static Edge *getEdgeInstance(std::string op){
		if ( edge_map.find(op) != edge_map.end() ){
			return edge_map[op]();
		}
		else{
			return NULL;
		}
	}
};
#endif

