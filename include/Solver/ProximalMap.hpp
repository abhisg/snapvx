#ifndef PROXIMAL_MAP_H_
#define PROXIMAL_MAP_H_

#include <unordered_map>
#include "Node.hpp"
#include "Edge.hpp"
#include "Square.hpp"
#include "ModSquare.hpp"
#include "NetLasso.hpp"

typedef enum ProximalOperator
{
        SQUARE,
        MOD_SQUARE,
        NETLASSO,
        NONE
}ProximalOperator;

//include all the solver functions here

//static map for node objectives
template<typename T> Node * createNodeInstance() { return new T; }
typedef std::map<ProximalOperator, Node*(*)()> node_map_type;
node_map_type nodemap;
nodemap[SQUARE] = &createInstance<Square>;
nodemap[MOD_SQUARE] = &createInstance<ModSquare>;

//static map for edge objectives
template<typename T> Edge * createNodeInstance() { return new T; }
typedef std::map<ProximalOperator, Edge*(*)()> node_map_type;
edge_map_type edgemap;
edgemap[NETLASSO] = &createInstance<NetLasso>;


Node *getNodeInstance(ProximalOperator op){
	return nodemap[op];
}

Edge *getEdgeInstance(ProximalOperator op){
	return edgemap[op];
}
#endif

