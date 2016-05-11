#include "ProximalMap.hpp"

std::map<std::string, Node*(*)()> ProximalMap::node_map;
std::map<std::string, Egde*(*)()> ProximalMap::edge_map;
