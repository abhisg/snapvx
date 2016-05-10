#ifndef PROXIMAL_MAP_H_
#define PROXIMAL_MAP_H_

#include <unordered_map>
//include all the solver functions here
typedef enum ProximalOperator
{
        SQUARE,
        MOD_SQUARE,
        NETLASSO,
        NONE
}ProximalOperator;

template<typename T> Base * createInstance() { return new T; }

typedef std::map<ProximalOperator, Base*(*)()> map_type;

map_type map;
map[SQUARE] = &createInstance<DerivedA>;
map[MOD_SQUARE] = &createInstance<DerivedB>;

#endif
