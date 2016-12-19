#ifndef _NODE_MAPS_H_
#define _NODE_MAPS_H_

#include "TabaSCo.decl.h"

#include <vector>

class KaswellMap : public CkArrayMap
{
	public:
		KaswellMap(std::vector<int> nodeList);
		inline int procNum(int, const CkArrayIndex &iIndex);
	protected:
		std::vector<int> nodeList;
};

class HaswellMap : public CkArrayMap
{
	public:
		HaswellMap(std::vector<int> nodeList);
		inline int procNum(int, const CkArrayIndex &iIndex);
	protected:
		std::vector<int> nodeList;
};

#endif
