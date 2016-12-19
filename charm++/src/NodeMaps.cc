#include "TabaSCo.decl.h"
#include "NodeMaps.h"

KaswellMap::KaswellMap(std::vector<int> procList)
{
	///TODO: Is this safe?
	this->nodeList = procList;
}

int KaswellMap::procNum(int, const CkArrayIndex &iIndex)
{
	int nelem = *(int*)iIndex.data();
	return this->nodeList[nelem % this->nodeList.size()];
	
	
	///TODO: This is all not usable
	/*
	//This is the Kaswell map. We want one per Kaswell node
	int *index=(int *) iIndex.data();
	
	int proc;
	//Use physnode API if available
	if(CmiCpuTopologyEnabled())
	{
		if(index[0] < this->nodeList.size())
		{
			//Get the first PE on the physical node
			proc = CmiGetFirstPeOnPhysicalNode(this->nodeList[index[0]]);
		}
		else
		{
			///TODO: Probably output a warning as this is likely not intended behavior
			int vecIndex = index[0] % this->nodeList.size();
			proc = CmiGetFirstPeOnPhysicalNode(this->nodeList[vecIndex]);
		}
	}
	else
	{
		//No physnode API, so we don't know what physical node is what rank, but we may be able to work around that
		int numNodes = CmiNumNodes();
		int nodeSize = CmiMyNodeSize();
		//Did we manually map what physical nodes should have DB Chares?
		//Yes, so just do that
		if(index[0] < this->nodeList.size())
		{
			proc = this->nodeList[index[0]] * nodeSize;
		}
		else
		{
			///TODO: Probably output a warning as this is not intended behavior
			int vecIndex = index[0] % this->nodeList.size();
			proc = this->nodeList[vecIndex] * nodeSize;
		}
	}
	*/
}

HaswellMap::HaswellMap(std::vector<int> procList)
{
	this->nodeList = procList;
	/*
	//We need to expand nodelist
	for(int i = 0; i < procList.size(); i++)
	{
		int * pelist;
		int num = 2;
		CmiGetPesOnPhysicalNode(procList[i], &pelist, &num);
		for(int j = 0; j < num; j++)
		{
			this->nodeList.push_back(pelist[j]);
		}
	}
	*/
}

int HaswellMap::procNum(int, const CkArrayIndex &iIndex)
{
	int nelem = *(int*)iIndex.data();
	return this->nodeList[nelem % this->nodeList.size()];
	/*
	int *index=(int *) iIndex.data();
	
	int proc;
	//Use physnode API if available
	if(CmiCpuTopologyEnabled())
	{
		if(index[0] < this->nodeList.size())
		{
			proc = this->nodeList[index[0]];
		}
		else
		{
			int vecIndex = index[0] % this->nodeList.size();
			proc = this->nodeList[vecIndex];
		}
	}
	else
	{
		//No physnode API, so we don't know what physical node is what rank, but we may be able to work around that
		int numNodes = CmiNumNodes();
		int nodeSize = CmiMyNodeSize();
		//Did we manually map what physical nodes should have DB Chares?
		//Yes, so just do that
		if(index[0] < this->nodeList.size())
		{
			proc = this->nodeList[index[0]];
		}
		else
		{
			int vecIndex = index[0] % this->nodeList.size();
			proc = this->nodeList[vecIndex];
		}
	}
	return proc;
	*/
}


