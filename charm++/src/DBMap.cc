#include "TabaSCo.decl.h"
#include "DBMap.h"
#include <cstring>
#include <cstddef>

/*
  Integration is fairly simple:

  First initialize the map
    CProxy_DBMap DBMapProxy = CProxy_DBMap::ckNew(intNumberOfDBNodes);
  or
    CProxy_DBMap DBMapProxy = CProxy_DBMap::ckNew(vectorWithDesiredPhysicalNodes);
  Then set up a CkArrayOptions to use the mapper for the desired number of DB Nodes
    CkArrayOptions mapOptions(intNumberOfDBNodes);
    mapOptions.setMap(DBMapProxy);
  Then init the Chare array as normal, replacing the dimensions/size with the DBMapProxy
    arrayProxy = CProxy_DBInterface::ckNew(mapOptions);
*/

DBMap::DBMap(int dbCount, int dbNodeCount)
{
    //Specify the number of Database Nodes/Chares
    this->dbNodeCount = dbNodeCount;
	this->dbCount = dbCount;
    this->useNodeList = false;
}

DBMap::DBMap(std::vector<int> nodeList)
{
    //Specify which Processor Ranks should have a DB Node/Chare
    this->dbNodeCount = nodeList.size();
    this->nodeList = nodeList;
    this->useNodeList = true;
}

int DBMap::procNum(int, const CkArrayIndex &iIndex)
{
    //*index is n-d array (1d for us) corresponding to array index
    int *index=(int *) iIndex.data();

    int proc;
    if(CmiCpuTopologyEnabled())
    { 
        // use charm++ physnode API if available
        //Do we know what ranks we want to have db chares on?
        if(this->useNodeList == true)
        {
            //Yes, so just do that
            if(index[0] < this->dbNodeCount)
            {
                //Get the first PE on the physical node
                proc = CmiGetFirstPeOnPhysicalNode(this->nodeList[index[0]]);
            }
            else
            {
                ///TODO: Probably output a warning as this is likely not intended behavior
                int vecIndex = index[0] % this->dbNodeCount;
                proc = CmiGetFirstPeOnPhysicalNode(this->nodeList[vecIndex]);
            }
        }
        else
        {
            //Just round-robin the physical nodes up to the desired count
			int nodeRange = std::min(CmiNumPhysicalNodes(), this->dbNodeCount);
            int physNode = index[0] % nodeRange;
            proc = CmiGetFirstPeOnPhysicalNode(physNode);
            ///WARNING: Is there an issue with pinning multiple chares to a single PE?
        }
    }
    else
    {
        //No physnode API, so we don't know what physical node is what rank, but we may be able to work around that
        int numNodes = CmiNumNodes();
        int nodeSize = CmiMyNodeSize();
        //Did we manually map what physical nodes should have DB Chares?
        if(this->useNodeList == true)
        {
            //Yes, so just do that
            if(index[0] < this->dbNodeCount)
            {
                proc = this->nodeList[index[0]] * nodeSize;
            }
            else
            {
                ///TODO: Probably output a warning as this is not intended behavior
                int vecIndex = index[0] % this->dbNodeCount;
                proc = this->nodeList[vecIndex] * nodeSize;
            }
        }
        else
        {
            //Just round-robin to assign
            proc = index[0] * nodeSize;;
            ///WARNING: Is there an issue with pinning multiple chares to a single PE?
        }
    }
    return proc;
}
