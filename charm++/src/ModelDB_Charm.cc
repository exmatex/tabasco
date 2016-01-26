#include "ModelDB_Charm.h"
#include "DBVecMessage.h"

#include "ResponsePoint.h"
#include "InterpolationModelFactory.h"
#include "InterpolationModel.h"

extern CProxy_DBInterface DBArray;

ModelDB_Charm::ModelDB_Charm(unsigned int index)
: dbVecIndex(index)
{

}

void ModelDB_Charm::insert(uint128_t & model_key, krigalg::InterpolationModelPtr krigingModel, krigcpl::ResponsePoint * point)
{
    std::vector<double> packedContainer;
    krigingModel->pack(*point, packedContainer);
    DBArray[dbVecIndex].push(model_key, packedContainer, point->size());
}

krigalg::InterpolationModelPtr ModelDB_Charm::extract(uint128_t & model_key, krigalg::InterpolationModelFactoryPointer  * newFact)
{
    ///TODO: Error checking would probably be smart
    DBVecMessage * msg = DBArray[dbVecIndex].pull(model_key);
	std::vector<double> packedContainer;
	///TODO: Determine the safe way to memcpy this as I have historically screwed it up
	for(int i = 0; i < msg->length; i++)
	{
		packedContainer.push_back(msg->data[i]);
	}
	krigalg::InterpolationModelPtr retPtr;
	retPtr = (*newFact)->build();
    retPtr->unpack(packedContainer);
    return retPtr;

}

void ModelDB_Charm::erase(uint128_t & model_key)
{
    //We had to comment this out with CoEVP/LULESH
    //DBArray[dbVecIndex].erase(model_key);
}
