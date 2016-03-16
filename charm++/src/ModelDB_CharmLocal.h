#ifndef _MODEL_DB_CHARMLOCAL_H_
#define _MODEL_DB_CHARMLOCAL_H_

#include "TabaSCo.decl.h"

#include "ModelDatabase.h"
#include "ModelDB_Enums.h"
#include "ModelDB_SharedDB.h"

class ModelDB_CharmLocal : public ModelDatabase {
private:
    SingletonDBBackendEnum dbType;
	ModelDB_SharedDB * dbRef;
public:
    ModelDB_CharmLocal(SingletonDBBackendEnum dbType);
    virtual void insert(uint128_t & model_key, krigalg::InterpolationModelPtr krigingModel, krigcpl::ResponsePoint * point);
    virtual krigalg::InterpolationModelPtr extract(uint128_t & model_key, krigalg::InterpolationModelFactoryPointer  * newFact);
    virtual void erase(uint128_t & model_key);
};

#endif
