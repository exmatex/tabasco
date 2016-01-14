#ifndef _MODEL_DB_INTERFACE_H_
#define _MODEL_DB_INTERFACE_H_

#include "TabaSCo.decl.h"

#include "ModelDatabase.h"

class ModelDBInterface : public ModelDatabase {
private:
    unsigned int dbVecIndex;
public:
    ModelDBInterface(unsigned int index);
    virtual void insert(uint128_t & model_key, krigalg::InterpolationModelPtr krigingModel, krigcpl::ResponsePoint * point);
    virtual krigalg::InterpolationModelPtr extract(uint128_t & model_key, krigalg::InterpolationModelFactoryPointer  * newFact);
    virtual void erase(uint128_t & model_key);
};

#endif
