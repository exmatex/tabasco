#ifndef _FINE_SCALE_MODEL_H_
#define _FINE_SCALE_MODEL_H_

#include "TabaSCo.decl.h"

#include "ElastoViscoPlasticity.h"
#include "ResponsePoint.h"
#include "tensor.h"

#ifndef NNS_AS_CHARE
#ifdef FLANN
#include "ApproxNearestNeighborsFLANN.h"
#else
#include "ApproxNearestNeighborsMTree.h"
#endif
#endif

// PUP operator for Tensor2Sym
inline void operator|(PUP::er &p, Tensor2Sym &tensor)
{
  PUParray(p, tensor.a, 6);
}

// PUP operator for Tensor2Gen
inline void operator|(PUP::er &p, Tensor2Gen &tensor)
{
  PUParray(p, tensor.a, 9);
}

// PUP operator for ConstitutiveData
inline void operator|(PUP::er &p, ConstitutiveData &cdata)
{
  p|cdata.sigma_prime;
  p|cdata.num_models;
  p|cdata.num_point_value_pairs;
  p|cdata.num_Newton_iters;

}

class FineScaleModel : public CBase_FineScaleModel {
  private:
    FineScaleModel_SDAG_CODE
    int newPt;
    int currentIter;
    int nbrCount;
    int nbrData;
    int useAdaptiveSampling;

    int nnsIndex;
    int interpIndex;
    int dbIndex;

    size_t stateSize;
       
  public:
  
  Constitutive* cm; 
  ElastoViscoPlasticity* em;

// Keep this for now, just to be compatible with CoEVP
//#ifndef NNS_AS_CHARE
  ApproxNearestNeighbors* ann;
//#endif
  ModelDatabase * modelDB;

  FineScaleModel();
  FineScaleModel(int state_size, bool use_adaptive_sampling, int nnsIndex, int interpIndex, int dbIndex);
  FineScaleModel(CkMigrateMessage *msg);
  ~FineScaleModel();
  void pup(PUP::er &p);


  // Entry methods
  void advance(const double delta_t, const Tensor2Gen& L_new, const double volume_change, int ssize, char* state);

  // Other methods
  void updateVbar_prime( const Tensor2Sym& Vbar_prime_old,
                         const Tensor2Sym& Dprime_new,
                         const Tensor2Gen& R_new,
                         const double      a_new,
                         const double      delta_t,
                         Tensor2Sym&       Vbar_prime_new,
                         Tensor2Sym&       Dbar_prime,
                         Tensor2Gen&       Wbar );

  void evaluateFineScaleModel( const Tensor2Sym& tau_bar_prime,
                               Tensor2Sym&       Dbar_prime,
                               Tensor4LSym&      Dbar_prime_deriv );

  void sample( const FineScale&           fine_scale_model,
               const std::vector<double>& point,
               std::vector<double>&       value );

  void sample( std::vector<double>&       value,
                const std::vector<double>& point,
                int&                       hint,
                std::vector<bool>&         flags,
                const FineScale&           fine_scale_model,
                double&                    error_estimate);

  bool interpolate(double            * value,
                   int               & hint,
                   const double      * point,
                   std::vector<bool> & flags,
                   double            & error_estimate);

  bool interpolate(double            * value,
                   double            * gradient,
                   int               & hint,
                   const double      * point,
                   std::vector<bool> & flags,
                   double            & error_estimate);

  void interpolate(double                                * value,
                    double                                * gradient,
                    InterpolationModelPtr                   krigingModel,
                    const ResponsePoint                   & queryPoint,
                    int                                     pointDimension,
                    int                                     valueDimension);

  void insert(int               & hint,
              const double      * point,
              const double      * value,
              const double      * gradient,
              std::vector<bool> & flags);

  Point getModelCenterMass(const InterpolationModel & krigingModel);

  std::pair<int, InterpolationModelPtr> 
      findClosestCoKrigingModel(const ResponsePoint        & point,
#ifndef NNS_AS_CHARE
                                ApproxNearestNeighbors     & ann,
#endif
                                krigalg::InterpolationModelFactoryPointer modelFactory,
                                ModelDatabase * modelDB,
                                double                       maxQueryPointModelDistance);

  std::pair<int, InterpolationModelPtr>
      findBestCoKrigingModel(bool &                       canInterpolateFlag,
                             const ResponsePoint &        point,
#ifndef NNS_AS_CHARE
                             ApproxNearestNeighbors     & ann,
#endif
                             ModelDatabase * modelDB,
                             const InterpolationModelFactoryPointer& _modelFactory,
                             double                       tolerance,
                             double                       meanErrorFactor,
                             double                       maxQueryPointModelDistance,
                             int                          maxNumberSearchModels,
                             int                          valueDimension);

  std::pair<int, InterpolationModelPtr>
      findBestCoKrigingModel(bool &                       canInterpolateFlag,
                             const ResponsePoint &        point,
#ifndef NNS_AS_CHARE
                             ApproxNearestNeighbors     & ann,
#endif
                             ModelDatabase * modelDB,
                             const InterpolationModelFactoryPointer& _modelFactory,
                             double                       tolerance,
                             double                       meanErrorFactor,
                             double                       maxQueryPointModelDistance,
                             int                          maxNumberSearchModels,
                             int                          maxKrigingModelSize,
                             int                          valueDimension);


  double checkError(InterpolationModelPtr   krigingModel,
                    const ResponsePoint   & queryPoint,
                    int                     valueDimension,
                    double                 _meanErrorFactor);

  double checkErrorInf(InterpolationModelPtr   krigingModel,
                       const ResponsePoint   & queryPoint,
                       int                     valueDimension,
                       double                 _meanErrorFactor);

  double compKrigingError(const InterpolationModel & krigingModel,
                          const Point              & queryPoint,
                          int                        valueId,
                          double                     meanErrorFactor);

  bool checkErrorAndInterpolate(double               * value,
                                InterpolationModelPtr  krigingModel,
                                const ResponsePoint  & queryPoint,
                                int                    valueDimension,
                                double                _tolerance,
                                double                _meanErrorFactor,
                                double               & errorEstimate);

  bool checkErrorAndInterpolate(double               * value,
                                double               * gradient,
                                InterpolationModelPtr  krigingModel,
                                const ResponsePoint  & queryPoint,
                                int                    pointDimension,
                                int                    valueDimension,
                                double                _tolerance,
                                double                _meanErrorFactor,
                                double               & errorEstimate );

  void addNewModel(
                    ModelDatabase *             modelDB,
#ifndef NNS_AS_CHARE
                    ApproxNearestNeighbors&                  ann,
#endif
                    const InterpolationModelFactoryPointer & _modelFactory,
                    int &                                    objectId,
                    const double *                           pointData,
                    const double *                           valueData,
                    const double *                           gradientData,
                    int                                      pointDimension,
                    int                                      valueDimension);


  std::vector<Value>
      copyValueData(const double * value,
                    const double * gradient,
                    int            pointDimension,
                    int            valueDimension);

  std::vector<Value>
      copyValueData(const double * value,
                    int            pointDimension,
                    int            valueDimension);

  uint128_t getKeyHash(const ResponsePoint& point);

};
#endif
