#include "TabaSCo.decl.h"
#include "FineScaleModel.h"
#include "NearestNeighborSearch.h"
#include "Interpolate.h"
#include "Evaluate.h"
#include "DBInterface.h"

#include "Taylor.h"
#include "MieGruneisen.h"
#include "KrigingDataBase.h"
#include "InterpolationDataBase.h"
#include "ResponsePoint.h"
#include "Database.h"
#include "MurmurHash3.h"

#include "ModelDB_HashMap.h"

#define MURMUR_SEED 42

extern CProxy_Main mainProxy;
extern CProxy_CoarseScaleModel coarseScaleArray;
extern CProxy_NearestNeighborSearch nnsArray;
extern CProxy_Interpolate interpolateArray;
extern CProxy_Evaluate evaluateArray;
extern CProxy_DBInterface DBArray;

FineScaleModel::FineScaleModel()
{}

FineScaleModel::FineScaleModel(int state_size, bool use_adaptive_sampling, int nnsIndex, int interpIndex, int dbIndex, int evalIndex) : stateSize(state_size), nnsIndex(nnsIndex), interpIndex(interpIndex), dbIndex(dbIndex), evalIndex(evalIndex) 
{
  // Ordering for a 2D array chare is x, y 
/*
  printf("FineScaleModel created on PE %d Index %d %d state size = %d\n", 
      CkMyPe(), thisIndex.x, thisIndex.y, state_size);
*/

  useAdaptiveSampling = (use_adaptive_sampling) ? 1 : 0;

  ConstitutiveGlobal cm_global;

  // Construct the fine-scale plasticity model
  double D_0 = 1.e-2;
  double m = 1./20.;
  double g = 2.e-3; // (Mbar)
  Plasticity* plasticity_model = (Plasticity*)(new Taylor(D_0, m, g));

  // Construct the approximate nearest neighbors search
  int point_dimension = plasticity_model->pointDimension();

#ifndef NNS_AS_CHARE
#ifdef FLANN

         int n_trees = 1;         // input this from somewhere
         int n_checks = 20;       // input this from somewhere

         ann = (ApproxNearestNeighbors*)(new ApproxNearestNeighborsFLANN(point_dimension,
                                                                         n_trees,
                                                                         n_checks));
#else
         std::string mtreeDirectoryName = ".";

         ann = (ApproxNearestNeighbors*)(new ApproxNearestNeighborsMTree(point_dimension,
                                                                         "kriging_model_database",
                                                                         mtreeDirectoryName,
                                                                         &(std::cout),
                                                                         false));
#endif
#endif

   modelDB = new ModelDB_HashMap(); 
  // Construct the equation of state
  EOS* eos_model;
  /* 
   From Table 1 (converted from GPa to Mbar) in P. J. Maudlin et al.,
   "On the modeling of the Taylor cylinder impact test for orthotropic
   textured materials: experiments and simulations", Inter. J.
   Plasticity 15 (1999), pp. 139-166.
  */
  double k1 = 1.968;  // Mbar
  double k2 = 2.598;  // Mbar
  double k3 = 2.566;  // Mbar
  double Gamma = 1.60;  // dimensionless
  eos_model = (EOS*)(new MieGruneisen(k1, k2, k3, Gamma));

  // Construct the constitutive model
  double bulk_modulus = 1.94; // Tantallum (Mbar)
  double shear_modulus = 6.9e-1; // Tantallum (Mbar)
 
  // Make empty Tensor
  Tensor2Gen L;

  // ann should be NULL
  cm = (Constitutive*)(new ElastoViscoPlasticity(cm_global, 
//#ifndef NNS_AS_CHARE
                                                 ann, 
//#endif
                                                 L, bulk_modulus, shear_modulus, eos_model,
                                                 plasticity_model, use_adaptive_sampling, stateSize));
  em = (ElastoViscoPlasticity*) cm;
}

FineScaleModel::FineScaleModel(CkMigrateMessage *msg)
{
  printf("Migrate Fine Scale chare %d %d\n", thisIndex.x, thisIndex.y);
}

FineScaleModel::~FineScaleModel()
{

}

void FineScaleModel::pup(PUP::er &p)
{
  CBase_FineScaleModel::pup(p);
  p|newPt;
  p|currentIter;
  p|nbrCount;
  p|nbrData;
  p|useAdaptiveSampling;
  p|nnsIndex;
  p|interpIndex;
  p|dbIndex;
  p|stateSize;
}

// Execute advance for element and return results to calling coarse model
void FineScaleModel::advance(const double delta_t, const Tensor2Gen& L_new, const double volume_change, int ssize, char* state)
{
  //ConstitutiveData cm_data = cm->advance(delta_t, L_new, volume_change, state);

   void* plasticity_model_state = em->setState(state);

   em->m_D_new = sym(L_new);
   em->m_W_new = skew(L_new);

   em->m_volume_change = volume_change;

   // Get the deviatoric strain rate at the new time
   Tensor2Sym Dprime_new( dev(em->m_D_new) );

   // Evaluate W^R per equation (29)
   Tensor2Gen WR;
   em->evaluateWR( em->m_Vbar_prime, em->m_Vbar_prime_dot, em->m_Dbar_prime, em->m_W_old,
                   em->m_Wbar, em->m_R, em->a(em->m_J), WR );
   
   // Update the rotation per equation (32)
   Tensor2Gen R;
   em->updateR( em->m_R, WR, delta_t, R );
   
   // Update the stretch determinant per equation (31)
   double J;
   em->updateJ( em->m_J, em->m_D_old, delta_t, J);
   
   Tensor2Sym Vbar_prime;
   // Use the current value of the stretch deviator to 
   // estimate the new time value
   Vbar_prime = em->m_Vbar_prime;
   
   // Update the stretch deviator per equation (30)
   Tensor2Sym Dbar_prime;
   Tensor2Gen Wbar_new;
   updateVbar_prime( em->m_Vbar_prime, Dprime_new, R, em->a(J), delta_t, Vbar_prime, Dbar_prime, Wbar_new );

   // Advance the fine-scale model
   em->m_plasticity_model->advance(delta_t, plasticity_model_state);

   // Update the internal state in preparation for the next call
   em->m_Vbar_prime_dot = Vbar_prime - em->m_Vbar_prime;
   em->m_Vbar_prime_dot /= delta_t;

   em->m_Vbar_prime = Vbar_prime;
   em->m_J = J;
   em->m_R = R;
   em->m_Dbar_prime = Dbar_prime;
   em->m_Wbar = Wbar_new;
   em->m_D_old = em->m_D_new;
   em->m_W_old = em->m_W_new;

   ConstitutiveData return_data;
   return_data.sigma_prime = em->stressDeviator();
   em->getModelInfo(return_data.num_models, return_data.num_point_value_pairs);
   return_data.num_Newton_iters = em->numNewtonIterations();

   em->getState(state);

  int num_samples = 0;
  int num_interpolations = 0;
  if (cm->adaptiveSamplingEnabled())
  {
    num_samples = cm->getNumSamples();
    num_interpolations = cm->getNumSuccessfulInterpolations();
  }
  //coarseScaleArray(thisIndex.x).receiveAdvanceResults(thisIndex.y, cm_data, ssize, state, num_samples, num_interpolations);
  coarseScaleArray(thisIndex.x).receiveAdvanceResults(thisIndex.y, return_data, ssize, state, num_samples, num_interpolations);
}

void FineScaleModel::updateVbar_prime( const Tensor2Sym& Vbar_prime_old,
                                       const Tensor2Sym& Dprime_new,
                                       const Tensor2Gen& R_new,
                                       const double      a_new,
                                       const double      delta_t,
                                       Tensor2Sym&       Vbar_prime_new,
                                       Tensor2Sym&       Dbar_prime_new,
                                       Tensor2Gen&       Wbar_new )
{
/*
 *  This function solves the Vbar_prime equation using the Dogleg Newton algorithm
 *  described in Jorge Nocedal and Stephen J. Wright, "Numerical Optimzation", Second Edition,
 *  Springer Series in Operation Research and Financial Engineering, Springer 2006.
 */

   double eta = 0.24;  // Tolerance for the actual merit function decrease relative
                       // to the prediction given by the local model.  A larger value
   // corresponds to a tighter tolerance, but the maximum value
   // allowed is 0.25;
   int max_iter = 40;  // Maximum number of iterations allowed, including those due only
                       // to trust region modifications
   double tol = 1.e-4; // Requested merit function tolerance.  The merit function is
                       // one-half the square of the residual norm
   
  // end of options (need to move this elsewhere)
  
  // For now, we assume the fine scale model spin is always zero, but we nevertheless
  // include an "update" stub here in case we want to add a non-zero spin in the future.
  Wbar_new = Tensor2Gen(0);

  assert(eta < 0.25);

   // Quantities saved for convergence failure post mortem
   double* diagnostic_merit_value = new double[max_iter+1];
   double* diagnostic_Delta = new double[max_iter+1];
   double* diagnostic_rho = new double[max_iter+1];
   double* diagnostic_p_norm = new double[max_iter+1];
   double* diagnostic_error = new double[max_iter+1];
   int*    diagnostic_model = new int[max_iter+1];

   double tol2 = 2. * tol;

   Tensor4LSym Dbar_prime_deriv;
   evaluateFineScaleModel( em->tauBarPrime(a_new, Vbar_prime_new), Dbar_prime_new, Dbar_prime_deriv );

   double saved_estimate = em->m_error_estimate;
   int current_model = em->m_hint;
   diagnostic_model[0] = em->m_hint;
   diagnostic_rho[0] = 0.;
   diagnostic_error[0] = em->m_error_estimate;

   Tensor2Sym r;
   em->computeResidual( Vbar_prime_new, Vbar_prime_old, Dprime_new, Dbar_prime_new,
                         R_new, a_new * delta_t, r );
   double rnorm2 = em->norm2(r);

   diagnostic_merit_value[0] = 0.5 * rnorm2;
   diagnostic_Delta[0] = em->m_Delta;

   Tensor4LSym jacobian;
   
   em->m_num_iters = 0;
   bool converged = rnorm2 < tol2;
   bool new_local_model = true;

   double Delta_saved = em->m_Delta;

   while ( !converged && em->m_num_iters < max_iter ) {

      if ( new_local_model ) {
         em->computeJacobian( Dbar_prime_deriv, a_new, delta_t, jacobian );
      }

      // Calculate dogleg step
      Tensor2Sym p;
      em->getDoglegStep(r, jacobian, em->m_Delta, p);

      Tensor2Sym proposed_solution = Vbar_prime_new + p;

      diagnostic_p_norm[em->m_num_iters] = norm(p);

      evaluateFineScaleModel( em->tauBarPrime(a_new, proposed_solution), Dbar_prime_new, Dbar_prime_deriv );

      diagnostic_model[em->m_num_iters+1] = em->m_hint;
      diagnostic_error[em->m_num_iters+1] = em->m_error_estimate;

      Tensor2Sym proposed_r;
      em->computeResidual( proposed_solution, Vbar_prime_old, Dprime_new, Dbar_prime_new,
                       R_new, a_new * delta_t, proposed_r );
      double rnorm2_proposed = em->norm2(proposed_r);

      // Evaluate the local linear model
      Tensor2Sym lm_proposed = r + jacobian * p;
      double lm_proposed_norm2 = em->norm2(lm_proposed);

      // Compute the ratio of the actual merit function reduction to the predicted reduction
      // to estimate the accuracy of the local model
      double rho_numerator = rnorm2 - rnorm2_proposed;
      double rho_denominator = rnorm2 - lm_proposed_norm2;

      if( rho_denominator < 0. ) {
         cout << "Negative rho denominator" << endl;
         exit(1);
      }

      double rho = rho_numerator / rho_denominator;

      // Evaluate the effectiveness of the current trust region
      if ( rho < 0.25 ) {  // Local model is not sufficiently predictive; shrink the trust region
         em->m_Delta *= 0.25;
      }
      else {
         if ( rho > 0.75 ) { // Local model works well; try increasing the trust region
            em->m_Delta *= 2.;
            if ( em->m_Delta > em->m_Delta_max ) em->m_Delta = em->m_Delta_max;
         }
         else {
            // Keep the current trust region
         }
      }
     
      if ( rho > eta ) {  // Norm decrease was satisfactory; accept the step 
         Vbar_prime_new = proposed_solution;
         r = proposed_r;
         rnorm2 = rnorm2_proposed;

         converged = rnorm2 < tol2;
         new_local_model = true;
         current_model = em->m_hint;

         saved_estimate = em->m_error_estimate;
         Delta_saved = em->m_Delta;
      }
      else {  // Step was rejected; keep the current solution, residual and Jacobian
         new_local_model = false;
      }

      em->m_num_iters++;

      diagnostic_rho[em->m_num_iters] = rho;
      diagnostic_Delta[em->m_num_iters] = em->m_Delta;
      diagnostic_merit_value[em->m_num_iters] = 0.5 * rnorm2;
   }

   if ( em->m_num_iters >= max_iter && !converged ) {
      cout << "ElastoViscoPlasticity::updateVbar_primeTR() failed to converge" << endl;

      for (int n=0; n<em->m_num_iters; ++n) {
         cout << "merit[" << n << "] = " << diagnostic_merit_value[n] << ", Delta = "
              << diagnostic_Delta[n] << ", model = " << diagnostic_model[n]
              << ", rho = " << diagnostic_rho[n] << ", p_norm = " << diagnostic_p_norm[n]
              << ", error estimate = " << diagnostic_error[n] << endl;
      }

      delete [] diagnostic_error;
      delete [] diagnostic_p_norm;
      delete [] diagnostic_rho;
      delete [] diagnostic_model;
      delete [] diagnostic_Delta;
      delete [] diagnostic_merit_value;
      exit(1);
   }

   delete [] diagnostic_error;
   delete [] diagnostic_p_norm;
   delete [] diagnostic_rho;
   delete [] diagnostic_model;
   delete [] diagnostic_Delta;
   delete [] diagnostic_merit_value;

}

void FineScaleModel::evaluateFineScaleModel( const Tensor2Sym& tau_bar_prime,
                                             Tensor2Sym&       Dbar_prime,
                                             Tensor4LSym&      Dbar_prime_deriv )
{
   if ( cm->adaptiveSamplingEnabled() ) {

      std::vector<double> point(em->m_plasticity_model->pointDimension());
      std::vector<double> value(em->m_plasticity_model->valueAndDerivativeDimension());

      em->m_plasticity_model->packInputVector( tau_bar_prime, point );

      sample( *(em->m_plasticity_model), point, value );

      em->m_plasticity_model->unpackOutputVector( value, Dbar_prime, Dbar_prime_deriv );
   }
   else {

#ifdef EVAL_AS_CHARE
      evaln_message* evaln_msg;
      evaluateArray(evalIndex).evalNative(tau_bar_prime, CkCallbackResumeThread((void*&)evaln_msg));
      evaln_msg->copyTo(Dbar_prime, Dbar_prime_deriv); 
      delete evaln_msg;
#else
      em->m_plasticity_model->evaluateNative( tau_bar_prime, Dbar_prime, Dbar_prime_deriv );
#endif

   }
}

void
FineScaleModel::sample( const FineScale&           fine_scale_model,
                        const std::vector<double>& point,
                        std::vector<double>&       value )
{
   std::vector<bool> interpolateFlags(cm->m_sampler->m_interp->InterpolationDataBase::NUMBER_FLAGS);

   cm->m_sampler->setVerbose(true);

   sample( value,
            point,
            cm->m_hint,
            interpolateFlags,
            fine_scale_model,
            cm->m_error_estimate );
}

void
FineScaleModel::sample( std::vector<double>&       value,
                        const std::vector<double>& point,
                        int&                       hint,
                        std::vector<bool>&         flags,
                        const FineScale&           fine_scale_model,
                        double&                    error_estimate)
{
   // The interpolation database takes pointers to doubles rather
   // than stl vectors, so we need to make temporaries for the
   // query point and returned value.
   
   int point_length = point.size();

   double* local_point = new double[point_length];
   for (int i=0; i<point_length; ++i) {
      local_point[i] = point[i] / cm->m_sampler->m_pointScaling[i];
   }

   double local_point_norm = cm->m_sampler->pointL2Norm(local_point);
   if ( local_point_norm > cm->m_sampler->m_point_norm_max ) {
      cm->m_sampler->m_point_norm_max = local_point_norm;
   }
   cm->m_sampler->m_point_norm_sum += local_point_norm;

   int value_length = value.size();
   double* local_value = new double[value_length];

   bool interpolationSuccess =
      //cm->m_sampler->m_interp->interpolate(local_value,
      interpolate(local_value,
                            hint,
                            local_point,
                            flags,
                            error_estimate);

   if (interpolationSuccess == false) {

#ifdef EVAL_AS_CHARE
      eval_message* eval_msg;
      evaluateArray(evalIndex).eval(point, CkCallbackResumeThread((void*&)eval_msg));
      eval_msg->copyTo(value);
      delete eval_msg;

#else
      fine_scale_model.evaluate(point, value);
#endif

      if (cm->m_sampler->m_verbose) {
         //         cout << "Interpolation failed: Adding key" << endl;
      }

      for (int i=0; i<value_length; ++i) {
         local_value[i] = value[i] / cm->m_sampler->m_valueScaling[i];
      }

      //cm->m_sampler->m_interp->insert(hint,
      insert(hint,
                       local_point,
                       local_value,
                       &(local_value[cm->m_sampler->m_valueDimension]),
                       flags);

      cm->m_sampler->m_num_fine_scale_evaluations++;

      hint = -1;
      error_estimate = 0.;

   }
   else {

      if (cm->m_sampler->m_verbose) {
         //         cout << "Interpolation succeeded" << endl;
      }

      interpolationSuccess =
      //cm->m_sampler->m_interp->interpolate(local_value,
      interpolate(local_value,
                            local_value + point_length,
                            hint,
                            local_point,
                            flags,
                            error_estimate);

      for (int i=0; i<value_length; ++i) {
         value[i] = local_value[i] * cm->m_sampler->m_valueScaling[i];
      }

      //      verifyInterpolationAccuracy(point, value, fineScaleModel);
    

      cm->m_sampler->m_num_successful_interpolations++;
   }

   cm->m_sampler->m_num_samples++;

   assert(cm->m_sampler->m_num_samples = cm->m_sampler->m_num_fine_scale_evaluations + cm->m_sampler->m_num_successful_interpolations);

   double local_value_norm = cm->m_sampler->valueL2Norm(local_value);
   if ( local_value_norm > cm->m_sampler->m_value_norm_max ) {
      cm->m_sampler->m_value_norm_max = local_value_norm;
   }
   cm->m_sampler->m_value_norm_sum += local_value_norm;

   delete [] local_value;
   delete [] local_point;
}

bool FineScaleModel::interpolate(double            * value,
                                 int               & hint,
                                 const double      * point,
                                 std::vector<bool> & flags,
                                 double            & error_estimate)
{
      //
      // make sure there is enough space in flags 
      //
       assert(flags.size() >= InterpolationDataBase::NUMBER_FLAGS);
      
      //
      // initialize flags container
      //

      std::fill(flags.begin(),
                flags.end(),
                false);

      //
      // shortcuts to frequently accesses data
      //
      
      const int pointDimension = cm->m_sampler->m_interp->getPointDimension();
      const int valueDimension = cm->m_sampler->m_interp->getValueDimension();
     
      //
      // instatiate point object from point data
      //
    
      const ResponsePoint queryPoint(pointDimension,
                                     point);
      
      //
      // use hint to get the most recently used model
      //
    
      if (cm->m_sampler->m_interp->_useHint && hint != -1) {

#ifdef NNS_AS_CHARE
         uint128_message* key_msg;
         nnsArray(nnsIndex).getKey(hint, CkCallbackResumeThread((void*&)key_msg));
         uint128_t model_key = key_msg->value;
         delete key_msg;
#else
         uint128_t model_key = cm->m_sampler->m_interp->_ann.getKey(hint);
#endif

         if (model_key == uint128_t_undefined) {

            flags[InterpolationDataBase::LOST_HINT_FLAG] = true;

         } else {

            InterpolationModelPtr hintKrigingModel = cm->m_sampler->m_interp->_modelDB->extract(model_key);

            //
            // check if can interpolate; need a valid model for this
            //
     
            if (hintKrigingModel->isValid() == true) {  
      
               //
               // check the distance between hintKrigingModel and point
               //
               const Point modelCenter = getModelCenterMass(*hintKrigingModel);
               const Vector pointRelativePosition = queryPoint - modelCenter;
               const double distanceSqr = krigalg::dot(pointRelativePosition,
                                                       pointRelativePosition);

               if (distanceSqr >
                   cm->m_sampler->m_interp->_maxQueryPointModelDistance*cm->m_sampler->m_interp->_maxQueryPointModelDistance) {
                  flags[InterpolationDataBase::LOST_HINT_FLAG] = true;
               } else {

                  const bool hintModelSuccess =
                     checkErrorAndInterpolate(value,
                                              hintKrigingModel,
                                              queryPoint,
                                              valueDimension,
                                              cm->m_sampler->m_interp->_tolerance,
                                              cm->m_sampler->m_interp->_meanErrorFactor,
                                              error_estimate );

                  if (hintModelSuccess == true) {
                     flags[InterpolationDataBase::USED_HINT_FLAG] = true;
                     return true;
                  }

               }

            }
         }

      }

      //
      // A kriging model based on hint did not produce a valid interpolation.
      // Find closest kriging model.
      //

      if (cm->m_sampler->m_interp->_maxNumberSearchModels == 1) {

        const std::pair<int, InterpolationModelPtr>
          closestKrigingModelData = findClosestCoKrigingModel(queryPoint,
#ifndef NNS_AS_CHARE
                                                              cm->m_sampler->m_interp->_ann,
#endif
                                                              cm->m_sampler->m_interp->_modelFactory,
                                                              cm->m_sampler->m_interp->_modelDB,
                                                              cm->m_sampler->m_interp->_maxQueryPointModelDistance);

        InterpolationModelPtr closestKrigingModel =
          closestKrigingModelData.second;
        hint = closestKrigingModelData.first;

        //
        // if no kriging model is available return
        //

        if (hint == id_undefined || closestKrigingModel->isValid() == false) {

          return false;

        }

        //
        // estimate error for all values from the kriging models
        //
        
        const bool interpolationSuccess =
          checkErrorAndInterpolate(value,
                                   closestKrigingModel,
                                   queryPoint,
                                   valueDimension,
                                   cm->m_sampler->m_interp->_tolerance,
                                   cm->m_sampler->m_interp->_meanErrorFactor,
                                   error_estimate);

        return interpolationSuccess;

      } else {

        //
        // search more than one kriging model
        //

        bool canInterpolateFlag;

        const std::pair<int, InterpolationModelPtr>
          bestKrigingModelData = findBestCoKrigingModel(canInterpolateFlag,
                                                        queryPoint,
#ifndef NNS_AS_CHARE
                                                        cm->m_sampler->m_interp->_ann,
#endif
                                                        cm->m_sampler->m_interp->_modelDB,
                                                        cm->m_sampler->m_interp->_modelFactory,
                                                        cm->m_sampler->m_interp->_tolerance,
                                                        cm->m_sampler->m_interp->_meanErrorFactor,
                                                        cm->m_sampler->m_interp->_maxQueryPointModelDistance,
                                                        cm->m_sampler->m_interp->_maxNumberSearchModels,
                                                        cm->m_sampler->m_interp->_maxKrigingModelSize,
                                                        valueDimension);

        InterpolationModelPtr bestKrigingModel =
          bestKrigingModelData.second;
        hint = bestKrigingModelData.first;

        //
        // if no kriging model is available return
        //

        if (canInterpolateFlag == false) {

          return false;

        }

        //
        // interpolate using the best kriging model available
        //
        
        return checkErrorAndInterpolate(value,
                                        bestKrigingModel,
                                        queryPoint,
                                        valueDimension,
                                        cm->m_sampler->m_interp->_tolerance,
                                        cm->m_sampler->m_interp->_meanErrorFactor,
                                        error_estimate);

        return true;

      }

      //
      // this should never be reached
      //
      
      assert(false);
      
      return false;
      
}

bool FineScaleModel::interpolate(double            * value,
                                 double            * gradient,
                                 int               & hint,
                                 const double      * point,
                                 std::vector<bool> & flags,
                                 double            & error_estimate)
{

      //
      // make sure there is enough space in flags 
      //

      assert(flags.size() >= InterpolationDataBase::NUMBER_FLAGS);

      //
      // initialize flags container
      //

      std::fill(flags.begin(),
                flags.end(),
                false);
      //
      // shortcuts to frequently accesses data
      //

      const int pointDimension = cm->m_sampler->m_interp->getPointDimension();
      const int valueDimension = cm->m_sampler->m_interp->getValueDimension();

      //
      // instatiate point object from point data
      //

      const ResponsePoint queryPoint(pointDimension,
                                     point);

      //
      // try hint (if valid)
      //

      if (cm->m_sampler->m_interp->_useHint && hint != -1) {

#ifdef NNS_AS_CHARE
         uint128_message* key_msg;
         nnsArray(nnsIndex).getKey(hint, CkCallbackResumeThread((void*&)key_msg));
         uint128_t model_key = key_msg->value;
         delete key_msg;
#else
         uint128_t model_key = cm->m_sampler->m_interp->_ann.getKey(hint);
#endif

        if (model_key == uint128_t_undefined) {

          flags[InterpolationDataBase::LOST_HINT_FLAG] = true;

        } else {

          const InterpolationModelPtr hintKrigingModel = cm->m_sampler->m_interp->_modelDB->extract(model_key);

          //
          // check the distance between hintKrigingModel and point
          //

          const Point modelCenter = getModelCenterMass(*hintKrigingModel);
          const Vector pointRelativePosition = queryPoint - modelCenter;
          const double distanceSqr = krigalg::dot(pointRelativePosition,
                                                    pointRelativePosition);

          if (distanceSqr >
                cm->m_sampler->m_interp->_maxQueryPointModelDistance*cm->m_sampler->m_interp->_maxQueryPointModelDistance) {
            flags[InterpolationDataBase::LOST_HINT_FLAG] = true;
          } else {

            const bool hintModelSuccess =
              checkErrorAndInterpolate(value,
                                       gradient,
                                       hintKrigingModel,
                                       queryPoint,
                                       pointDimension,
                                       valueDimension,
                                       cm->m_sampler->m_interp->_tolerance,
                                       cm->m_sampler->m_interp->_meanErrorFactor,
                                       error_estimate);

            if (hintModelSuccess == true) {
              flags[InterpolationDataBase::USED_HINT_FLAG] = true;
              return true;
            }

          }

        }

      }

      //
      // find closest kriging model
      //
      
      if (cm->m_sampler->m_interp->_maxNumberSearchModels == 1) {

        const std::pair<int, InterpolationModelPtr>
          closestKrigingModelData = findClosestCoKrigingModel(queryPoint,
#ifndef NNS_AS_CHARE
                                                              cm->m_sampler->m_interp->_ann,
#endif
                                                              cm->m_sampler->m_interp->_modelFactory,
                                                              cm->m_sampler->m_interp->_modelDB,
                                                              cm->m_sampler->m_interp->_maxQueryPointModelDistance);

        InterpolationModelPtr closestKrigingModel =
          closestKrigingModelData.second;
        hint = closestKrigingModelData.first;

        //
        // if no kriging model is available return
        //

        if (hint == id_undefined) {

          return false;

        }

        //
        // check error and interpolate
        //

        const bool interpolationSuccess =
          checkErrorAndInterpolate(value,
                                   gradient,
                                   closestKrigingModel,
                                   queryPoint,
                                   pointDimension,
                                   valueDimension,
                                   cm->m_sampler->m_interp->_tolerance,
                                   cm->m_sampler->m_interp->_meanErrorFactor,
                                   error_estimate);

        return interpolationSuccess;

      } else {

        //
        // search more than one kriging model
        //

        bool canInterpolateFlag;

        const std::pair<int, InterpolationModelPtr>
          bestKrigingModelData = findBestCoKrigingModel(canInterpolateFlag,
                                                        queryPoint,
#ifndef NNS_AS_CHARE
                                                        cm->m_sampler->m_interp->_ann,
#endif
                                                        cm->m_sampler->m_interp->_modelDB,
                                                        cm->m_sampler->m_interp->_modelFactory,
                                                        cm->m_sampler->m_interp->_tolerance,
                                                        cm->m_sampler->m_interp->_meanErrorFactor,
                                                        cm->m_sampler->m_interp->_maxQueryPointModelDistance,
                                                        cm->m_sampler->m_interp->_maxNumberSearchModels,
                                                        cm->m_sampler->m_interp->_maxKrigingModelSize,
                                                        valueDimension);

        InterpolationModelPtr bestKrigingModel =
          bestKrigingModelData.second;
        hint = bestKrigingModelData.first;

        //
        // if no kriging model is available return
        //

        if (canInterpolateFlag == false) {


          return false;

        }

        //
        // interpolate using the best kriging model available
        //

        interpolate(value,
                     gradient,
                     bestKrigingModel,
                     queryPoint,
                     pointDimension,
                     valueDimension);

        return true;

      }

      //
      // this should never be reached
      //

      assert(false);

      return false;

}

void FineScaleModel::interpolate(
                    double                                * value,
                    double                                * gradient,
                    InterpolationModelPtr                   krigingModel,
                    const ResponsePoint                   & queryPoint,
                    int                                     pointDimension,
                    int                                     valueDimension)
{
        //
        // firewalls
        //

        assert(krigingModel->hasGradient() == true);

        for (int iValue = 0; iValue < valueDimension; ++iValue) {

          //
          // compute the value
          //

          const Value valueEstimate =
            krigingModel->interpolate(iValue,
                                      queryPoint);

          //
          // put the value of the function into value (valueEstimate
          // contains the value of the function follwed by the gradient;
          //

          value[iValue] = valueEstimate[0];

          //
          // store gradient data
          //

          for (int i = 0; i < pointDimension; ++i)
            gradient[i*valueDimension + iValue] = valueEstimate[1 + i];

        }

        return;
}

void FineScaleModel::insert(int               & hint,
              const double      * point,
              const double      * value,
              const double      * gradient,
              std::vector<bool> & flags)
{
      // 
      // make sure there is enough space in flags 
      //

      assert(flags.size() >= InterpolationDataBase::NUMBER_FLAGS);

      //
      // initialize flags container
      //

      std::fill(flags.begin(),
                flags.end(),
                false);

      //
      // update number of point/value pairs
      //
      
      ++cm->m_sampler->m_interp->_numberPointValuePairs;

      //
      // shortcuts to frequently accessed data
      //

      const int pointDimension = cm->m_sampler->m_interp->getPointDimension();
      const int valueDimension = cm->m_sampler->m_interp->getValueDimension();

      //
      // check if the list of kriging models in non-empty
      //

      if (hint == id_undefined) {

        //
        // create and add new model
        //

        addNewModel(
                     cm->m_sampler->m_interp->_modelDB,
#ifndef NNS_AS_CHARE
                     cm->m_sampler->m_interp->_ann,
#endif
                     cm->m_sampler->m_interp->_modelFactory,
                     hint,
                     point,
                     value,
                     gradient,
                     pointDimension,
                     valueDimension);

        //
        // update number of kriging models
        //

        ++cm->m_sampler->m_interp->_numberKrigingModels;

      } else {

        //
        // get a handle to the right kriging model
        //

#ifdef NNS_AS_CHARE
         uint128_message* key_msg;
         nnsArray(nnsIndex).getKey(hint, CkCallbackResumeThread((void*&)key_msg));
         uint128_t model_key = key_msg->value;
         delete key_msg;
#else
        uint128_t model_key = cm->m_sampler->m_interp->_ann.getKey(hint);
#endif

        InterpolationModelPtr krigingModel = cm->m_sampler->m_interp->_modelDB->extract(model_key);

        //
        // check the size of the model; if the next point would put
        // the model above _maxKrigingModelSize start a new model;
        // otherwise just add point/value pair
        //

        if (krigingModel->getNumberPoints() == cm->m_sampler->m_interp->_maxKrigingModelSize) {

           addNewModel(
                       cm->m_sampler->m_interp->_modelDB,
#ifndef NNS_AS_CHARE
                       cm->m_sampler->m_interp->_ann,
#endif
                       cm->m_sampler->m_interp->_modelFactory,
                       hint,
                       point,
                       value,
                       gradient,
                       pointDimension,
                       valueDimension);

          //
          // update number of kriging models
          //

          ++cm->m_sampler->m_interp->_numberKrigingModels;

          //
          // record event
          //

          flags[InterpolationDataBase::MODEL_SIZE_LIMIT_FLAG] = true;

        } else {

          //
          // copy value and gradient data
          //

          std::vector<Value> pointValue;


          if (krigingModel->hasGradient() == true)
            pointValue = copyValueData(value,
                                       gradient,
                                       pointDimension,
                                       valueDimension);
          else
            pointValue = copyValueData(value,
                                       pointDimension,
                                       valueDimension);

          //
          // copy point data into a Point object
          //

          const Point pointObject(pointDimension,
                                  point);

          //
          // insert point/value data into the model
          //

          const bool addPointSuccess = krigingModel->addPoint(pointObject,
                                                              pointValue);

          if (addPointSuccess == true) {

            //
            // compute the center of mass for the updated model
            //

            const Point centerMass = getModelCenterMass(*krigingModel);

            //
            // remove old kriging model from the database
            //

#ifdef NNS_AS_CHARE
            nnsArray(nnsIndex).remove(hint);
#else
            cm->m_sampler->m_interp->_ann.remove(hint);
#endif

            cm->m_sampler->m_interp->_modelDB->erase(model_key);

            //
            // insert updated kriging model into database
            //

            const ResponsePoint centerMassRP(pointDimension,
                                             &(centerMass[0]));

            // Create a key string corresponding to the center of mass point
            // and insert it into the approximate nearest neighbor database
            
            std::vector<double> point_data;
            point_data.resize(centerMassRP.size());
            for (int i=0; i<centerMassRP.size(); ++i) {
               point_data[i] = centerMassRP[i];
            }

            uint128_t new_model_key = getKeyHash(centerMassRP);
#ifdef NNS_AS_CHARE
            int_message* insert_msg;
            nnsArray(nnsIndex).insert(point_data, new_model_key, CkCallbackResumeThread((void*&)insert_msg));
            hint = insert_msg->value;
            delete insert_msg;
#else
            hint = cm->m_sampler->m_interp->_ann.insert(point_data, new_model_key);
#endif

            // Insert the interpolation model into the interpolation model database

            cm->m_sampler->m_interp->_modelDB->insert( new_model_key, krigingModel, (ResponsePoint *) &centerMassRP );

          } else {

            //
            // point insertion failed-add new model
            //

                         addNewModel(
                         cm->m_sampler->m_interp->_modelDB,
#ifndef NNS_AS_CHARE
                         cm->m_sampler->m_interp->_ann,
#endif
                         cm->m_sampler->m_interp->_modelFactory,
                         hint,
                         point,
                         value,
                         gradient,
                         pointDimension,
                         valueDimension);

            //
            // update number of kriging models
            //

            ++cm->m_sampler->m_interp->_numberKrigingModels;

            //
            // record event
            //

            flags[InterpolationDataBase::MODEL_INSERT_LIMIT_FLAG] = true;

          }

        }

      }

      return;

}

Point FineScaleModel::getModelCenterMass(const InterpolationModel & krigingModel)
{
    const std::vector<Point> & points = krigingModel.getPoints();
    assert(points.empty() == false);
    Point centerMass(points.front().size(),
                         0.0);

    std::vector<Point>::const_iterator pointsIter;
        std::vector<Point>::const_iterator pointsEnd = points.end();

        for (pointsIter  = points.begin();
             pointsIter != pointsEnd;
             ++pointsIter) {
            const Point & point = *pointsIter;
            centerMass += static_cast<Vector>(point);
        }

        mtl::scale(centerMass,
                   1.0/points.size());

        return centerMass;
}

std::pair<int, InterpolationModelPtr>
FineScaleModel::findClosestCoKrigingModel(
                const ResponsePoint        & point,
#ifndef NNS_AS_CHARE
                ApproxNearestNeighbors     & ann,
#endif
                krigalg::InterpolationModelFactoryPointer modelFactory,
                ModelDatabase * modelDB,
                double                       maxQueryPointModelDistance)
{
        //
        // query tree for the closest model
        //

        int k_neighbors = 1;

        std::vector<double> x;
        x.resize(point.size());
        for (int i=0; i<point.size(); ++i) {
           x[i] = point[i];
        }
#ifdef NNS_AS_CHARE
        std::vector<int> ids;
        std::vector<uint128_t> keys;
        std::vector<double> dists;

        knn_message* knn_msg;
        nnsArray(nnsIndex).knn(x, k_neighbors, CkCallbackResumeThread((void*&)knn_msg));
        knn_msg->copyTo(ids, keys, dists);
        delete knn_msg;
#else
        std::vector<int> ids(k_neighbors);
        std::vector<uint128_t> keys(k_neighbors);
        std::vector<double> dists(k_neighbors);

        ann.knn(x, k_neighbors, ids, keys, dists);
#endif

        bool found_neighbors = (ids.size() == k_neighbors);

        InterpolationModelPtr closestKrigingModel;
        int closestKrigingModelId;

        //
        // If a neighbor is found, compare the distance between the query point and the model
        // with the value of maxQueryPointModelDistance
        //
        
        if ( found_neighbors &&
             dists[0] <= maxQueryPointModelDistance ) {

           //
           // get handle to the located object
           //
           
           uint128_t model_key = keys[0];

           closestKrigingModel = modelDB->extract(model_key);

           closestKrigingModelId = ids[0];
        }
        else {
           closestKrigingModelId = -1;
        }

        return std::make_pair(closestKrigingModelId,
                              closestKrigingModel);

}

std::pair<int, InterpolationModelPtr>
FineScaleModel::findBestCoKrigingModel(
                bool &                       canInterpolateFlag,
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
                int                          valueDimension)
{
        printf("in findBest\n");
exit(1);
        canInterpolateFlag = false;

        //
        // query the tree for the maxNumberSearchModels closest models
        //

        int k_neighbors = maxNumberSearchModels;

        std::vector<double> x;
        x.resize(point.size());
        for (int i=0; i<point.size(); ++i) {
           x[i] = point[i];
        }
#ifdef NNS_AS_CHARE
        std::vector<int> ids;
        std::vector<uint128_t> keys;
        std::vector<double> dists;

        //knn_message* knn_msg = new knn_message();
        knn_message* knn_msg;
        nnsArray(nnsIndex).knn(x, k_neighbors, CkCallbackResumeThread((void*&)knn_msg));
        knn_msg->copyTo(ids, keys, dists);
        delete knn_msg;
#else
        std::vector<int> ids(k_neighbors);
        std::vector<uint128_t> keys(k_neighbors);
        std::vector<double> dists(k_neighbors);

        ann.knn(x, k_neighbors, ids, keys, dists);
#endif

        bool found_neighbors = (ids.size() > 0);

        InterpolationModelPtr bestKrigingModel;
        int bestKrigingModelId;

        std::pair<int, InterpolationModelPtr> bestCoKrigingModel;

        if ( found_neighbors ) {

           //
           // iterate through the search results
           //

           double minError = std::numeric_limits<double>::max();

           for (int iter=0; iter<ids.size(); ++iter) {

              uint128_t model_key = keys[iter];

              InterpolationModelPtr krigingModel = modelDB->extract(model_key);
              //
              // skip invalid models
              //

              if (krigingModel->isValid() == false)
                 continue;

              //
              // compute error
              //

              const double errorEstimate = checkError(krigingModel,
                                                      point,
                                                      valueDimension,
                                                      meanErrorFactor);

              //
              // check if a better model encountered
              //

              if (errorEstimate < minError) {

                 bestCoKrigingModel = std::make_pair(ids[iter], krigingModel);

                 minError = errorEstimate;

              }

           }

           //
           // check if the best model satisfies tolerance requirement
           //

           if (minError <= tolerance*tolerance)
              canInterpolateFlag = true;

        }
        else {
          bestCoKrigingModel = make_pair(-1, InterpolationModelPtr());
        }

        return bestCoKrigingModel;
}

std::pair<int, InterpolationModelPtr>
FineScaleModel::findBestCoKrigingModel(
                bool &                       canInterpolateFlag,
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
                int                          valueDimension)
{
        canInterpolateFlag = false;

        //
        // query the tree for the maxNumberSearchModels closest models
        //

        int k_neighbors = maxNumberSearchModels;

        std::vector<double> x;
        x.resize(point.size());
        for (int i=0; i<point.size(); ++i) {
           x[i] = point[i];
        }
#ifdef NNS_AS_CHARE
        std::vector<int> ids;
        std::vector<uint128_t> keys;
        std::vector<double> dists;

        knn_message* knn_msg;
        nnsArray(nnsIndex).knn(x, k_neighbors, CkCallbackResumeThread((void*&)knn_msg));
        knn_msg->copyTo(ids, keys, dists);
        delete knn_msg;
#else
        std::vector<int> ids(k_neighbors);
        std::vector<uint128_t> keys(k_neighbors);
        std::vector<double> dists(k_neighbors);

        ann.knn(x, k_neighbors, ids, keys, dists);
#endif

        bool found_neighbors = (ids.size() > 0);

        if ( found_neighbors ) {

           //
           // iterate through the search results
           //
           
           for (int iter=0; iter<ids.size(); ++iter) {

              uint128_t model_key = keys[iter];

              InterpolationModelPtr krigingModel = modelDB->extract(model_key);

              //
              // skip if invalid model
              //
              
              if (krigingModel->isValid() == false)
                 continue;

              //
              // skip if too far away; as results come back closest-first, we know
              // all of the other results are too far away as well, so can break
              // could probably 
              //

              if (dists[iter] > maxQueryPointModelDistance)
                 break;

              //
              // compute error
              //

              const double errorEstimate = checkError(krigingModel,
                                                      point,
                                                      valueDimension,
                                                      meanErrorFactor);
              if (errorEstimate <= tolerance*tolerance) {
                 canInterpolateFlag = true;

                 return std::make_pair(ids[iter], krigingModel);

              }
           }

           //
           // a model suitable for interpolation has not been
           // found-return the closest model
           //
           
           uint128_t model_key = keys[0];

           InterpolationModelPtr krigingModel = modelDB->extract(model_key);

           return std::make_pair(ids[0], krigingModel);

        }
        else {
           return std::make_pair(-1, InterpolationModelPtr());
        }
}

double
FineScaleModel::checkError(InterpolationModelPtr   krigingModel,
                           const ResponsePoint   & queryPoint,
                           int                     valueDimension,
                           double                 _meanErrorFactor)
{

        return checkErrorInf(krigingModel,
                             queryPoint,
                             valueDimension,
                             _meanErrorFactor);

}

double
FineScaleModel::checkErrorInf(InterpolationModelPtr   krigingModel,
                              const ResponsePoint   & queryPoint,
                              int                     valueDimension,
                              double                 _meanErrorFactor)
{
        double maxError = 0.0;

        //
        // iterate over values
        //

        for (int iValue = 0; iValue < valueDimension; ++iValue) {

          //
          // compute the error estimate
          //

          const double errorEstimate =
            compKrigingError(*krigingModel,
                             queryPoint,
                             iValue,
                             _meanErrorFactor);

          //
          // save error
          //

          maxError = std::max(maxError,
                              errorEstimate);


        }

        return maxError;
}

double
FineScaleModel::compKrigingError(const InterpolationModel & krigingModel,
                                 const Point              & queryPoint,
                                 int                        valueId,
                                 double                     meanErrorFactor)
{
        const int numberPoints = krigingModel.getNumberPoints();

        //
        // firewalls
        //

        assert(numberPoints >= 1);

        //
        // compute min number of points to attempt meaningful
        // interpolation
        //

        const int minNumberPoints = krigingModel.hasGradient() ? 1 :
          2*(krigingModel.getPointDimension() + 1) - 1;

        //
        // compute the error if the kriging model contains a single point;
        // otherwise, simply return the kriging prediction
        //

        if (numberPoints <= minNumberPoints ) {

          //
          // get the kriging estimate at query point
          //

          const Value queryValue = krigingModel.interpolate(valueId,
                                                            queryPoint);

          //
          // get all points in the model; there should really only be
          // ONE if we have gradient information
          //

          const std::vector<Point> & points = krigingModel.getPoints();
          assert( krigingModel.hasGradient() ? (points.size() == 1) : true );

          //
          // get the kriging estimate at the origin of the kriging model
          //

          const Value originValue =
            krigingModel.interpolate(valueId,
                                     points.front());

          //
          // compute the error as the difference between the value at the
          // queryPoint and origin point
          //

          return (queryValue[0] - originValue[0])*
            (queryValue[0] - originValue[0]);


        } else
          return meanErrorFactor*meanErrorFactor*
            krigingModel.getMeanSquaredError(valueId, queryPoint)[0];

        //
        // can never be reached
        //
        
        return 0.0;
}

std::vector<Value>
FineScaleModel::copyValueData(const double * value,
                              const double * gradient,
                              int            pointDimension,
                              int            valueDimension)
{

        //
        // instatiate return object
        //

        std::vector<Value> pointValues;

        //
        // iterate over all values
        //

        for (int iValue = 0; iValue < valueDimension; ++iValue) {

          //
          // copy value data
          //

          Value pointValue(pointDimension + 1);

          pointValue[0] = value[iValue];

          //
          // copy gradient data
          //

          for (int i = 0; i < pointDimension; ++i)
            pointValue[1 + i] = gradient[i*valueDimension + iValue];

          //
          // add pointValue
          //

          pointValues.push_back(pointValue);

        }

        return pointValues;

}

std::vector<Value>
FineScaleModel::copyValueData(const double * value,
                               int            pointDimension,
                               int            valueDimension)
{

        //
        // instatiate return object
        //

        std::vector<Value> pointValues;

        //
        // iterate over all values
        //

        for (int iValue = 0; iValue < valueDimension; ++iValue) {

          //
          // copy value data
          //

          Value pointValue(1);

          pointValue[0] = value[iValue];

          //
          // add pointValue
          //

          pointValues.push_back(pointValue);

        }

        return pointValues;

}

uint128_t 
FineScaleModel::getKeyHash(const ResponsePoint& point)
       {
          int point_size = point.size();

          std::vector<double> data;
          data.resize(point_size);
          for (int i=0; i<point_size; ++i) {
             data[i] = point[i];
          }

          uint128_t hash;
          MurmurHash3_x64_128(&data[0], point_size*sizeof(double)/sizeof(char), MURMUR_SEED, &hash);

          return hash;
       }

bool 
FineScaleModel::checkErrorAndInterpolate(double               * value,
                                         InterpolationModelPtr  krigingModel,
                                         const ResponsePoint  & queryPoint,
                                         int                    valueDimension,
                                         double                _tolerance,
                                         double                _meanErrorFactor,
                                         double               & errorEstimate)
{
        const double toleranceSqr = _tolerance*_tolerance;

        errorEstimate = 0.;

        for (int iValue = 0; iValue < valueDimension; ++iValue) {

          //
          // compute the error estimate
          //

          const double iErrorEstimate =
            compKrigingError(*krigingModel,
                             queryPoint,
                             iValue,
                             _meanErrorFactor);

          errorEstimate = std::max(errorEstimate, sqrt(iErrorEstimate));


          //
          // check the errorEstimate against the tolerance; if the 
          // estimate is greater than tolerance simpoy return failure
          //

          if ( fabs(iErrorEstimate) > toleranceSqr) {

            return false;

          }

          //
          // compute the value
          //

          const Value valueEstimate =
            krigingModel->interpolate(iValue,
                                      queryPoint);

          //
          // put the value of the function into value (valueEstimate
          // contains the value of the function follwed by the gradient;
          // here we are interested only in the value of the function so
          // the gradient is simply discarded).
          //

          value[iValue] = valueEstimate[0];

        }

        return true;

}

bool
FineScaleModel::checkErrorAndInterpolate(double               * value,
                                         double               * gradient,
                                         InterpolationModelPtr  krigingModel,
                                         const ResponsePoint  & queryPoint,
                                         int                    pointDimension,
                                         int                    valueDimension,
                                         double                _tolerance,
                                         double                _meanErrorFactor,
                                         double               & errorEstimate )
{
       //
       // firewall
       //

       assert(krigingModel->hasGradient() == true);

        //
        // estimate error for all values from the kriging models
        //

        const double toleranceSqr = _tolerance*_tolerance;

        errorEstimate = 0.;

        for (int iValue = 0; iValue < valueDimension; ++iValue) {

          //
          // compute the error estimate
          //

const double iErrorEstimate =
            compKrigingError(*krigingModel,
                             queryPoint,
                             iValue,
                             _meanErrorFactor);

          errorEstimate = std::max(errorEstimate, sqrt(iErrorEstimate));

          //
          // check the errorEstimate against the tolerance; if the 
          // estimate is greater than tolerance simpoy return failure
          //

          if ( fabs(iErrorEstimate) > toleranceSqr) {

            return false;

          }

          //
          // compute the value
          //

          const Value valueEstimate =
            krigingModel->interpolate(iValue,
                                      queryPoint);

          //
          // put the value of the function into value (valueEstimate
          // contains the value of the function follwed by the gradient;
          //

          value[iValue] = valueEstimate[0];

          //
          // store gradient data
          //
          
          for (int i = 0; i < pointDimension; ++i)
            gradient[i*valueDimension + iValue] = valueEstimate[1 + i];

        }

        return true;

}

void FineScaleModel::addNewModel(
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
                    int                                      valueDimension)
{  
        //
        // start new kriging model
        //

        InterpolationModelPtr krigingModel = _modelFactory->build();

        //
        // copy value and gradient data
        //

        std::vector<Value> pointValue;


        if (krigingModel->hasGradient() == true)
          pointValue = copyValueData(valueData,
                                     gradientData,
                                     pointDimension,
                                     valueDimension);
        else
          pointValue = copyValueData(valueData,
                                     pointDimension,
                                     valueDimension);

        //
        // copy point data into a Point object
        //

        const ResponsePoint point(pointDimension,
                                  pointData );

        //
        // insert point/value data into the model
        //

        krigingModel->addPoint(point,
                               pointValue);

        // Create a key string corresponding to the new point at
        // which the new interpolation model is centered
        //
        // Insert the model key into the approximate nearest neighbor database
        //

        std::vector<double> point_data;
        point_data.resize(point.size());
        for (int i=0; i<point.size(); ++i) {
           point_data[i] = point[i];
        }

        uint128_t model_key = getKeyHash(point);
#ifdef NNS_AS_CHARE
        int_message* insert_msg;
        nnsArray(nnsIndex).insert(point_data, model_key, CkCallbackResumeThread((void*&)insert_msg));
        objectId = insert_msg->value;
        delete insert_msg;
#else
        objectId = ann.insert(point_data, model_key);
#endif

        // Insert the interpolation model into the interpolation model database
        
        modelDB->insert(model_key,krigingModel, (ResponsePoint *)&point);

        return;

}

