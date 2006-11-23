#ifndef FastSimulation_ShowerDevelopment_HSParameters_H
#define FastSimulation_ShowerDevelopment_HSParameters_H

#include "FWCore/ParameterSet/interface/ParameterSet.h"
/** 
 *  Parameters used in the hadron fast simulation 
 */ 

class HSParameters
{
 public:
  HSParameters(){;}
  HSParameters(const edm::ParameterSet& params);
  ~HSParameters(){;}
    
  // Methods to provide FamosHDShower with constants 
  inline int    getHDlossesOpt     ()  const {return lossesOpt_;}
  inline int    getHDnDepthSteps   ()  const {return nDepthSteps_;}
  inline int    getHDnTRsteps      ()  const {return nTRsteps_;}
  inline double getHDtransParam    ()  const {return transParam_;}
  inline double getHDeSpotSize     ()  const {return eSpotSize_;}
  inline double getHDdepthStep     ()  const {return depthStep_;}
  inline double getHDcriticalEnergy()  const {return criticalEnergy_;}
  inline double getHDmaxTRfactor   ()  const {return maxTRfactor_;}
  inline double getHDbalanceEH     ()  const {return balanceEH_;}



 private:
  //FamosHDshower-related group of parameters
  int    lossesOpt_,  nDepthSteps_,    nTRsteps_;
  double transParam_, eSpotSize_, depthStep_, criticalEnergy_, maxTRfactor_;
  double balanceEH_;
  

};

#endif
