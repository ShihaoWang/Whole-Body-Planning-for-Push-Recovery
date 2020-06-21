#include "CommonHeader.h"
#include "NonlinearOptimizerInfo.h"

static double PositivePosNVel(  const double & PosDiff,
                                const double & InitVelocity,  double & GoalVelocity,
                                const double & VelocityBound, const double & AccBound){
   // Two procedure => accelerate to velocity limit and remains there
   double AccDist = (VelocityBound * VelocityBound - InitVelocity * InitVelocity)/(2.0 * AccBound);
   if(AccDist>PosDiff){ // Does not need to accelerate to limit to reach GoalPos
     double AccTime = (-InitVelocity + sqrt(InitVelocity * InitVelocity + 2.0 * AccBound * PosDiff))/(1.0 * AccBound);
     GoalVelocity = InitVelocity + AccTime * AccBound;
     return AccTime;
   }
   else { // This means that DOF accelerates to limit and then remain on maximum speed
     double AccTime1 = (VelocityBound - InitVelocity)/AccBound;
     double AccTime2 = (PosDiff - AccDist)/VelocityBound;
     double AccTime = AccTime1 + AccTime2;
     GoalVelocity = VelocityBound;
     return AccTime;
   }
 }

 static double NegativePosNVel(  const double & PosDiff,
                                 const double & InitVelocity,  double & GoalVelocity,
                                 const double & VelocityBound, const double & AccBound){
  double AccDist = (VelocityBound * VelocityBound - InitVelocity * InitVelocity)/(2.0 * AccBound);
  if(AccDist>-PosDiff){
    double AccTime = (InitVelocity + sqrt(InitVelocity * InitVelocity - 2.0 * AccBound * PosDiff))/(1.0 * AccBound);
    GoalVelocity = InitVelocity - AccTime * AccBound;
    return AccTime;
  }
  else { // This means that DOF accelerates to limit and then remain on maximum speed
    double AccTime1 = (VelocityBound + InitVelocity)/AccBound;
    double AccTime2 = -(PosDiff + AccDist)/VelocityBound;
    double AccTime = AccTime1 + AccTime2;
    GoalVelocity = -VelocityBound;
    return AccTime;
  }
}

static double AccPhaseTimeInner(const double & PosDiff,
                                const double & InitVelocity,  double & GoalVelocity,
                                const double & VelocityBound, const double & AccBound){
  // The Velocity/Acceleartion bound is assumed to be bidirectional.
  if(PosDiff>0.0){
    if(InitVelocity>0.0){
      return PositivePosNVel(PosDiff, InitVelocity, GoalVelocity, VelocityBound, AccBound);
    }
    else {
      double AccTime1 = -InitVelocity/AccBound;
      double AccDist1 = InitVelocity * InitVelocity/(2.0 * AccBound);
      double AccTime2 = PositivePosNVel(PosDiff + AccDist1, 0.0, GoalVelocity, VelocityBound, AccBound);
      return AccTime1 + AccTime2;
    }
  }
  else{
    if(InitVelocity<0.0){
      return NegativePosNVel(PosDiff, InitVelocity, GoalVelocity, VelocityBound, AccBound);
    } else {
      double AccTime1 = InitVelocity/AccBound;
      double AccDist1 = InitVelocity * InitVelocity/(2.0 * AccBound);
      double AccTime2 = NegativePosNVel(PosDiff - AccDist1, 0.0, GoalVelocity, VelocityBound, AccBound);
      return AccTime1 + AccTime2;
    }
  }
}

static void Velocity2Pos(   const double & PosDiff,
                            const double & InitVelocity,  double & GoalVelocity,
                            const double & VelocityBound, const double & AccBound,
                            const double & DurationTime){
  // Figure out the velocity that this system reaches GoalPos with Duration time.
  if(PosDiff>0.0){
    double AccEst = 2.0 * (PosDiff - InitVelocity * DurationTime)/(DurationTime * DurationTime);  // This value is sure to be less than the bound.
    double VelocityEst = InitVelocity + AccEst * DurationTime;
    if(VelocityEst<VelocityBound) GoalVelocity = VelocityEst;
    else GoalVelocity = VelocityBound;
  }
  else{
    double AccEst = -2.0 * (PosDiff - InitVelocity * DurationTime)/(DurationTime * DurationTime);  // This value is sure to be less than the bound.
    double VelocityEst = InitVelocity - AccEst * DurationTime;
    if(VelocityEst>-VelocityBound) GoalVelocity = VelocityEst;
    else GoalVelocity = -VelocityBound;
  }
}

double AccPhaseTimePathMethod(  const std::vector<double> & CurConfig,      const std::vector<double> & NextConfig,
                                const std::vector<double> & CurVelocity,    std::vector<double> & NextVelocity,
                                const std::vector<double> & VelocityBound,  const std::vector<double> & AccelerationBound,
                                const std::vector<int> SwingLinkChain){
                       // This function solves for the time in acceleration phase.
  /*
    For each DOF, the following compuation procedure is conducted to figure the time.
    1. Determine whether should accelerate in postive or negative direction.
    2. Determine
  */
  std::vector<double> AccPhaseTimeTotal(SwingLinkChain.size());
  std::vector<double> SwingLinkChainVelocity(SwingLinkChain.size());
  for (int i = 0; i < SwingLinkChain.size(); i++) {
    double InitPos = CurConfig[SwingLinkChain[i]];
    double GoalPos = NextConfig[SwingLinkChain[i]];
    double PosDiff = GoalPos - InitPos;
    double InitVelocity = CurVelocity[SwingLinkChain[i]];
    double GoalVelocity;
    double VelBound = VelocityBound[SwingLinkChain[i]];
    double AccBound = AccelerationBound[SwingLinkChain[i]];
    double AccPhaseTime = AccPhaseTimeInner(PosDiff, InitVelocity, GoalVelocity, VelBound, AccBound);
    printf("Link: %d PosDiff: %f,   InitVelocity: %f,   GoalVelocity: %f,   AccTime: %f and Valid: %d\n",
                                  SwingLinkChain[i], PosDiff, InitVelocity, GoalVelocity, AccPhaseTime, GoalVelocity * PosDiff>0.0);
    AccPhaseTimeTotal[i] = AccPhaseTime;
    SwingLinkChainVelocity[i] = GoalVelocity;
  }
  double AccTime = *max_element(AccPhaseTimeTotal.begin(), AccPhaseTimeTotal.end());
  double AccLinkIndex = std::distance(AccPhaseTimeTotal.begin(), std::max_element(AccPhaseTimeTotal.begin(), AccPhaseTimeTotal.end()));
  for (int i = 0; i < SwingLinkChain.size(); i++) {
    double InitPos = CurConfig[SwingLinkChain[i]];
    double GoalPos = NextConfig[SwingLinkChain[i]];
    double PosDiff = GoalPos - InitPos;
    double InitVelocity = CurVelocity[SwingLinkChain[i]];
    double GoalVelocity;
    double VelBound = VelocityBound[SwingLinkChain[i]];
    double AccBound = AccelerationBound[SwingLinkChain[i]];
    Velocity2Pos(PosDiff, InitVelocity, GoalVelocity, VelBound, AccBound, AccTime);
    printf("Updated Link: %d,   PosDiff: %f,  InitVelocity: %f,   GoalVelocity: %f and Valid: %d\n",
            SwingLinkChain[i], PosDiff, InitVelocity, GoalVelocity,  GoalVelocity * PosDiff>0.0);
    NextVelocity[SwingLinkChain[i]] = GoalVelocity;
  }
  return AccTime;
}

double DecPhaseTimePathMethod(  const std::vector<double> & CurConfig,      const std::vector<double> & NextConfig,
                                const std::vector<double> & CurVelocity,    std::vector<double> & NextVelocity,
                                const std::vector<double> & VelocityBound,  const std::vector<double> & AccelerationBound,
                                const std::vector<int> SwingLinkChain){
   // This function solves for the time in deceleration phase.


}
