#include "RobotInfo.h"
#include "CommonHeader.h"
#include <ode/ode.h>
#include "Control/PathController.h"
#include "Control/JointTrackingController.h"
#include "NonlinearOptimizerInfo.h"

LinearPath InitialSimulation(WorldSimulation & Sim, const SimPara & SimParaObj){
  const char *FailureStateTrajStr_Name  = SimParaObj.FailureStateTrajStr.c_str();
  const char *CtrlStateTrajStr_Name     = SimParaObj.CtrlStateTrajStr.c_str();
  const char *PlanStateTrajStr_Name     = SimParaObj.PlanStateTrajStr.c_str();
  LinearPath InitTraj;
  while(Sim.time <= SimParaObj.InitDuration){
    InitTraj.Append(Sim.time, Sim.world->robots[0]->q);
    std::printf("Initial Simulation Time: %f\n", Sim.time);
    StateTrajAppender(FailureStateTrajStr_Name, Sim.time, Sim.world->robots[0]->q);
    StateTrajAppender(CtrlStateTrajStr_Name,    Sim.time, Sim.world->robots[0]->q);
    StateTrajAppender(PlanStateTrajStr_Name,    Sim.time, Sim.world->robots[0]->q);
    Sim.Advance(SimParaObj.TimeStep);
    Sim.UpdateModel();
  }
  return InitTraj;
}

void PushImposer(WorldSimulation & Sim, const double & CurTime, const SimPara & SimParaObj, const bool & FailureFlag){
  if((CurTime<SimParaObj.PushDuration)&&(!FailureFlag)){
    int LinkIndex = 19;
    double ImpulseScale = 1.0 * CurTime/SimParaObj.PushDuration;
    Vector3 ImpulseForce = ImpulseScale * SimParaObj.ImpulseForceMax;
    dBodyAddForceAtRelPos(Sim.odesim.robot(0)->body(LinkIndex), ImpulseForce.x, ImpulseForce.y, ImpulseForce.z, 0.0, 0.0, 0.0);
    PushInfoFileAppender(Sim.time, ImpulseForce.x, ImpulseForce.y, ImpulseForce.z, SimParaObj.CurrentCasePath);
  }
}

std::vector<double> ConfigReferenceGene(const Robot & SimRobotObj,  double & InnerTime,
                                        ReachabilityMap & RMObject, SelfLinkGeoInfo & SelfLinkGeoObj,
                                        ControlReferenceInfo & ControlReference, SimPara & SimParaObj){
  // This function generates robot's reference configuration at each time.
  double TouchDownTol  = 0.005;                      //  5 mm as a Touch Down Terminal Tolerance.
  std::vector<double> qDes;
  int SwingLinkInfoIndex = ControlReference.getSwingLinkInfoIndex();
  std::vector<double> EndEffectorSDVec;
  for (Vector3 & LocalContact : NonlinearOptimizerInfo::RobotLinkInfo[SwingLinkInfoIndex].LocalContacts) {
    Vector3 SwingLinkContactPos;
    SimRobotObj.GetWorldPosition(LocalContact, NonlinearOptimizerInfo::RobotLinkInfo[SwingLinkInfoIndex].LinkIndex, SwingLinkContactPos);
    double CurrentDist = NonlinearOptimizerInfo::SDFInfo.SignedDistance(SwingLinkContactPos);
    EndEffectorSDVec.push_back(CurrentDist);
  }
  double SwingContactDist = *min_element(EndEffectorSDVec.begin(), EndEffectorSDVec.end());
  Config qDesConfig;
  ControlReference.PlannedConfigTraj.Eval(InnerTime, qDesConfig);     // A problem with this interpolation is a bad visualization due to Euler Angle singularity.
  // double Yaw, Pitch, Roll;
  // ControlReference.EulerAngleInterpolator(InnerTime, Yaw, Pitch, Roll);
  // qDesConfig[3] = Yaw;
  // qDesConfig[4] = Pitch;
  // qDesConfig[5] = Roll;

  double TimeRatio = 0.5;

  if(!ControlReference.getTouchDownFlag()){
    for (int i = 0; i < SimRobotObj.q.size(); i++)
        qDes.push_back(qDesConfig[i]);
  } else qDes = ControlReference.getTouchDownConfig();

  if(!ControlReference.getTouchDownFlag()){
    if(ControlReference.ControlReferenceType==1){
      // Contact Addition Case
      if(SwingContactDist<TouchDownTol){
        ControlReference.setTouchDownFlag(true);
        ControlReference.setTouchDownConfig(qDes);
      }
    }
    else{
      // Contact Modification Case
      if(InnerTime>TimeRatio * ControlReference.TimeTraj.back()){
        if(SwingContactDist<TouchDownTol){
          ControlReference.setTouchDownFlag(true);
          ControlReference.setTouchDownConfig(qDes);
        }
      }
    }
  }
  return qDes;
}
