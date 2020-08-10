// This function is used to extensively simulate the result from the four failure metric
#include "RobotInfo.h"
#include "CommonHeader.h"
#include <ode/ode.h>
#include "Control/PathController.h"
#include "Control/JointTrackingController.h"
#include "NonlinearOptimizerInfo.h"

int SimulationTest(WorldSimulation & Sim, const std::vector<ContactStatusInfo> & InitContactInfo, ReachabilityMap & RMObject, SelfLinkGeoInfo & SelfLinkGeoObj, SimPara & SimParaObj){
  /* Simulation parameters */
  int     DOF             = Sim.world->robots[0]->q.size();
  double  DetectionWait   = SimParaObj.DetectionWait;
  double  DetectionCount  = DetectionWait;
  int     PlanStageIndex  = 0;
  double  RHPDuration     = 0.25;                 // Duration for MPC executation until next planning
  double  RHPCounter      = RHPDuration;
  bool    RHPFlag         = false;                // True means MPC planning is working.
  bool    TouchDownFlag   = false;
  bool    CtrlFlag        = false;
  double  CtrlStartTime   = 0.0;

  /* Override the default controller with a PolynomialPathController */
  auto NewControllerPtr = std::make_shared<PolynomialPathController>(*Sim.world->robots[0]);
  Sim.SetController(0, NewControllerPtr);
  NewControllerPtr->SetConstant(Sim.world->robots[0]->q);

  // Initial Simulation
  LinearPath InitTraj, FailureStateTraj, CtrlStateTraj, PlanStateTraj;
  InitTraj  = InitialSimulation(Sim, SimParaObj);
  FailureStateTraj = InitTraj;
  CtrlStateTraj = InitTraj;
  PlanStateTraj = InitTraj;
  double TotalDuration = SimParaObj.TotalDuration + Sim.time;
  std::printf("Initial Simulation Done!\n");
  std::vector<double> qDes = PlanStateTraj.milestones[PlanStateTraj.milestones.size()-1];               // This is commanded robot configuration to the controller.
  double InitTime = Sim.time;

  ControlReferenceInfo  ControlReferenceObj;                            // Used for control reference generation.
  FailureStateInfo      FailureStateObj;
  Robot                 SimRobot;
  bool                  FailureFlag = false;
  bool                  OverallFailureFlag = false;
  double                FailureMetric;
  double                SimTime;
  std::vector<ContactStatusInfo> curContactInfo =  InitContactInfo;
  while(Sim.time <= TotalDuration){
    SimRobot = *Sim.world->robots[0];
    SimTime = Sim.time;
    if(!FailureFlag) PushImposer(Sim,  SimTime - InitTime, SimParaObj, FailureFlag);

    Vector3 COMPos, COMVel;
    getCentroidalState(SimRobot, COMPos, COMVel);
    std::vector<Vector3> ActContactPos = ActiveContactFinder(SimRobot, curContactInfo);
    std::vector<PIPInfo> PIPTotal = PIPGenerator(ActContactPos, COMPos, COMVel);
    ContactPolytopeWriter(ActContactPos, PIPTotal, SimParaObj);

    SelfLinkGeoObj.LinkBBsUpdate(SimRobot);
    if(!OverallFailureFlag){
      if(!CtrlFlag){
        FailureMetric = FailureMetricEval(PIPTotal);
        if(DetectionCount>=DetectionWait){
          std::printf("Simulation Time: %f, and Failure Metric Value: %f\n", Sim.time, FailureMetric);
          if(FailureMetric < 0.0){
            if(!FailureStateObj.FailureStateFlag)  FailureStateObj.FailureStateUpdate(SimTime, SimRobot.q, SimRobot.dq);
            FailureFlag = true;
          }
          if(FailureFlag){
            SimParaObj.setPlanStageIndex(PlanStageIndex);
            SimParaObj.setSimTime(SimTime);
            ControlReferenceObj = ControlReferenceGene(SimRobot, curContactInfo, RMObject, SelfLinkGeoObj, SimParaObj);
            if(ControlReferenceObj.getReadyFlag()){
              CtrlFlag = true;
              CtrlStartTime = SimTime;
              PlanStageIndex++;
            }
          }
        } else{
          DetectionCount+=SimParaObj.TimeStep;
          std::printf("Simulation Time: %f, and Failure Metric Value: %f\n", Sim.time, FailureMetric);
        }
       } else {
        double InnerTime = SimTime - CtrlStartTime;
        qDes =  ConfigReferenceGene(SimRobot, InnerTime, RMObject, SelfLinkGeoObj, ControlReferenceObj, SimParaObj);
        if (ControlReferenceObj.getTouchDownFlag()){
          CtrlFlag = false;
          DetectionCount = 0.0;
          curContactInfo = ControlReferenceObj.GoalContactStatus;
          FailureFlag = false;
        }
      }
    }

    OverallFailureFlag = FailureChecker(SimRobot, RMObject);

    NewControllerPtr->SetConstant(Config(qDes));
    StateLogger(Sim, FailureStateObj, CtrlStateTraj, PlanStateTraj, FailureStateTraj, qDes, SimParaObj);
    Sim.Advance(SimParaObj.TimeStep);
    Sim.UpdateModel();
  }
  int PushRecovSuccFlag = 0;
  if(FailureChecker(SimRobot, RMObject)) PushRecovSuccFlag = 0;
  else PushRecovSuccFlag = 1;

  ofstream PlanStateTrajFile;
  PlanStateTrajFile.open (SimParaObj.PlanStateTrajStr.c_str());
  PlanStateTraj.Save(PlanStateTrajFile);
  PlanStateTrajFile.close();

  ofstream CtrlStateTrajFile;
  CtrlStateTrajFile.open (SimParaObj.CtrlStateTrajStr.c_str());
  CtrlStateTraj.Save(CtrlStateTrajFile);
  CtrlStateTrajFile.close();

  // Then we gonna have to simulate the robot's falling trajectory
  if(FailureStateObj.FailureStateFlag){
    Sim.time = FailureStateObj.FailureTime;
    Sim.world->robots[0]->UpdateConfig(FailureStateObj.FailureConfig);
    Sim.world->robots[0]->dq = FailureStateObj.FailureVelocity;
    Sim.controlSimulators[0].oderobot->SetConfig(FailureStateObj.FailureConfig);
    Sim.controlSimulators[0].oderobot->SetVelocities(FailureStateObj.FailureVelocity);
    NewControllerPtr->SetConstant(FailureStateObj.FailureConfig);
    while(Sim.time <= CtrlStateTraj.EndTime()){
      FailureStateTraj.Append(Sim.time,    Sim.world->robots[0]->q);
      StateTrajAppender(SimParaObj.FailureStateTrajStr.c_str(), Sim.time, Sim.world->robots[0]->q);
      Sim.Advance(SimParaObj.TimeStep);
      Sim.UpdateModel();
    }
  }
  int FailureSuccFlag = 0;
  if(FailureChecker(*Sim.world->robots[0], RMObject)) FailureSuccFlag = 1;

  // Write these three trajectories into files.
  ofstream FailureStateTrajFile;
  FailureStateTrajFile.open (SimParaObj.FailureStateTrajStr.c_str());
  FailureStateTraj.Save(FailureStateTrajFile);
  FailureStateTrajFile.close();

  if((FailureSuccFlag) && (PushRecovSuccFlag)) return 1;
  else return 0;
}
