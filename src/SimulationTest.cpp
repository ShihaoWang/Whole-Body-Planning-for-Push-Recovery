// This function is used to extensively simulate the result from the four failure metric
#include "RobotInfo.h"
#include "CommonHeader.h"
#include <ode/ode.h>
#include "Control/PathController.h"
#include "Control/JointTrackingController.h"
#include "NonlinearOptimizerInfo.h"

int SimulationTest(WorldSimulation & Sim, const std::vector<ContactStatusInfo> & InitContactInfo, ReachabilityMap & RMObject, SelfLinkGeoInfo & SelfLinkGeoObj, const SimPara & SimParaObj){
  /* Simulation parameters */
  int     DOF             = Sim.world->robots[0]->q.size();
  double  DetectionWait   = SimParaObj.DetectionWait;
  int     PlanStageIndex  = 0;
  double  MPCDuration     = 0.1;                // Duration for MPC executation until next planning
  double  MPCCounter      = MPCDuration;
  bool    MPCFlag         = false;              // True means MPC planning is working.
  bool    TouchDownFlag   = false;

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
  Vector3 COMPos(0.0, 0.0, 0.0), COMVel(0.0, 0.0, 0.0);
  double InitTime = Sim.time;

  ControlReferenceInfo  ControlReferenceObj;                            // Used for control reference generation.
  FailureStateInfo      FailureStateObj;
  Robot                 SimRobot;
  bool                  FailureFlag = false;
  double                FailureMetric;
  double                SimTime;
  std::vector<ContactStatusInfo> curContactInfo =  InitContactInfo;

  while(Sim.time <= TotalDuration){
    SimRobot = *Sim.world->robots[0];
    SimTime = Sim.time;
    PushImposer(Sim,  SimTime - InitTime, SimParaObj, FailureFlag);

    getCentroidalState(SimRobot, COMPos, COMVel);
    std::vector<Vector3> ActContactPos = ActiveContactFinder(SimRobot, curContactInfo);
    std::vector<PIPInfo> PIPTotal = PIPGenerator(ActContactPos, COMPos, COMVel);
    ContactPolytopeWriter(ActContactPos, PIPTotal, SimParaObj);

    FailureMetric = FailureMetricEval(PIPTotal);
    std::printf("Simulation Time: %f, and Failure Metric Value: %f\n", Sim.time, FailureMetric);
    if(FailureMetric < 0.0) FailureFlag = true;
    if((FailureFlag)&&(MPCCounter>=MPCDuration)){
      ControlReferenceObj = ControlReferenceGene(SimRobot, curContactInfo, RMObject, SelfLinkGeoObj, SimParaObj);
    }

    // if((!RHPFlag)&&(FailureFlag)){
    //   qDes = RawOnlineConfigReference(Sim, InitTime, ControlReference, TerrColGeom, SelfLinkGeoObj, DetectionWait, MPCFlag, RobotContactInfo, RMObject);
    //   ControlReference.RunningTime+=TimeStep;
    //   TouchDownFlag = ControlReference.TouchDownTerminalFlag;
    // }
    //
    // if((MPCFlag)&&(MPCCount<MPCDuration)&&(RHPFlag)){
    //   qDes = OnlineConfigReference(Sim, InitTime, ControlReference, TerrColGeom, SelfLinkGeoObj, DetectionWait, MPCFlag, RobotContactInfo, RMObject, MPCCount);
    //   ControlReference.RunningTime+=TimeStep;
    //   TouchDownFlag = ControlReference.TouchDownTerminalFlag;
    //   MPCCount+=TimeStep;
    //   printf("MPC count: %f\n", MPCCount);
    //   if(!MPCFlag) MPCCount = 0.0;
    // }
    // else{
    //   switch (CPPIPIndex){
    //     case -1:
    //     {
    //       if((MPCCount>=MPCDuration)&&(RHPFlag)&&(!ControlReference.TouchDownTerminalFlag)){
    //         MPCCount = 0.0;
    //       }
    //     }
    //     break;
    //     default:{
    //       FailureFlag = true;
    //       if(TouchDownFlag){
    //         if(DetectionWait>=DetectionWait){
    //           InitTime = Sim.time;
    //           ContactStatusOptionRef = -1;
    //           if(!FailureStateObj.FailureInitFlag)  FailureStateObj.FailureStateUpdate(InitTime, SimRobot.q, SimRobot.dq);
    //
    //           double PlanTime;
    //           SelfLinkGeoObj.LinkBBsUpdate(SimRobot);
    //           ControlReference = ControlReferenceGeneration(SimRobot, COMPos, COMVel, RefFailureMetric, RobotContactInfo, RMObject, SelfLinkGeoObj, TimeStep, PlanTime, SpecificPath, PlanStageIndex, DisTol, ContactStatusOptionRef, PreviousContactStatusIndex, Sim.time);
    //           if(ControlReference.ControlReferenceFlag){
    //             PlanTimeRecorder(PlanTime, SpecificPath);
    //             ContactStatusOptionRef = ControlReference.ContactStatusOptionIndex;
    //             DetectionWait = 0.0;
    //             PlanStageIndex++;
    //             MPCFlag = true;
    //             MPCCount = 0.0;
    //             TouchDownFlag = false;
    //           }
    //         }
    //         else DetectionWait+=TimeStep;
    //       }else{  // Then this is the MPC planning
    //         if(RHPFlag&&!ControlReference.TouchDownPhaseFlag){
    //           double PlanTime;
    //           SelfLinkGeoObj.LinkBBsUpdate(SimRobot);
    //           ControlReferenceInfo ControlReferenceMPC = ControlReferenceGeneration(SimRobot, COMPos, COMVel, RefFailureMetric, RobotContactInfo, RMObject, SelfLinkGeoObj, TimeStep, PlanTime, SpecificPath, PlanStageIndex, DisTol, ContactStatusOptionRef, PreviousContactStatusIndex, Sim.time);
    //           if(ControlReferenceMPC.ControlReferenceFlag){
    //             double RunningTime = ControlReference.RunningTime;
    //             ControlReference = ControlReferenceMPC;
    //             ControlReference.RunningTime+=RunningTime;
    //             PlanTimeRecorder(PlanTime, SpecificPath);
    //             ContactStatusOptionRef = ControlReference.ContactStatusOptionIndex;
    //             DetectionWait = 0.0;
    //             PlanStageIndex++;
    //             MPCFlag = true;
    //             MPCCount = 0.0;
    //             TouchDownFlag = false;
    //           }else{
    //             MPCFlag = true;
    //             MPCCount = 0.0;
    //           }
    //         }
    //       }
    //     }
    //   }
    // }
    NewControllerPtr->SetConstant(Config(qDes));
    StateLogger(Sim, FailureStateObj, CtrlStateTraj, PlanStateTraj, FailureStateTraj, qDes, SimParaObj);
    Sim.Advance(SimParaObj.TimeStep);
    Sim.UpdateModel();
  }
  // if(FailureChecker(SimRobot, TerrColGeom, RMObject, DisTol)) PushRecovSuccFlag = 0;
  // else PushRecovSuccFlag = 1;
  //
  // ofstream CtrlStateTrajFile;
  // CtrlStateTrajFile.open (CtrlStateTrajStr_Name);
  // CtrlStateTraj.Save(CtrlStateTrajFile);
  // CtrlStateTrajFile.close();
  //
  // ofstream PlanStateTrajFile;
  // PlanStateTrajFile.open (PlanStateTrajStr_Name);
  // PlanStateTraj.Save(PlanStateTrajFile);
  // PlanStateTrajFile.close();
  //
  // // Then we gonna have to simulate the robot's falling trajectory
  // if(FailureStateObj.FailureInitFlag){
  //   Sim.time = FailureStateObj.FailureTime;
  //
  //   Sim.world->robots[0]->UpdateConfig(FailureStateObj.FailureConfig);
  //   Sim.world->robots[0]->dq = FailureStateObj.FailureVelocity;
  //
  //   Sim.controlSimulators[0].oderobot->SetConfig(FailureStateObj.FailureConfig);
  //   Sim.controlSimulators[0].oderobot->SetVelocities(FailureStateObj.FailureVelocity);
  //
  //   NewControllerPtr->SetConstant(FailureStateObj.FailureConfig);
  //   while(Sim.time <= CtrlStateTraj.EndTime()){
  //     FailureStateTraj.Append(Sim.time,    Sim.world->robots[0]->q);
  //     StateTrajAppender(FailureStateTrajStr_Name, Sim.time, Sim.world->robots[0]->q);
  //     Sim.Advance(TimeStep);
  //     Sim.UpdateModel();
  //   }
  //   if(FailureChecker(*Sim.world->robots[0], TerrColGeom, RMObject, DisTol)) ActualFailureFlag = 1;
  //   else ActualFailureFlag = 0;
  // }
  // else ActualFailureFlag = 0;
  // // Write these three trajectories into files.
  // ofstream FailureStateTrajFile;
  // FailureStateTrajFile.open (FailureStateTrajStr_Name);
  // FailureStateTraj.Save(FailureStateTrajFile);
  // FailureStateTrajFile.close();

  return -1;
}
