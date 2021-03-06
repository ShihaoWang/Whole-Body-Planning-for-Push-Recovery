#ifndef COMMON_HEADER_H
#define COMMON_HEADER_H
#include <iostream>
#include <limits.h>
#include <string>
#include <sys/stat.h>
#include <unistd.h>
#include <string>
#include <ctime>
#include <KrisLibrary/geometry/Conversions.h>
#include <KrisLibrary/geometry/MultiVolumeGrid.h>
#include <KrisLibrary/geometry/CollisionMesh.h>
#include <KrisLibrary/geometry/AnyGeometry.h>
#include <KrisLibrary/meshing/VolumeGrid.h>
#include <KrisLibrary/meshing/PointCloud.h>
#include "RobotInfo.h"
#include "Spline.h"

/* 0. Robot Info Initiaization */
std::vector<LinkInfo> ContactInfoLoader(const string & ContactLinkFile, int & ContactPointNo);
std::vector<ContactStatusInfo> ContactStatusInfoLoader(const string & ContactStatusFile);

/* 1. Environment Geometry */
SignedDistanceFieldInfo SignedDistanceFieldGene(const string & SDFPath, const RobotWorld& WorldObj, const int& GridsNo);
SignedDistanceFieldInfo SignedDistanceFieldLoader(const string & SDFPath, const int GridsNo);
ReachabilityMap ReachabilityMapGenerator(Robot & SimRobot, const std::vector<LinkInfo> & RobotLinkInfo, const std::vector<int> & TorsoLink);

/* 2. Robot State File Operations */
void RobotConfigLoader(Robot &SimRobot, const string &user_path, const string &file_name);
void RobotConfigWriter(const std::vector<double> & Config, const string &user_path, const string &config_file_name);
std::vector<int> TorsoLinkReader(const string & TorsoLinkFilePath);
void StateTrajAppender(const char *stateTrajFile_Name, const double & Time_t, const std::vector<double> & Configuration);
void PushInfoFileAppender(const double & SimTime, const double & Fx_t, const double & Fy_t, const double & Fz_t, const string & SpecificPath);
void StateLogger(WorldSimulation & Sim, FailureStateInfo & FailureStateObj, LinearPath & CtrlStateTraj, LinearPath & PlanStateTraj, LinearPath & FailureStateTraj, std::vector<double> & qDes, const SimPara & SimParaObj);

/* 3. Robot Utilities*/
int FileIndexFinder(bool UpdateFlag);
void FilePathManager(const string & SpecificPath);
std::vector<string>  EdgeFileNamesGene(const string & CurrentCasePath);
Vector3 ImpulseDirectionGene(Robot & SimRobotObj, const std::vector<ContactStatusInfo> & RobotContactInfo, const int & Option);
Vector3 FlatRandomDirection();
void getCentroidalState(const Robot & SimRobot, Vector3 & COMPos, Vector3 & COMVel);
std::vector<Vector3> ActiveContactFinder(const Robot & SimRobot, const std::vector<ContactStatusInfo> & RobotContactInfo);
// std::vector<double> YPRShifter(std::vector<double> OptConfig);
void Vector3Writer(const std::vector<Vector3> & ContactPoints, const std::string & ContactPointFileName);
std::vector<Vector3> BoxVertices(const Box3D & Box3DObj);
std::vector<double> ConfigReferenceGene(const Robot & SimRobotObj,  double & InnerTime,
                                        ReachabilityMap & RMObject, SelfLinkGeoInfo & SelfLinkGeoObj,
                                        ControlReferenceInfo & ControlReference, SimPara & SimParaObj);
bool FailureChecker(Robot & SimRobot, ReachabilityMap & RMObject);
void PlanTimeRecorder(const double & PlanTimeVal, const string & CurrentCasePath);
void PlanningInfoFileAppender(const int & PlanStageIndex, const int & TotalLinkNo, const string & CurrentCasePath, const double & CurTime);
void PlanResWriter(const string & SpecificPath, const int & PushRecovFlag);
void SwingLinkStatePrint(const std::vector<double> & Config, const std::vector<int> & SwingLinkChain);
bool PenetrationTester(const Robot & SimRobotObj, const int & SwingLinkInfoIndex);
void getEndEffectorXYAxes(const Robot & SimRobotInner, const int & SwingLinkInfoIndex, Vector3 & EndEffectorInitxDir, Vector3 & EndEffectorInityDir);
double EstimatedFailureMetric(const Robot & SimRobotInner, const std::vector<ContactStatusInfo> & GoalContactInfo, const Vector3 & COMPos, const Vector3 & COMVel);
void WholeBodyConfigAppender(std::vector<Config> & WholeBodyConfigTraj, const Config & UpdatedConfig);
Vector3 ContactForceFinder(WorldSimulation & Sim);
void ContactForceAppender(const char *CFFile_Name, double Time_t, Vector3 Force);
void KineticEnergyAppender(const char *KEFile_Name, double Time_t, double KE);

/* 4. Main Simulation*/
LinearPath InitialSimulation(WorldSimulation & Sim, const SimPara & SimParaObj);
void PushImposer(WorldSimulation & Sim, const double & CurTime, const SimPara & SimParaObj, const bool & FailureFlag);
int SimulationTest(WorldSimulation & Sim, const std::vector<ContactStatusInfo> & _InitContactInfo, ReachabilityMap & RMObject, SelfLinkGeoInfo & SelfLinkGeoObj, SimPara & _SimParaObj);

/* 5. Convex Polytope Functions*/
FacetInfo FlatConvexHullGeneration(const std::vector<Vector3> & ContactPoints);
std::vector<PIPInfo> PIPGenerator(const std::vector<Vector3> & ContactPoints, const Vector3 & COMPos, const Vector3 & COMVel);
void ContactPolytopeWriter(const std::vector<Vector3> & ActiveContact, const std::vector<PIPInfo> & PIPTotal, const SimPara & SimParaObj);
double FailureMetricEval(const std::vector<PIPInfo> & PIPTotal);
PIPInfo TipOverPIPGenerator(const std::vector<Vector3> & ActiveContacts, const Vector3 & COMPos, const Vector3 & COMVel, bool & ValidFlag);

/* 6. Control Reference Generation*/
ControlReferenceInfo ControlReferenceGene(Robot & SimRobot, const std::vector<ContactStatusInfo> & RobotContactInfo, ReachabilityMap & RMObject, SelfLinkGeoInfo & SelfLinkGeoObj, SimPara & SimParaObj);

/* 7. Whole-Body Estimation */
Config WholeBodyDynamicsIntegrator(Robot & SimRobot, InvertedPendulumInfo & InvertedPendulumObj, const double & TimeDuration);
std::vector<double> CurrentBaseDeltaCal(const Robot & _SimRobot, const InvertedPendulumInfo & _InvertedPendulumObj, const double & TimeDuration);

/* 8. Transient Path Generation*/
CubicSplineInfo TransientPathGene(const Robot & SimRobot, SelfLinkGeoInfo & SelfLinkGeoObj, SimPara & SimParaObj);
ControlReferenceInfo TrajectoryPlanning(Robot & SimRobotInner, const InvertedPendulumInfo & InvertedPendulumObj, ReachabilityMap & RMObject,SelfLinkGeoInfo & SelfLinkGeoObj,
                                        EndEffectorPathInfo & EndEffectorPathObj, SimPara & SimParaObj);
std::vector<double> TrajConfigOptimazation(const Robot & SimRobot, ReachabilityMap & RMObject, SelfLinkGeoInfo & _SelfLinkGeoObj, SimPara & SimParaObj, const double & EndEffectorProjx, const double & EndEffectorProjy, double sVal);
std::vector<double> LastStageConfigOptimazation(const Robot & SimRobot, ReachabilityMap & RMObject, SelfLinkGeoInfo & _SelfLinkGeoObj, SimPara & SimParaObj, const int & StageIndex);
ControlReferenceInfo TaskTrajectoryPlanning(Robot & SimRobotInner, const InvertedPendulumInfo & InvertedPendulumInner, ReachabilityMap & RMObject,SelfLinkGeoInfo & SelfLinkGeoObj,
                                            EndEffectorPathInfo & EndEffectorPathObj, SimPara & SimParaObj);
bool TaskTrajectoryPlanningInner( const double & _sVal, double & _sNew,
                                  Robot & _SimRobotInner,
                                  const Config & _PreVelocity,                      const std::vector<double> & _CurrentBaseDelta,
                                  const Config & _CurrentConfig,                    const Config & _CurrentVelocity,
                                  std::vector<double> & _NextConfig,                std::vector<double> & _NextVelocity,
                                  const EndEffectorPathInfo & EndEffectorPathInner, const std::vector<int> & _SwingLinkChain,
                                  const Vector3 & _EndEffectorInitxDir,             const Vector3 & _EndEffectorInityDir,
                                  SimPara & _SimParaObj,
                                  const double & _StageTime,                        const double & _DampingRatio);
#endif
