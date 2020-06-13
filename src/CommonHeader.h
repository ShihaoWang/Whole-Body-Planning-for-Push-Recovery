#ifndef COMMON_HEADER_H
#define COMMON_HEADER_H
#include <iostream>
#include <fstream>
#include <limits.h>
#include <string>
#include <KrisLibrary/geometry/Conversions.h>
#include <KrisLibrary/geometry/MultiVolumeGrid.h>
#include <KrisLibrary/geometry/CollisionMesh.h>
#include <KrisLibrary/geometry/AnyGeometry.h>
#include <KrisLibrary/meshing/VolumeGrid.h>
#include <KrisLibrary/meshing/PointCloud.h>
#include "RobotInfo.h"

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

/* 3. Robot Utilities*/
int FileIndexFinder(bool UpdateFlag);
void FilePathManager(const string & SpecificPath);
std::vector<string>  EdgeFileNamesGene(const string & CurrentCasePath);
Vector3 ImpulseDirectionGene(Robot & SimRobotObj, const std::vector<ContactStatusInfo> & RobotContactInfo, const int & Option);
Vector3 FlatRandomDirection();
void getCentroidalState(const Robot & SimRobot, Vector3 & COMPos, Vector3 & COMVel);
std::vector<Vector3> ActiveContactFinder(const Robot & SimRobot, const std::vector<ContactStatusInfo> & RobotContactInfo);

/* 4. Main Simulation*/
LinearPath InitialSimulation(WorldSimulation & Sim, const SimPara & SimParaObj);
void PushImposer(WorldSimulation & Sim, const double & CurTime, const SimPara & SimParaObj, const bool & FailureFlag);
int SimulationTest(WorldSimulation & Sim, const std::vector<ContactStatusInfo> & _InitContactInfo, ReachabilityMap & RMObject, SelfLinkGeoInfo & SelfLinkGeoObj, const SimPara & _SimParaObj);

/* 5. Convex Polytope Functions*/
FacetInfo FlatConvexHullGeneration(const std::vector<Vector3> & ContactPoints);
std::vector<PIPInfo> PIPGenerator(const std::vector<Vector3> & ContactPoints, const Vector3 & COMPos, const Vector3 & COMVel);
#endif
