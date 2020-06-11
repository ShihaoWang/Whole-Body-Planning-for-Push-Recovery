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

/* 3. Robot Utilities*/
int FileIndexFinder(bool UpdateFlag);
void FilePathManager(const string & SpecificPath);
std::vector<string>  EdgeFileNamesGene(const string & CurrentCasePath);

/* 4. Main Simulation*/
LinearPath InitialSimulation(WorldSimulation & Sim, const SimPara & SimParaObj);
int SimulationTest(WorldSimulation & Sim, const std::vector<ContactStatusInfo> & _InitContactInfo, ReachabilityMap & RMObject, SelfLinkGeoInfo & SelfLinkGeoObj, const SimPara & _SimParaObj);

#endif
