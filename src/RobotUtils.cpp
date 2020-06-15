// This function is used to calculate certain robot utility functions
#include "RobotInfo.h"
#include "CommonHeader.h"
#include "NonlinearOptimizerInfo.h"
#include <KrisLibrary/robotics/Inertia.h>
#include <random>
#include <limits>
#include <sys/stat.h>

int FileIndexFinder(bool UpdateFlag){
  string FileIndexName = "AtNo.txt";         // This file should be located in the "build" folder.
  ifstream FileIndexReader(FileIndexName);
  int FileIndex;
  string str_line;
  if (FileIndexReader.is_open()){
    while (getline (FileIndexReader,str_line))
    FileIndex = stoi(str_line);
    FileIndexReader.close();
  }
  else std:cerr<< "Unable to open FileIndex file";
  if(UpdateFlag){
    const char *FileIndexWriter_Name = FileIndexName.c_str();
    std::ofstream FileIndexWriter;
    FileIndexWriter.open(FileIndexWriter_Name);
    FileIndexWriter<<std::to_string(FileIndex + 1)<<"\n";
    FileIndexWriter.close();
  }
  return FileIndex;
}

static bool IsPathExist(const std::string &s){
  struct stat buffer;
  return (stat (s.c_str(), &buffer) == 0);
}

std::vector<double> YPRShifter(std::vector<double> OptConfig){
  for (int i = 3; i < 6; i++){
    if(OptConfig[i]>M_PI)
      OptConfig[i]-=2.0 * M_PI;
    if(OptConfig[i]<-M_PI)
      OptConfig[i]+=2.0 * M_PI;
  }
  return OptConfig;
}

void FilePathManager(const string & SpecificPath){
  if(IsPathExist(SpecificPath))
    printf("%s exist!\n", SpecificPath.c_str());
  else {
    string str = "mkdir " + SpecificPath;
    const char *command = str.c_str();
    system(command);
  }
  // Let them be internal objects
  string str = "cd " + SpecificPath + " && ";

  str+="rm -f *Traj.txt && ";
  str+="rm -f *.path && ";
  str+="rm -f *InfoFile.txt && ";
  str+="rm -f PlanTime.txt && ";
  str+="rm -f *.bin && ";
  str+="rm -f *OptConfig*.config && ";
  str+="rm -f PlanRes.txt";
  const char *command = str.c_str();
  system(command);
}

Vector3 ImpulseDirectionGene(Robot & SimRobotObj, const std::vector<ContactStatusInfo> & RobotContactInfo, const int & Option){
  std::vector<LinkInfo> RobotLinkInfo = NonlinearOptimizerInfo::RobotLinkInfo;
  Vector3 ImpulseDirection(0.0, 0.0, 0.0);
  if(Option == 1){
    Vector3 A, B; // Towards the direction where the foot is on air.
    SimRobotObj.GetWorldPosition(RobotLinkInfo[0].AvgLocalContact, RobotLinkInfo[0].LinkIndex, A);
    SimRobotObj.GetWorldPosition(RobotLinkInfo[1].AvgLocalContact, RobotLinkInfo[1].LinkIndex, B);
    if(RobotContactInfo[0].LocalContactStatus[0])
          ImpulseDirection = B - A;
    else  ImpulseDirection = A - B;
  } else {
    if(Option == 2){
      std::vector<Vector3> SPVertices;
      for (int i = 0; i < RobotLinkInfo.size(); i++){
        int LinkiPNo = RobotLinkInfo[i].LocalContacts.size();
        for (int j = 0; j < LinkiPNo; j++){
          if(RobotContactInfo[i].LocalContactStatus[j]){
            Vector3 LinkiPjPos;
            SimRobotObj.GetWorldPosition(RobotLinkInfo[i].LocalContacts[j], RobotLinkInfo[i].LinkIndex, LinkiPjPos);
            LinkiPjPos.z = 0.0;
            SPVertices.push_back(LinkiPjPos);
          }
        }
      }
      Vector3 COM_Pos = SimRobotObj.GetCOM();
      FacetInfo SPObj = FlatConvexHullGeneration(SPVertices);    // This is the support polygon
      COM_Pos.z = 0.0;
      std::vector<double> DistVec = SPObj.ProjPoint2EdgeDistVec(COM_Pos);
      std::vector<int> DistVecIndices(DistVec.size());
      for (int i = 0; i < DistVec.size(); i++)
      {
        DistVec[i] = DistVec[i] * DistVec[i];
        DistVecIndices[i] = i;
      }
      int MinIndex = std::distance(DistVec.begin(), std::min_element(DistVec.begin(), DistVec.end()));
      ImpulseDirection = -SPObj.EdgeNorms[MinIndex];
    } else ImpulseDirection = FlatRandomDirection();
  }
  ImpulseDirection.setNormalized(ImpulseDirection);
  return ImpulseDirection;
}

static double RandomBoundedValue(const double &bound){
  std::uniform_real_distribution<double> unif(-1.0 * bound, 1.0 * bound);
  std::random_device rand_dev;          // Use random_device to get a random seed.
  std::mt19937 rand_engine(rand_dev()); // mt19937 is a good pseudo-random number generator.
  double boundval = unif(rand_engine);
  return boundval;
}

Vector3 FlatRandomDirection(){
  double xDir = RandomBoundedValue(1.0);
  double yDir = RandomBoundedValue(1.0);
  Vector3 Dir(xDir, yDir, 0.0);
  Dir.setNormalized(Dir);
  return Dir;
}

void getCentroidalState(const Robot & SimRobot, Vector3 & COMPos, Vector3 & COMVel){
  // This function is used to get the centroidal position and velocity
  Vector3 COM = SimRobot.GetCOM();
  Matrix pCOMpq;
  SimRobot.GetCOMJacobian(pCOMpq);
  double COMVel_x =0.0, COMVel_y =0.0, COMVel_z =0.0;
  for (int i = 0; i < pCOMpq.n; i++){
    COMVel_x = COMVel_x + pCOMpq(0,i) * SimRobot.dq[i];
    COMVel_y = COMVel_y + pCOMpq(1,i) * SimRobot.dq[i];
    COMVel_z = COMVel_z + pCOMpq(2,i) * SimRobot.dq[i];
  }
  Vector3 COMVel_i(COMVel_x, COMVel_y, COMVel_z);
  COMPos = COM;
  COMVel = COMVel_i;
}

std::vector<Vector3> ActiveContactFinder(const Robot & SimRobot, const std::vector<ContactStatusInfo> & RobotContactInfo){
  std::vector<Vector3> ActContacts;
  for (int i = 0; i < NonlinearOptimizerInfo::RobotLinkInfo.size(); i++){
    for (int j = 0; j < NonlinearOptimizerInfo::RobotLinkInfo[i].LocalContacts.size(); j++){
      if(RobotContactInfo[i].LocalContactStatus[j]){
        Vector3 LinkiPjPos;
        SimRobot.GetWorldPosition(NonlinearOptimizerInfo::RobotLinkInfo[i].LocalContacts[j], NonlinearOptimizerInfo::RobotLinkInfo[i].LinkIndex, LinkiPjPos);
        ActContacts.push_back(LinkiPjPos);
      }
    }
  }
  return ActContacts;
}
