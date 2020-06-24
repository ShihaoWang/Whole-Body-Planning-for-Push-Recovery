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

void Vector3Writer(const std::vector<Vector3> & ContactPoints, const std::string & ContactPointFileName){
  if(!ContactPoints.size()) return;
  int NumberOfContactPoints = ContactPoints.size();
  std::vector<double> FlatContactPoints(3 * NumberOfContactPoints);
  int FlatContactPointIndex = 0;
  for (int i = 0; i < NumberOfContactPoints; i++){
    FlatContactPoints[FlatContactPointIndex] = ContactPoints[i].x;
    FlatContactPointIndex++;
    FlatContactPoints[FlatContactPointIndex] = ContactPoints[i].y;
    FlatContactPointIndex++;
    FlatContactPoints[FlatContactPointIndex] = ContactPoints[i].z;
    FlatContactPointIndex++;
  }
  FILE * FlatContactPointsFile = NULL;
  string ContactPointFile = ContactPointFileName + ".bin";
  const char *ContactPointFile_Name = ContactPointFile.c_str();
  FlatContactPointsFile = fopen(ContactPointFile_Name, "wb");
  fwrite(&FlatContactPoints[0], sizeof(double), FlatContactPoints.size(), FlatContactPointsFile);
  fclose(FlatContactPointsFile);
  return;
}

static vector<Vector3> Interpolation(const Vector3 & a, const Vector3 & b, const int & No){
  vector<Vector3> Pts(No);
  double ratio = (1.0)/(1.0 * No + 1.0);
  for (int i = 0; i < No; i++) {
    Pts[i] = a + (b-a) * ratio * (1.0 * i + 1.0);
  }return Pts;
}

std::vector<Vector3> BoxVertices(const Box3D & Box3DObj){
  // Here 4 more points will be sampled along each edge
  std::vector<Vector3> Vertices;
  Vector3 GlobalPoint;
  Vector3 LocalPoint(0.0, 0.0, 0.0);
  Box3DObj.fromLocal(LocalPoint, GlobalPoint);
  Vector3 xBox3DOffset = Box3DObj.dims.x * Box3DObj.xbasis;
  Vector3 yBox3DOffset = Box3DObj.dims.y * Box3DObj.ybasis;
  Vector3 zBox3DOffset = Box3DObj.dims.z * Box3DObj.zbasis;

  Vector3 A = GlobalPoint;
  Vector3 B = GlobalPoint + xBox3DOffset;
  Vector3 C = GlobalPoint + xBox3DOffset + zBox3DOffset;
  Vector3 D = GlobalPoint + zBox3DOffset;

  Vector3 E = GlobalPoint + yBox3DOffset;
  Vector3 F = GlobalPoint + xBox3DOffset + yBox3DOffset;
  Vector3 G = GlobalPoint + xBox3DOffset + zBox3DOffset + yBox3DOffset;
  Vector3 H = GlobalPoint + zBox3DOffset + yBox3DOffset;

  const int InterPtNo = 0;
  vector<Vector3> ABEdge = Interpolation(A, B, InterPtNo);
  vector<Vector3> BCEdge = Interpolation(B, C, InterPtNo);
  vector<Vector3> CDEdge = Interpolation(C, D, InterPtNo);
  vector<Vector3> DAEdge = Interpolation(D, A, InterPtNo);

  vector<Vector3> EFEdge = Interpolation(E, F, InterPtNo);
  vector<Vector3> FGEdge = Interpolation(F, G, InterPtNo);
  vector<Vector3> GHEdge = Interpolation(G, H, InterPtNo);
  vector<Vector3> HEEdge = Interpolation(H, E, InterPtNo);

  Vertices.push_back(A);
  Vertices.insert(Vertices.end(), ABEdge.begin(), ABEdge.end());

  Vertices.push_back(B);
  Vertices.insert(Vertices.end(), BCEdge.begin(), BCEdge.end());

  Vertices.push_back(C);
  Vertices.insert(Vertices.end(), CDEdge.begin(), CDEdge.end());

  Vertices.push_back(D);
  Vertices.insert(Vertices.end(), DAEdge.begin(), DAEdge.end());

  Vertices.push_back(E);
  Vertices.insert(Vertices.end(), EFEdge.begin(), EFEdge.end());

  Vertices.push_back(F);
  Vertices.insert(Vertices.end(), FGEdge.begin(), FGEdge.end());

  Vertices.push_back(G);
  Vertices.insert(Vertices.end(), GHEdge.begin(), GHEdge.end());

  Vertices.push_back(H);
  Vertices.insert(Vertices.end(), HEEdge.begin(), HEEdge.end());
  return Vertices;
}

bool FailureChecker(Robot & SimRobot, ReachabilityMap & RMObject){
  double DistTol = 0.01;
  std::vector<double> LinkTerrDistVec;
  for (int i = 5; i < SimRobot.q.size(); i++){
    if(!RMObject.EndEffectorIndices.count(i)){
      if(SimRobot.geometry[i]->Distance(NonlinearOptimizerInfo::TerrColGeom)<DistTol)
        return true;
    }
  }
  return false;
}

void PlanResWriter(const string & SpecificPath, const int & PushRecovFlag){
  const string PlanResStr = SpecificPath + "PlanRes.txt";
  const char *PlanResStr_Name = PlanResStr.c_str();
  std::ofstream PlanResWriter;
  PlanResWriter.open(PlanResStr_Name);
  PlanResWriter<<std::to_string(PushRecovFlag)<<"\n";
  PlanResWriter.close();
}
