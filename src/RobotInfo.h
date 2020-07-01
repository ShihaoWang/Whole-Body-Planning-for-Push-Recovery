#ifndef ROBOTINFO_H
#define ROBOTINFO_H
#include <KrisLibrary/robotics/RobotDynamics3D.h>
#include <KrisLibrary/robotics/NewtonEuler.h>
#include <Interface/SimulationGUI.h>
#include <unsupported/Eigen/CXX11/Tensor>
#include <KrisLibrary/meshing/TriMeshTopology.h>
#include <KrisLibrary/math3d/geometry3d.h>
#include <KrisLibrary/geometry/CollisionMesh.h>
#include <KrisLibrary/geometry/PQP/src/PQP.h>
#include "Modeling/Paths.h"
#include "Splines.h"
#include <fstream>

struct LinkInfo {
  LinkInfo(){LinkIndex = -1;}
  LinkInfo(const int &link_index){ LinkIndex = link_index; }
  void AddLocalConact(const Vector3& local_contact){ LocalContacts.push_back(local_contact); }
  void AvgContactUpdate(){
    switch (LocalContacts.size()){
      case 0: throw std::invalid_argument( "LocalContacts should have been initialized!" );
      break;
      default:{
        Vector3 SumLocalContacts(0.0, 0.0, 0.0);
        for (int i = 0; i < LocalContacts.size(); i++){
          SumLocalContacts.x = SumLocalContacts.x + LocalContacts[i].x;
          SumLocalContacts.y = SumLocalContacts.y + LocalContacts[i].y;
          SumLocalContacts.z = SumLocalContacts.z + LocalContacts[i].z;
        }
        AvgLocalContact.x = SumLocalContacts.x/LocalContacts.size();
        AvgLocalContact.y = SumLocalContacts.y/LocalContacts.size();
        AvgLocalContact.z = SumLocalContacts.z/LocalContacts.size();
      }
      break;
    }
  }
  int LinkIndex;
  std::vector<Vector3> LocalContacts;
  std::vector<Vector3> ContactPositions;
  std::vector<Vector3> ContactVelocities;
  std::vector<double> ContactDists;
  Vector3 AvgLocalContact;
};

struct ContactStatusInfo{
  ContactStatusInfo(){ LinkIndex = -1;}
  ContactStatusInfo(const int & link_index){ LinkIndex = link_index; }
  void AddLocalConactStatus(const int & _contactstatus){ LocalContactStatus.push_back(_contactstatus); }
  void StatusSwitch(const int & Val){
    for(int i = 0; i<LocalContactStatus.size(); i++){
      LocalContactStatus[i] = Val;
    }
  }
  int LinkIndex;
  std::vector<int> LocalContactStatus;
};

struct TerrainInfo{
  TerrainInfo(){num_tris = -1;}
  int num_tris;
  std::vector<Tri> Tris;
  std::vector<Vector3> TriNormals;
  std::vector<int> Indices;
};

struct SignedDistanceFieldInfo{
  SignedDistanceFieldInfo(){
    Envi_x_min = 0;           Envi_x_max = 0;
    Envi_y_min = 0;           Envi_y_max = 0;
    Envi_z_min = 0;           Envi_z_max = 0;
    Envi_x_unit = 0;          Envi_y_unit = 0;          Envi_z_unit = 0;
    Envi_x_length = 0;        Envi_y_length = 0;        Envi_z_length = 0;
    GridNo = 0;
  }
  SignedDistanceFieldInfo(const Eigen::Tensor<double, 3>& _SDFTensor, const std::vector<double> &_SDFSpecs){
    SDFTensor = _SDFTensor;
    Envi_x_min = _SDFSpecs[0];          Envi_x_max = _SDFSpecs[1];
    Envi_y_min = _SDFSpecs[2];          Envi_y_max = _SDFSpecs[3];
    Envi_z_min = _SDFSpecs[4];          Envi_z_max = _SDFSpecs[5];
    Envi_x_unit = _SDFSpecs[6];         Envi_y_unit = _SDFSpecs[7];         Envi_z_unit = _SDFSpecs[8];
    Envi_x_length = _SDFSpecs[9];       Envi_y_length = _SDFSpecs[10];      Envi_z_length = _SDFSpecs[11];
    GridNo = (int)_SDFSpecs[12];
  }
  double SignedDistance(const Vector3 &Point) const{
    double x_FloatIndex = (Point.x - Envi_x_min)/Envi_x_unit * 1.0;
    double y_FloatIndex = (Point.y - Envi_y_min)/Envi_y_unit * 1.0;
    double z_FloatIndex = (Point.z - Envi_z_min)/Envi_z_unit * 1.0;

    int x_leftindex = std::floor(x_FloatIndex);
    int y_leftindex = std::floor(y_FloatIndex);
    int z_leftindex = std::floor(z_FloatIndex);

    if(x_leftindex<0){
      x_leftindex = 0;
    } else {
      if(x_leftindex>GridNo-2) {
        x_leftindex = GridNo-2;
      }
    }

    if(y_leftindex<0) {
      y_leftindex = 0;
    } else {
      if(y_leftindex>GridNo-2) {
        y_leftindex = GridNo-2;
      }
    }

    if(z_leftindex<0) {
      z_leftindex = 0;
    } else {
      if(z_leftindex>GridNo-2) {
        z_leftindex = GridNo-2;
      }
    }

    int x_rightindex = x_leftindex + 1;
    int y_rightindex = y_leftindex + 1;
    int z_rightindex = z_leftindex + 1;

    double valA = SDFTensor(x_leftindex, y_leftindex, z_leftindex);
    double valB = SDFTensor(x_rightindex, y_leftindex, z_leftindex);
    double valC = SDFTensor(x_rightindex, y_rightindex, z_leftindex);
    double valD = SDFTensor(x_leftindex, y_rightindex, z_leftindex);

    double valE = SDFTensor(x_leftindex, y_leftindex, z_rightindex);
    double valF = SDFTensor(x_rightindex, y_leftindex, z_rightindex);
    double valG = SDFTensor(x_rightindex, y_rightindex, z_rightindex);
    double valH = SDFTensor(x_leftindex, y_rightindex, z_rightindex);

    // Since this is a tri-linear interpolation, there are three ways to do the interpolation

    // Type1:
    // Along x-direction
    double valMAB = (x_FloatIndex - x_leftindex*1.0) * (valB - valA) + valA;
    double valMDC = (x_FloatIndex - x_leftindex*1.0) * (valC - valD) + valD;
    double valMEF = (x_FloatIndex - x_leftindex*1.0) * (valF - valE) + valE;
    double valMHG = (x_FloatIndex - x_leftindex*1.0) * (valG - valH) + valH;

    // Along y-drection
    double valMABDC = (y_FloatIndex - y_leftindex*1.0) * (valMDC - valMAB) + valMAB;
    double valMEFGH = (y_FloatIndex - y_leftindex*1.0) * (valMHG - valMEF) + valMEF;

    // Along z-direction
    double valMABCDEFHG = (z_FloatIndex - z_leftindex*1.0) * (valMEFGH - valMABDC) + valMABDC;

    // // Type2:
    // // Along y-drection
    // double valMAD = (y_FloatIndex - y_leftindex*1.0) * (valD - valA) + valA;
    // double valMBC = (y_FloatIndex - y_leftindex*1.0) * (valC - valB) + valB;
    // double valMEH = (y_FloatIndex - y_leftindex*1.0) * (valH - valE) + valE;
    // double valMFG = (y_FloatIndex - y_leftindex*1.0) * (valG - valF) + valF;
    //
    // //Along x-direction
    // double valMADBC = (x_FloatIndex - x_leftindex * 1.0) * (valMBC - valMAD) + valMAD;
    // double valMEHFG = (x_FloatIndex - x_leftindex * 1.0) * (valMFG - valMEH) + valMEH;
    //
    // //Along z-direction
    // double valMADBCEHFG = (z_FloatIndex - z_leftindex*1.0) * (valMEHFG - valMADBC) + valMADBC;

    return valMABCDEFHG;
  }
  Vector3 SignedDistanceNormal(const Vector3 &Point) const {
    // This function is used to calculate the (1 x 3) Jacobian matrix given the current Position
    // This function is used to compute the distance from a 3D point to the environment terrain
    // The first job is to figure out the nearest neighbours of the Points

    double x_FloatIndex = (Point.x - Envi_x_min)/Envi_x_unit * 1.0;
    double y_FloatIndex = (Point.y - Envi_y_min)/Envi_y_unit * 1.0;
    double z_FloatIndex = (Point.z - Envi_z_min)/Envi_z_unit * 1.0;

    int x_leftindex = std::floor(x_FloatIndex);
    int y_leftindex = std::floor(y_FloatIndex);
    int z_leftindex = std::floor(z_FloatIndex);

    if(x_leftindex<0) {
      x_leftindex = 0;
    } else {
      if(x_leftindex>GridNo-2) {
        x_leftindex = GridNo-2;
      }
    }

    if(y_leftindex<0) {
      y_leftindex = 0;
    } else {
      if(y_leftindex>GridNo-2) {
        y_leftindex = GridNo-2;
      }
    }

    if(z_leftindex<0) {
      z_leftindex = 0;
    } else {
      if(z_leftindex>GridNo-2) {
        z_leftindex = GridNo-2;
      }
    }

    int x_rightindex = x_leftindex + 1;
    int y_rightindex = y_leftindex + 1;
    int z_rightindex = z_leftindex + 1;

    double valA = SDFTensor(x_leftindex, y_leftindex, z_leftindex);
    double valB = SDFTensor(x_rightindex, y_leftindex, z_leftindex);
    double valC = SDFTensor(x_rightindex, y_rightindex, z_leftindex);
    double valD = SDFTensor(x_leftindex, y_rightindex, z_leftindex);

    double valE = SDFTensor(x_leftindex, y_leftindex, z_rightindex);
    double valF = SDFTensor(x_rightindex, y_leftindex, z_rightindex);
    double valG = SDFTensor(x_rightindex, y_rightindex, z_rightindex);
    double valH = SDFTensor(x_leftindex, y_rightindex, z_rightindex);

    // Since this is a tri-linear interpolation, there are three ways to do the interpolation

    /*
      Jacobian matrix in the x-direction
    */

    // Cut the point with a plane orthgonal to z axis
    double valMAE =  (z_FloatIndex - z_leftindex*1.0) * (valE - valA) + valA;
    double valMBF =  (z_FloatIndex - z_leftindex*1.0) * (valF - valB) + valB;
    double valMDH =  (z_FloatIndex - z_leftindex*1.0) * (valH - valD) + valD;
    double valMCG =  (z_FloatIndex - z_leftindex*1.0) * (valG - valC) + valC;
    // Cut the point with a plane orthgonal to y axis
    double valMAEDH = (y_FloatIndex - y_leftindex*1.0) * (valMDH - valMAE) + valMAE;
    double valMBFCG = (y_FloatIndex - y_leftindex*1.0) * (valMCG - valMBF) + valMBF;
    // The values at the edge give the jacobian to x
    double JacDistTo_x = (valMBFCG - valMAEDH)/Envi_x_unit;

    /*
      Jacobian matrix in the y-direction
    */
    double valMABEF = (x_FloatIndex - x_leftindex*1.0) * (valMBF - valMAE) + valMAE;
    double valMDHCG = (x_FloatIndex - x_leftindex*1.0) * (valMCG - valMDH) + valMDH;

    double JacDistTo_y = (valMBFCG - valMAEDH)/Envi_y_unit;

    /*
      Jacobian matrix in the z-direction
    */
    // Cut the point with a plane orthgonal to x axis
    double valMAB = (x_FloatIndex - x_leftindex*1.0) * (valB - valA) + valA;
    double valMDC = (x_FloatIndex - x_leftindex*1.0) * (valC - valD) + valD;
    double valMEF = (x_FloatIndex - x_leftindex*1.0) * (valF - valE) + valE;
    double valMHG = (x_FloatIndex - x_leftindex*1.0) * (valG - valH) + valH;
    // Cut the point with a plane orthgonal to y axis
    double valMABDC = (y_FloatIndex - y_leftindex*1.0) * (valMDC - valMAB) + valMAB;
    double valMEFHG = (y_FloatIndex - y_leftindex*1.0) * (valMHG - valMEF) + valMHG;
    // The values at the edge give the jacobian to z
    double JacDistTo_z = (valMEFHG - valMABDC)/Envi_z_unit;

    Vector3 RawNormal(JacDistTo_x, JacDistTo_y, JacDistTo_z);
    RawNormal.setNormalized(RawNormal);
    return RawNormal;
  }
  Eigen::Tensor<double, 3> SDFTensor;
  double Envi_x_min, Envi_x_max;
  double Envi_y_min, Envi_y_max;
  double Envi_z_min, Envi_z_max;
  double Envi_x_unit, Envi_y_unit, Envi_z_unit;
  double Envi_x_length, Envi_y_length, Envi_z_length;
  int GridNo;
};

struct RMPoint {
  // This struct is used to save the information of a point in Reachability Map.
  RMPoint();
  RMPoint(const double & r, const Vector3 & Pos) {
    // Constructor
    Radius = r;
    Position = Pos;
    Direction = Position;
    double PositionLength = std::sqrt(Position.x * Position.x + Position.y * Position.y + Position.z * Position.z);
    Direction.x = Direction.x/PositionLength;
    Direction.y = Direction.y/PositionLength;
    Direction.z = Direction.z/PositionLength;
  }
  double Radius;
  Vector3 Position;
  Vector3 Direction;
};

struct ReachabilityMap {
  // This struct is used to save the reachability map
  ReachabilityMap();
  ReachabilityMap(const std::map<int, std::vector<RMPoint>> & _RMLayers):RMLayers(_RMLayers){};
  void ReachabilityMapPara(const double _MaxRadius, const int _LayerNumber, const int _PointNumberOnInner, const double _LayerDiff, const double _MinRadius){
    MaxRadius = _MaxRadius;
    LayerNumber = _LayerNumber;
    PointNumberOnInner = _PointNumberOnInner;
    LayerDiff = _LayerDiff;
    MinRadius = _MinRadius;
  }
  std::vector<Vector3> IdealReachablePointsFinder(const Robot & SimRobot, const int & LinkInfoIndex){
    std::vector<Vector3> ReachablePoints;
    ReachablePoints.reserve(TotalPoint);
    double PivotalLinkIndex = EndEffectorPivotalIndex[LinkInfoIndex];
    Vector3 RefPoint, ZeroPos(0.0, 0.0, 0.0);
    SimRobot.GetWorldPosition(ZeroPos, PivotalLinkIndex, RefPoint);
    for (int i = 0; i < LayerNumber; i++){
      std::vector<RMPoint> RMLayer_i = RMLayers[i];
      for (int j = 0; j < RMLayer_i.size(); j++){
        Vector3 RMPointPos = RMLayer_i[j].Position + RefPoint;
        ReachablePoints.push_back(RMPointPos);
      }
    }
    return ReachablePoints;
  }
  std::vector<Vector3> ReachablePointsFinder(const Robot & SimRobot, const int & LinkInfoIndex, SignedDistanceFieldInfo & SDFInfo, const Vector3 & COMVel){
    double Radius = EndEffectorRadius[LinkInfoIndex];
    double PivotalLinkIndex = EndEffectorPivotalIndex[LinkInfoIndex];

    Vector3 RefPoint, ZeroPos(0.0, 0.0, 0.0);
    SimRobot.GetWorldPosition(ZeroPos, PivotalLinkIndex, RefPoint);
    int ReachablePointNo = 0;
    return ReachablePointsGene(RefPoint, Radius, SDFInfo, ReachablePointNo, COMVel);
  }
  std::vector<Vector3> ReachablePointsFinder(const Robot & SimRobot, const int & LinkInfoIndex, SignedDistanceFieldInfo & SDFInfo){
    double Radius = EndEffectorRadius[LinkInfoIndex];
    double PivotalLinkIndex = EndEffectorPivotalIndex[LinkInfoIndex];

    Vector3 RefPoint, ZeroPos(0.0, 0.0, 0.0);
    SimRobot.GetWorldPosition(ZeroPos, PivotalLinkIndex, RefPoint);
    Vector3 COMVel(0.0, 0.0, 0.0);
    int ReachablePointNo = 0;
    return ReachablePointsGene(RefPoint, Radius, SDFInfo, ReachablePointNo, COMVel);
  }
  std::vector<Vector3> ReachablePointsGene(const Vector3 & RefPoint, const double & Radius, SignedDistanceFieldInfo & SDFInfo, int & ReachablePointNo, const Vector3 & COMVel){
    // Here MinRadius is used for comparison while Radius is the maximum allowed radius of robot's end effector.
    // This function is used to get presumably active ReachablePointsGene() from all sampled points.
    const double DisTol = 0.01;        // 1cm as a signed distance tolerance.
    std::vector<Vector3> ReachablePoints;
    ReachablePoints.reserve(TotalPoint);
    ReachablePointNo = 0;
    int LayerIndex = 0;
    double LayerRadius = MinRadius;
    if(Radius>MaxRadius){
      LayerIndex = LayerNumber-1;
    } else{
      while (LayerRadius<Radius){
        LayerRadius+=LayerDiff;
        LayerIndex++;
      }
      LayerRadius-=LayerDiff;
    }

    for (int i = 0; i < LayerIndex; i++){
      std::vector<RMPoint> RMLayer_i = RMLayers[i];
      for (int j = 0; j < RMLayer_i.size(); j++){
        Vector3 RMPointPos = RMLayer_i[j].Position + RefPoint;
        double CurrentDist = SDFInfo.SignedDistance(RMPointPos);
        // if(CurrentDist*CurrentDist<DisTol*DisTol){
          if((CurrentDist>0.0)&&(CurrentDist<DisTol)){
          ReachablePoints.push_back(RMPointPos);
          ReachablePointNo++;
        }
      }
    }
    return ReachablePoints;
  }
  std::vector<Vector3> ContactFreePointsFinder(const double & radius, const std::vector<Vector3> & ReachablePoints,const std::vector<std::pair<Vector3, double>> & ContactFreeInfo){
    // This function can only be called after ReachablePointsFinder() to reduce the extra point further.
    std::vector<Vector3> ContactFreePoints;
    ContactFreePoints.reserve(ReachablePoints.size());
    int ContactFreeNo = 0;
    for (int i = 0; i < ReachablePoints.size(); i++){
      Vector3 ReachablePoint = ReachablePoints[i];
      bool ContactFreeFlag = true;
      int ContactFreeInfoIndex = 0;
      while (ContactFreeInfoIndex<ContactFreeInfo.size()){
        Vector3 RefPoint = ContactFreeInfo[ContactFreeInfoIndex].first;
        double Radius = ContactFreeInfo[ContactFreeInfoIndex].second;
        Vector3 PosDiff = ReachablePoint - RefPoint;
        double PosDiffDis = sqrt(PosDiff.x * PosDiff.x + PosDiff.y * PosDiff.y + PosDiff.z * PosDiff.z);
        if(PosDiffDis<=(Radius + radius)){
          ContactFreeFlag = false;
          break;
        }
        ContactFreeInfoIndex++;
      }
      if(ContactFreeFlag){
        ContactFreePoints.push_back(ReachablePoints[i]);
        ContactFreeNo++;
      }
    }
    return ContactFreePoints;
  };
  std::map<int, std::vector<RMPoint>> RMLayers;
  std::vector<double> EndEffectorCollisionRadius;
  std::vector<double> EndEffectorRadius;
  std::vector<int> EndEffectorLinkIndex;              //
  std::vector<int> EndEffectorPivotalIndex;
  std::map<int, std::vector<int>> EndEffectorLink2Pivotal;       // This map saves intermediate joint from End Effector Joint to Pivotal Joint.
  std::map<int, int> EndEffectorIndices;
  double MaxRadius;
  int LayerNumber;
  int PointNumberOnInner;
  double LayerDiff;
  double MinRadius;
  int TotalPoint;
};

struct SelfLinkGeoInfo{
  SelfLinkGeoInfo(){
    SelfLinkGeoFlag = false;
  };
  SelfLinkGeoInfo(const Robot & SimRobot, const std::map<int, std::vector<int>> & EndEffectorLink2Pivotal, const std::vector<int> & SelfCollisionFreeLink){
    for (int i = 5; i < SimRobot.q.size(); i++){
      AABB3D AABB3D_i = SimRobot.geometry[i]->GetAABBTight();
      Frame3D LinkTransforms_i = SimRobot.links[i].T_World;
      LinkBBs.push_back(AABB3D_i);
      LinkTransforms.push_back(LinkTransforms_i);
    }

    for (int i = 0; i < EndEffectorLink2Pivotal.size(); i++){
      std::vector<int> EndLink2PivIndices = EndEffectorLink2Pivotal.at(i);
      for (int j = 5; j < SimRobot.q.size(); j++){
        if(std::find(EndLink2PivIndices.begin(), EndLink2PivIndices.end(), j) == EndLink2PivIndices.end()){
          if(std::find(SelfCollisionFreeLink.begin(), SelfCollisionFreeLink.end(), j) == SelfCollisionFreeLink.end()){
            SelfCollisionLinkMap[i].push_back(j-5);
          }
        }
      }
    }
    SelfLinkGeoFlag = true;
  };
  void LinkBBsUpdate(const Robot& SimRobot){
    for (int i = 5; i < SimRobot.q.size(); i++){
      AABB3D AABB3D_i = SimRobot.geometry[i]->GetAABBTight();
      LinkBBs[i-5] = AABB3D_i;
      Frame3D LinkTransforms_i = SimRobot.links[i].T_World;
      LinkTransforms[i-5] = LinkTransforms_i;
    }
  }
  double SingleLinkDist(const int & LinkCountIndex, const Vector3 & GlobalPoint){
    return LinkBBs[LinkCountIndex].signedDistance(GlobalPoint);
  }

  void SingleLinkDistNGrad(const int & LinkCountIndex, const Vector3 & GlobalPoint, double & Dist, Vector3 & DistGrad){
    const int GridNo = 100;
    double dx = LinkBBs[LinkCountIndex].size().x/(1.0 * GridNo);
    double dy = LinkBBs[LinkCountIndex].size().y/(1.0 * GridNo);
    double dz = LinkBBs[LinkCountIndex].size().z/(1.0 * GridNo);

    Dist = LinkBBs[LinkCountIndex].signedDistance(GlobalPoint);

    Vector3 GlobalPointx = GlobalPoint;
    GlobalPointx.x += dx;
    double Distx = LinkBBs[LinkCountIndex].signedDistance(GlobalPointx);

    Vector3 GlobalPointy = GlobalPoint;
    GlobalPointy.y += dy;
    double Disty = LinkBBs[LinkCountIndex].signedDistance(GlobalPointy);

    Vector3 GlobalPointz = GlobalPoint;
    GlobalPointz.z += dz;
    double Distz = LinkBBs[LinkCountIndex].signedDistance(GlobalPointz);

    DistGrad.x = (Distx - Dist)/dx;
    DistGrad.y = (Disty - Dist)/dy;
    DistGrad.z = (Distz - Dist)/dz;
    DistGrad.getNormalized(DistGrad);
  }

  std::vector<Vector3> BBVertices(const int & BBIndex){
    std::vector<Vector3> Vertices(8);
    // This function calculates the vertices for current bounding box.
    AABB3D CurBB = LinkBBs[BBIndex];
    Vector3 CurBBSize = CurBB.size();
    RigidTransform CurBBtoWorld = LinkTransforms[BBIndex];

    Vector3 bmin = CurBB.bmin;
    Vector3 bmax = CurBB.bmax;

    Vector3 bmin2bmax = bmax - bmin;

    Vector3 xAxis(CurBBtoWorld.R.data[0][0], CurBBtoWorld.R.data[0][1], CurBBtoWorld.R.data[0][2]);
    Vector3 yAxis(CurBBtoWorld.R.data[1][0], CurBBtoWorld.R.data[1][1], CurBBtoWorld.R.data[1][2]);
    Vector3 zAxis(CurBBtoWorld.R.data[2][0], CurBBtoWorld.R.data[2][1], CurBBtoWorld.R.data[2][2]);

    double xbmin2bmax = xAxis.dot(bmin2bmax);
    double ybmin2bmax = yAxis.dot(bmin2bmax);
    double zbmin2bmax = zAxis.dot(bmin2bmax);

    // Now we can have eight points explicitly.
    Vertices[0] = bmin;
    Vertices[1] = bmin + xbmin2bmax * xAxis;
    Vertices[2] = bmin + zbmin2bmax * zAxis;
    Vertices[3] = bmin + xbmin2bmax * xAxis + zbmin2bmax * zAxis;

    Vertices[4] = bmin + ybmin2bmax * yAxis;
    Vertices[5] = bmin + ybmin2bmax * yAxis + xbmin2bmax * xAxis;
    Vertices[6] = bmin + ybmin2bmax * yAxis + zbmin2bmax * zAxis;
    Vertices[7] = bmin + ybmin2bmax * yAxis + xbmin2bmax * xAxis + zbmin2bmax * zAxis;
    return Vertices;
  }

  void Vector3Writer(const std::vector<Vector3> & ContactPoints, const std::string &ContactPointFileName){
    if(ContactPoints.size()==0){
      std::cerr<< " ContactPoints has zero element!\n" << endl;
    }
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

  double SelfCollisionDist(const int & LinkIndex, const Vector3 & GlobalPoint){
    const int ActLinkNo = SelfCollisionLinkMap[LinkIndex].size();
    std::vector<double> DistVec;
    DistVec.reserve(ActLinkNo);
    for (int i = 0; i < ActLinkNo; i++){
      int SelfLinkIndex = SelfCollisionLinkMap[LinkIndex][i];
      double Dist_i = LinkBBs[SelfLinkIndex].signedDistance(GlobalPoint);
      DistVec.push_back(Dist_i);
    }
    double Dist = *std::min_element(DistVec.begin(), DistVec.end());
    return Dist;
  }

  void SelfCollisionDistNGrad(const int & LinkIndex, const Vector3 & GlobalPoint, double & Dist, Vector3 & Grad){
    // This function is used to calculate robot's self-collision distance given a point
    const int ActLinkNo = SelfCollisionLinkMap[LinkIndex].size();
    std::vector<double> DistVec;
    DistVec.reserve(ActLinkNo);
    std::vector<Vector3> GradVec;
    GradVec.reserve(ActLinkNo);
    std::vector<double> DistWeights;
    DistWeights.reserve(ActLinkNo);
    double Dist_i;
    Vector3 Grad_i;
    for (int i = 0; i < ActLinkNo; i++){
      int SelfLinkIndex = SelfCollisionLinkMap[LinkIndex][i];
      SingleLinkDistNGrad(SelfLinkIndex, GlobalPoint, Dist_i, Grad_i);
      DistVec.push_back(Dist_i);
      GradVec.push_back(Grad_i);
    }
    Dist = *std::min_element(DistVec.begin(), DistVec.end());
    double Scale = abs(Dist);
    for (int i = 0; i < ActLinkNo; i++){
      DistWeights.push_back(exp(-1.0 * DistVec[i]/Scale));
    }
    // Set its value to be zero!
    Grad.x = 0.0;
    Grad.y = 0.0;
    Grad.z = 0.0;
    for (int i = 0; i < ActLinkNo; i++)
      Grad+=DistWeights[i] * GradVec[i];
    Grad.getNormalized(Grad);
  }
  bool SelfLinkGeoFlag;
  std::vector<AABB3D> LinkBBs;
  std::vector<RigidTransform> LinkTransforms;
  std::map<int, std::vector<int>> SelfCollisionLinkMap;       // This map saves intermediate joint from End Effector Joint to Pivotal Joint.
};

struct DataRecorderInfo{
  DataRecorderInfo(){
    PlanStageIndex = -1;
    LinkNo = -1;
  };
  void setPlanStageIndexNLinkNo(const int & _PlanStageIndex, const int & _LinkNo){
    PlanStageIndex = _PlanStageIndex;
    LinkNo = _LinkNo;
  }
  void setRCSData( const std::vector<Vector3> & _ReachableContacts,
                    const std::vector<Vector3> & _CollisionFreeContacts,
                    const std::vector<Vector3> & _SupportiveContacts){
                      ReachableContacts = _ReachableContacts;
                      CollisionFreeContacts = _CollisionFreeContacts;
                      SupportiveContacts = _SupportiveContacts;
  }
  void setCCSData(  const std::vector<Vector3> & _CandidateContacts,
                    const std::vector<Vector3> & _CandidateContactWeights,
                    const std::vector<Vector3> & _SelectedContacts){
                      CandidateContacts = _CandidateContacts;
                      CandidateContactWeights = _CandidateContactWeights;
                      SelectedContacts = _SelectedContacts;
  }
  void setPathWaypoints(const std::vector<Vector3> & _PathWaypoints){ PathWaypoints = _PathWaypoints; }
  void setTrajs(const LinearPath & _PlannedConfigTraj, const LinearPath & _EndEffectorTraj){
    PlannedConfigTraj = _PlannedConfigTraj;
    EndEffectorTraj = _EndEffectorTraj;
  }
  void Vector3Writer(const std::vector<Vector3> & ContactPoints, const std::string & ContactPointFileName){
    if(ContactPoints.size() ==0) return;
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
  void Write2File(const string & CurrentCasePath){
    // This function will only be called if planning is successful!
    const string InnerPath = CurrentCasePath + std::to_string(PlanStageIndex) + "_" + std::to_string(LinkNo) + "_";
    Vector3Writer(ReachableContacts,      InnerPath + "ReachableContacts");
    Vector3Writer(CollisionFreeContacts,  InnerPath + "CollisionFreeContacts");
    Vector3Writer(SupportiveContacts,     InnerPath + "SupportiveContacts");
    Vector3Writer(CandidateContacts,      InnerPath + "CandidateContacts");
    Vector3Writer(CandidateContactWeights,InnerPath + "CandidateContactWeights");
    Vector3Writer(SelectedContacts,       InnerPath + "SelectedContacts");
    Vector3Writer(PathWaypoints,          InnerPath + "PathWaypoints");

    // Write these two trajectories into files.
    std::ofstream PlannedConfigTrajFile;
    const string PlannedConfigTrajName = InnerPath + "PlannedConfigTraj.path";
    PlannedConfigTrajFile.open (PlannedConfigTrajName.c_str());
    PlannedConfigTraj.Save(PlannedConfigTrajFile);
    PlannedConfigTrajFile.close();

    std::ofstream EndEffectorTrajFile;
    const string EndEffectorTrajName = InnerPath + "EndEffectorTraj.path";
    EndEffectorTrajFile.open(EndEffectorTrajName.c_str());
    EndEffectorTraj.Save(EndEffectorTrajFile);
    EndEffectorTrajFile.close();
  }
  std::vector<Vector3> ReachableContacts;
  std::vector<Vector3> CollisionFreeContacts;
  std::vector<Vector3> SupportiveContacts;
  std::vector<Vector3> CandidateContacts;
  std::vector<Vector3> CandidateContactWeights;
  std::vector<Vector3> SelectedContacts;
  std::vector<Vector3> PathWaypoints;

  LinearPath PlannedConfigTraj;
  LinearPath EndEffectorTraj;

  int PlanStageIndex;
  int LinkNo;
};

struct SimPara{
  SimPara(){
    ForceMax = -1.0;
    PushDuration = -1.0;
    DetectionWait = -1.0;
    TimeStep = -1.0;
    InitDuration = -1.0;
    TotalDuration = -1.0;
    ForwardDuration = -1.0;
    PhaseRatio = -1.0;
    PhaseTimeStep = 0.05;
    PlanStageIndex = -1;
    SimTime = -1.0;
    TransPathFeasiFlag = false;
    SwingLinkInfoIndex =-1;
    CurrentContactPos.setZero();
    TrajConfigOptFlag = false;
  };
  SimPara(const double & _ForceMax,
          const double & _PushDuration,
          const double & _DetectionWait,
          const double & _TimeStep,
          const double & _InitDuration,
          const double & _TotalDuration,
          const double & _ForwardDuration,
          const double & _PhaseRatio,
          const double & _PhaseTimeStep,
          const double & _ReductionRatio):    ForceMax(_ForceMax),
                                              PushDuration(_PushDuration),
                                              DetectionWait(_DetectionWait),
                                              TimeStep(_TimeStep),
                                              InitDuration(_InitDuration),
                                              TotalDuration(_TotalDuration),
                                              ForwardDuration(_ForwardDuration),
                                              PhaseRatio(_PhaseRatio),
                                              PhaseTimeStep(_PhaseTimeStep),
                                              ReductionRatio(_ReductionRatio){}
  void CurrentCasePathUpdate(const string _CurrentCasePath){
    CurrentCasePath = _CurrentCasePath;

    string fedge_aFile = CurrentCasePath + "EdgeATraj.txt";
    // const char *fedge_aFile_Name = fedge_aFile.c_str();
    string fedge_bFile = CurrentCasePath + "EdgeBTraj.txt";
    // const char *fedge_bFile_Name = fedge_bFile.c_str();
    string fEdgeCOMFile = CurrentCasePath + "EdgeCOMTraj.txt";
    // const char *fEdgeCOMFile_Name = fEdgeCOMFile.c_str();
    string fEdgexTrajFile = CurrentCasePath + "EdgexTraj.txt";
    // const char *fEdgexTrajFile_Name = fEdgexTrajFile.c_str();
    string fEdgeyTrajFile = CurrentCasePath + "EdgeyTraj.txt";
    // const char *fEdgeyTrajFile_Name = fEdgeyTrajFile.c_str();
    string fEdgezTrajFile = CurrentCasePath + "EdgezTraj.txt";
    // const char *fEdgezTrajFile_Name = fEdgezTrajFile.c_str();
    string fVertexTrajFile = CurrentCasePath + "EdgeVertexTraj.txt";

    EdgeFileNames.push_back(fedge_aFile);
    EdgeFileNames.push_back(fedge_bFile);
    EdgeFileNames.push_back(fEdgeCOMFile);
    EdgeFileNames.push_back(fEdgexTrajFile);
    EdgeFileNames.push_back(fEdgeyTrajFile);
    EdgeFileNames.push_back(fEdgezTrajFile);
    EdgeFileNames.push_back(fVertexTrajFile);

    FailureStateTrajStr =  CurrentCasePath + "FailureStateTraj.path";
    // const char *FailureStateTrajStr_Name = FailureStateTrajStr.c_str();
    CtrlStateTrajStr    =  CurrentCasePath + "CtrlStateTraj.path";
    // const char *CtrlStateTrajStr_Name = CtrlStateTrajStr.c_str();
    PlanStateTrajStr = CurrentCasePath + "PlanStateTraj.path";
    // const char *PlanStateTrajStr_Name = PlanStateTrajStr.c_str();
  }
  string getCurrentCasePath(){ return CurrentCasePath; }
  void setImpulseForceMax(const Vector3 & ImpulseDirection){ ImpulseForceMax = ForceMax * ImpulseDirection; }
  void setPlanStageIndex(const int & _PlanStageIndex) {PlanStageIndex = _PlanStageIndex; }
  int  getPlanStageIndex(){ return PlanStageIndex; }
  void setPlanEndEffectorIndex( const int & _PlanEndEffectorIndex) { PlanEndEffectorIndex = _PlanEndEffectorIndex; }
  int  getPlanEndEffectorIndex(){ return PlanEndEffectorIndex; }
  void setSimTime(const double & _SimTime) { SimTime = _SimTime; }
  double getSimTime() { return SimTime; }
  void setContactInit(const Vector3 _ContactInit){ ContactInit = _ContactInit; }
  Vector3 getContactInit() { return ContactInit; }
  void setContactGoal(const Vector3 _ContactGoal){ ContactGoal = _ContactGoal;}
  Vector3 getContactGoal(){ return ContactGoal;}
  void setDirectionInit(const Vector3 & _DirectionInit ){ DirectionInit = _DirectionInit; }
  void setDirectionGoal(const Vector3 & _DirectionGoal ){ DirectionGoal = _DirectionGoal; }
  Vector3 getGoalDirection() { return DirectionGoal; }
  void setTransPathFeasiFlag(const bool & _TransPathFeasiFlag){ TransPathFeasiFlag = _TransPathFeasiFlag; }
  bool getTransPathFeasiFlag(){ return TransPathFeasiFlag;}
  void setSwingLinkInfoIndex(const int _SwingLinkInfoIndex) {SwingLinkInfoIndex = _SwingLinkInfoIndex; }
  int  getSwingLinkInfoIndex() { return SwingLinkInfoIndex;}
  void setCurrentContactPos(const Vector3 & _CurrentContactPos) {CurrentContactPos = _CurrentContactPos; }
  Vector3 getCurrentContactPos(){ return CurrentContactPos; }
  void setTrajConfigOptFlag(const bool & _TrajConfigOptFlag){ TrajConfigOptFlag = _TrajConfigOptFlag;}
  bool getTrajConfigOptFlag() {return TrajConfigOptFlag;}
  void setFixedContactStatusInfo(const std::vector<ContactStatusInfo> & _FixedContactStatusInfo){ FixedContactStatusInfo =_FixedContactStatusInfo;}


  double  ForceMax;
  double  PushDuration;
  double  DetectionWait;
  double  TimeStep;
  double  InitDuration;
  double  TotalDuration;
  double  ForwardDuration;
  double  PhaseRatio;             // This ratio determines the boundary between acceleration and deceleration.
  double  PhaseTimeStep;
  double  ReductionRatio;
  int     PlanStageIndex;
  int     PlanEndEffectorIndex;    // This PlanEndEffectorIndex saves successful end effector for push recovery.
  double  SimTime;
  bool    TransPathFeasiFlag;
  int     SwingLinkInfoIndex;
  Vector3 CurrentContactPos;
  bool    TrajConfigOptFlag;

  DataRecorderInfo DataRecorderObj;
  Vector3 ImpulseForceMax;
  std::string CurrentCasePath;
  std::vector<string >EdgeFileNames;
  string FailureStateTrajStr, CtrlStateTrajStr, PlanStateTrajStr;
  Vector3 ContactInit, ContactGoal;
  Vector3 DirectionInit, DirectionGoal;
  std::vector<ContactStatusInfo> FixedContactStatusInfo;
};

struct ControlReferenceInfo{
  ControlReferenceInfo(){
    ReadyFlag = false;
    TouchDownFlag = false;
    ControlReferenceType = -1;
    SwingLinkInfoIndex = -1;           // Used for RobotLinkInfo
    ContactStatusInfoIndex = -1;
    GoalContactPos.setZero();
    GoalContactGrad.setZero();
    }
  void setSwingLinkInfoIndex(const int & _SwingLinkInfoIndex) {SwingLinkInfoIndex = _SwingLinkInfoIndex;}
  int  getSwingLinkInfoIndex() { return SwingLinkInfoIndex; }
  bool getReadyFlag(){ return ReadyFlag;}
  void setReadyFlag(const bool & _ReadyFlag ){ ReadyFlag = _ReadyFlag; }
  void setTouchDownFlag(const bool & _TouchDownFlag){ TouchDownFlag=_TouchDownFlag; }
  bool getTouchDownFlag(){ return TouchDownFlag; }
  int  getControlReferenceType(){ return ControlReferenceType; }
  void setControlReferenceType(const int &_ControlReferenceType) { ControlReferenceType = _ControlReferenceType; }
  void setGoalContactPosNGrad(const Vector3 & _GoalContactPos, const Vector3 & _GoalContactGrad){
    GoalContactPos = _GoalContactPos;
    GoalContactGrad = _GoalContactGrad;
  }
  Vector3 getGoalContactPos(){  return GoalContactPos; }
  Vector3 getGoalContactGrad(){ return GoalContactGrad;}

  void SetInitContactStatus(const std::vector<ContactStatusInfo> &_InitContactStatus){ InitContactStatus = _InitContactStatus; }
  std::vector<ContactStatusInfo> getInitContactStatus(){return InitContactStatus;}

  void SetGoalContactStatus(const std::vector<ContactStatusInfo> & _GoalContactStatus) {GoalContactStatus = _GoalContactStatus;}
  std::vector<ContactStatusInfo> getGoalContactStatus() {return GoalContactStatus;}

  void TrajectoryUpdate(const std::vector<double> & timeTraj, const std::vector<Config> & configTraj, const std::vector<Vector3> & endeffectorTraj){
    TimeTraj = timeTraj;
    PlannedConfigTraj = LinearPath(timeTraj, configTraj);
    std::vector<Vector> endeffectorPath;
    for (Vector3 EndEffectorPos: endeffectorTraj){
      Vector EndEffectorPosVec;
      EndEffectorPosVec.resize(3);
      EndEffectorPosVec[0] = EndEffectorPos[0];
      EndEffectorPosVec[1] = EndEffectorPos[1];
      EndEffectorPosVec[2] = EndEffectorPos[2];
      endeffectorPath.push_back(EndEffectorPosVec);
    }
    EndEffectorTraj = LinearPath(timeTraj, endeffectorPath);
  }
  void setWaitTime(const double & _WaitTime) { WaitTime = _WaitTime; }
  double getWaitTime() { return WaitTime; }
  void setTouchDownConfig(const std::vector<double> _TouchDownConfig){ TouchDownConfig = _TouchDownConfig; }
  std::vector<double> getTouchDownConfig(){ return TouchDownConfig;}

  bool    ReadyFlag;
  bool    TouchDownFlag;
  int     ControlReferenceType;
  int     SwingLinkInfoIndex;
  int     ContactStatusInfoIndex;
  double  WaitTime;

  Vector3 GoalContactPos;
  Vector3 GoalContactGrad;

  std::vector<double> TouchDownConfig;

  LinearPath PlannedConfigTraj;
  LinearPath EndEffectorTraj;

  std::vector<double> TimeTraj;
  std::vector<ContactStatusInfo> InitContactStatus;
  std::vector<ContactStatusInfo> GoalContactStatus;
};

struct FailureStateInfo{
  FailureStateInfo(){
    FailureTime = 0.0;
    FailureStateFlag = false;
  };

  Config getFailureStateConfig() { return FailureConfig; }
  Config getFailureStateVelocity(){ return FailureVelocity; }
  bool   getFailureStateFlag() { return FailureStateFlag; }
  void FailureStateUpdate(const double & _FailureTime, const Config & _FailureConfig, const Config & _FailureVelocity){
    FailureTime = _FailureTime;
    FailureConfig = _FailureConfig;
    FailureVelocity = _FailureVelocity;
    FailureStateFlag = true;
  };
  double  FailureTime;
  Config  FailureConfig;
  Config  FailureVelocity;
  bool    FailureStateFlag;
};

struct FacetInfo{
  FacetInfo(){
    FacetValidFlag = false;
  };
  void setFacetValidFlag(const bool & _FacetValidFlag){FacetValidFlag = _FacetValidFlag;}
  bool getFacetValidFlag(){ return FacetValidFlag;}
  void setFacetEdges(const std::vector<std::pair<Vector3, Vector3>> & _FacetEdges) { FacetEdges = _FacetEdges; }
  void setFacetNorm(const Vector3& _FacetNorm){ FacetNorm = _FacetNorm;}
  double ProjPoint2EdgeDist(const Vector3& _Point){
    std::vector<double> ProjPoint2Edge_vec(EdgeNorms.size());
    Vector3 Vertex2Point = _Point - FacetEdges[0].first;
    double Point2Facet = Vertex2Point.dot(FacetNorm);
    Vector3 Facet2Point = Point2Facet * FacetNorm;
    for (int i = 0; i < EdgeNorms.size(); i++){
      Vertex2Point = _Point - FacetEdges[i].first;
      Vector3 Vertex2ProjPoint = Vertex2Point - Facet2Point;
      double ProjPoint2Edge_i = Vertex2ProjPoint.dot(EdgeNorms[i]);
      ProjPoint2Edge_vec[i] = ProjPoint2Edge_i;
    }
    return *min_element(ProjPoint2Edge_vec.begin(), ProjPoint2Edge_vec.end());
  }
  std::vector<double> ProjPoint2EdgeDistVec(const Vector3& _Point){
    std::vector<double> ProjPoint2Edge_vec(EdgeNorms.size());
    Vector3 Vertex2Point = _Point - FacetEdges[0].first;
    double Point2Facet = Vertex2Point.dot(FacetNorm);
    Vector3 Facet2Point = Point2Facet * FacetNorm;
    for (int i = 0; i < EdgeNorms.size(); i++){
      Vertex2Point = _Point - FacetEdges[i].first;
      Vector3 Vertex2ProjPoint = Vertex2Point - Facet2Point;
      double ProjPoint2Edge_i = Vertex2ProjPoint.dot(EdgeNorms[i]);
      ProjPoint2Edge_vec[i] = ProjPoint2Edge_i;
    }
    return ProjPoint2Edge_vec;
  }
  void EdgesUpdate(){
    Edges.reserve(EdgeNorms.size());
    EdgesDirection.reserve(EdgeNorms.size());
    for (int i = 0; i < EdgeNorms.size(); i++){
      Vector3 Edge_i = FacetEdges[i].second - FacetEdges[i].first;
      Edges.push_back(Edge_i);
      Vector3 Edge_i_normalized;
      Edge_i.getNormalized(Edge_i_normalized);
      EdgesDirection.push_back(Edge_i_normalized);
    }
  }
  std::vector<std::pair<Vector3, Vector3>> FacetEdges;
  std::vector<Vector3> EdgeNorms;
  Vector3 FacetNorm;
  std::vector<Vector3> Edges;
  std::vector<Vector3> EdgesDirection;
  bool FacetValidFlag;
};

struct PIPInfo{
  // This struct saves the information of the projected inverted pendulum from the CoM to the edge of convex polytope
  PIPInfo(){
    L = 0.25;         // The reference bound range is [0.25, 0.85]
    Ldot = 0.0;
    theta = 0.0;
    thetadot = 0.0;
    g = 9.81;
    g_angle = 0.0;
    speed = -1.0;
    onFlag = false;
  }
  PIPInfo(double _L, double _Ldot, double _theta, double _thetadot, double _g, double _g_angle){
    L = _L;
    Ldot = _Ldot;
    theta = _theta;
    thetadot = _thetadot;
    g = _g;
    g_angle = _g_angle;
  }
  void setPrimeUnits(const Vector3 & x_prime_unit_,const Vector3 & y_prime_unit_,const Vector3 & z_prime_unit_){
    x_prime_unit = x_prime_unit_;
    y_prime_unit = y_prime_unit_;
    z_prime_unit = z_prime_unit_;
  }
  void setUnits(const Vector3 & x_unit_,const Vector3 & y_unit_,const Vector3 & z_unit_){
    x_unit = x_unit_;
    y_unit = y_unit_;
    z_unit = z_unit_;
  }
  void setEdgeAnB(const Vector3 & edge_a_, const Vector3 & edge_b_){
    edge_a = edge_a_;
    edge_b = edge_b_;
  }
  void setIntersection(const Vector3 & intersection_){ intersection = intersection_;}
  void setSpeed(const double & _speed ) {speed = _speed;}
  double getSpeed(){ return speed;}

  double  L, Ldot, theta, thetadot;
  double  g, g_angle;
  double  speed;                            // This value indicates the horizontal velocity.
  bool    onFlag;                           // Whether the origin is at intersection or not?!
  Vector3 x_prime_unit, y_prime_unit, z_prime_unit;
  Vector3 x_unit, y_unit, z_unit;
  Vector3 edge_a, edge_b;                     // The Edge points from edge_a to edge_b.
  Vector3 intersection;                     // The point where the COM intersects the edge.
};

struct ContactForm{
  ContactForm();
  ContactForm(      const std::vector<ContactStatusInfo> & _ContactStatusInfoObj,
                    const int & _SwingLinkInfoIndex,
                    const int & _ContactType):  FixedContactStatusInfo(_ContactStatusInfoObj),
                                                SwingLinkInfoIndex(_SwingLinkInfoIndex),
                                                ContactType(_ContactType){};
  std::vector<ContactStatusInfo> FixedContactStatusInfo;
  int SwingLinkInfoIndex;
  int ContactType;
};

struct InvertedPendulumInfo{
  InvertedPendulumInfo(){
    L = -1.0;
    g = -1.0;
    Theta = -1.0;
    Thetadot = -1.0;
  };
  InvertedPendulumInfo( const double & _L,
                        const double & _g,
                        const double & _Theta,
                        const double & _Thetadot,
                        const Vector3 & _COMPos,
                        const Vector3 & _COMVel): L(_L),
                                                  g(_g),
                                                  Theta(_Theta),
                                                  Thetadot(_Thetadot),
                                                  COMPos(_COMPos),
                                                  COMVel(_COMVel){};
  void setEdges(const Vector3 & _edge_a, const Vector3 & _edge_b){
    edge_a = _edge_a;
    edge_b = _edge_b;
  }
  Vector3 edge_a, edge_b;
  double L, g;
  double Theta;
  double Thetadot;
  Vector3 COMPos;
  Vector3 COMVel;
};

struct SplineInfo{
  SplineInfo(){
    sStart = -1.0;
    sEnd = -1.0;
  };

  SplineInfo( const double & _sStart, const double & _sEnd,
              const Vector3 &_a, const Vector3 &_b,
              const Vector3 &_c, const Vector3 &_d):  sStart(_sStart), sEnd(_sEnd),
                                                      a(_a), b(_b), c(_c), d(_d){};
  Vector3 SplinePosVector(const double & s){
    Vector3 Pos = a * s * s * s + b * s * s + c * s + d;
    return Pos;
  }

  Vector3 SplineVelVector(const double & s){
    Vector3 Vel = 3.0 * a * s *s + 2.0 * b * s + c;
    return Vel;
  }
  double sStart, sEnd;
  Vector3 a, b, c, d;
};

struct EndEffectorPathInfo{
  // This info saves robot's Spline Object.
  EndEffectorPathInfo(){};
  EndEffectorPathInfo(const std::vector<SplineLib::cSpline3> & _SplineObj){
    SplineObj = _SplineObj;
    SplineNumber = SplineObj.size();
    SplineIndLength.reserve(SplineNumber);
    SplineSumLength.reserve(SplineNumber);
    TotalLength = 0.0f;
    for (int i = 0; i < SplineNumber; i++){
        float len = Length(SplineObj[i], 0.01f);
        SplineIndLength.push_back(len);
        TotalLength += len;
        SplineSumLength.push_back(TotalLength);
    }
    SplineLib::Vec3f GoalPos = Position(SplineObj[SplineNumber - 1], 1.0);
    GoalContactPos.x = GoalPos.x;
    GoalContactPos.y = GoalPos.y;
    GoalContactPos.z = GoalPos.z;
  }
  void s2Pos(const double & s, Vector3 & Pos){
    // Here this function is used to get the corresponding position given the total s.
    double sLength = s * TotalLength;
    // Enumerate SplineSumLength to get the segment index
    if(sLength>TotalLength) sLength = TotalLength;
    if(s<0) sLength = 0.0;
    int SegmentNo = 0;
    while (SegmentNo<SplineNumber-1){
      if(SplineSumLength[SegmentNo]>=sLength)
        break;
      SegmentNo++;
    }
    double ResLength = 0.0;
    switch (SegmentNo){
      case 0:
      {
        ResLength = sLength;
      }
      break;
      default:
      {
        ResLength = sLength - SplineSumLength[SegmentNo-1];
      }
      break;
    }
    float sLength_i = ResLength/SplineIndLength[SegmentNo];
    SplineLib::Vec3f ps = Position(SplineObj[SegmentNo], sLength_i);
    Pos.x = ps.x;
    Pos.y = ps.y;
    Pos.z = ps.z;
  }
  void PosNTang(const double & s, Vector3 & Pos, Vector3 & Tang)
  {
    // Here this function is used to get the corresponding position given the total s.
    double sLength = s * TotalLength;
    // Enumerate SplineSumLength to get the segment index
    if(sLength>TotalLength) sLength = TotalLength;
    if(s<0) sLength = 0.0;
    int SegmentNo = 0;
    while (SegmentNo<SplineNumber-1)
    {
      if(SplineSumLength[SegmentNo]>=sLength)
      {
        break;
      }
      SegmentNo++;
    }

    double ResLength = 0.0;
    switch (SegmentNo)
    {
      case 0:
      {
        ResLength = sLength;
      }
      break;
      default:
      {
        ResLength = sLength - SplineSumLength[SegmentNo - 1];
      }
      break;
    }

    float sLength_i = ResLength/SplineIndLength[SegmentNo];

    SplineLib::Vec3f ps = Position(SplineObj[SegmentNo], sLength_i);
    SplineLib::Vec3f vs = Velocity(SplineObj[SegmentNo], sLength_i);

    // This part is used for debugging.
    int SplineIndex;
    SplineLib::cSpline3 splines[SplineNumber];
    for (int i = 0; i < SplineNumber; i++)  splines[i] = SplineObj[i];
    float s_i = FindClosestPoint(ps, SplineNumber, splines, &SplineIndex);

    Pos.x = ps.x;
    Pos.y = ps.y;
    Pos.z = ps.z;

    Tang.x = vs.x;
    Tang.y = vs.y;
    Tang.z = vs.z;

    // double TangMag = sqrt(Tang.x * Tang.x + Tang.y * Tang.y + Tang.z * Tang.z);
    //
    // Tang.x = Tang.x/TangMag;
    // Tang.y = Tang.y/TangMag;
    // Tang.z = Tang.z/TangMag;
  }

  double Pos2s(const Vector3 & Pos){
    // This function is used to compute the corresponding sTotal from Euclidean position.
    SplineLib::cSpline3 splines[SplineNumber];
    for (int i = 0; i < SplineNumber; i++)
    {
      splines[i] = SplineObj[i];
    }
    SplineLib::Vec3f qp(Pos.x, Pos.y, Pos.z);
    int SplineIndex;
    float s_i = FindClosestPoint(qp, SplineNumber, splines, &SplineIndex);
    double s = 0.0;
    switch (SplineIndex)
    {
      case 0:
      {
        double sLength = s_i * SplineIndLength[0];
        s = sLength/TotalLength;
      }
      break;
      default:
      {
        double sLength = s_i * SplineIndLength[SplineIndex] + SplineSumLength[SplineIndex-1];
        s = sLength/TotalLength;
      }
      break;
    }
    return s;
  }
  std::vector<SplineLib::cSpline3> SplineObj;
  std::vector<double> SplineIndLength;    // This saves the individual length for each piecewise curve.
  std::vector<double> SplineSumLength;    // This saves the accumulated length from the starting point along the curve.
  double TotalLength;
  int SplineNumber;
  Vector3 GoalContactPos;                 // This saves the goal contact position.
};

#endif
