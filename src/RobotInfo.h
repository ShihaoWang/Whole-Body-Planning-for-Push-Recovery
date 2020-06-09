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

struct LinkInfo
{
  LinkInfo(){LinkIndex = -1;}
  LinkInfo(const int &link_index){ LinkIndex = link_index; }
  void AddLocalConact(const Vector3& local_contact){ LocalContacts.push_back(local_contact); }
  void AvgContactUpdate()
  {
    switch (LocalContacts.size())
    {
      case 0:
      {
        throw std::invalid_argument( "LocalContacts should have been initialized!" );
      }
      break;
      default:
      {
        Vector3 SumLocalContacts(0.0, 0.0, 0.0);
        for (int i = 0; i < LocalContacts.size(); i++)
        {
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

struct ContactStatusInfo
{
  // This struct saves the information of the contact status of each link end effector
  ContactStatusInfo(){ LinkIndex = -1;}
  ContactStatusInfo(const int & link_index){ LinkIndex = link_index; }
  void AddLocalConactStatus(const int & _contactstatus){ LocalContactStatus.push_back(_contactstatus); }
  void StatusSwitch(const int & Val)
  {
    for(int i = 0; i<LocalContactStatus.size(); i++)
    {
      LocalContactStatus[i] = Val;
    }
  }
  int LinkIndex;
  std::vector<int> LocalContactStatus;
};

struct TerrainInfo
{
  // Each terrain has a TerrainInfo struct
  TerrainInfo(){num_tris = -1;}
  int num_tris;
  std::vector<Tri> Tris;
  std::vector<Vector3> TriNormals;
  std::vector<int> Indices;
};

struct SignedDistanceFieldInfo
{
  // This struct is used to save the information of the signed distance field
  SignedDistanceFieldInfo()
  {
    Envi_x_min = 0;           Envi_x_max = 0;
    Envi_y_min = 0;           Envi_y_max = 0;
    Envi_z_min = 0;           Envi_z_max = 0;
    Envi_x_unit = 0;          Envi_y_unit = 0;          Envi_z_unit = 0;
    Envi_x_length = 0;        Envi_y_length = 0;        Envi_z_length = 0;
    GridNo = 0;
  }
  SignedDistanceFieldInfo(const Eigen::Tensor<double, 3>& _SDFTensor, const std::vector<double> &_SDFSpecs)
  {
    SDFTensor = _SDFTensor;
    Envi_x_min = _SDFSpecs[0];          Envi_x_max = _SDFSpecs[1];
    Envi_y_min = _SDFSpecs[2];          Envi_y_max = _SDFSpecs[3];
    Envi_z_min = _SDFSpecs[4];          Envi_z_max = _SDFSpecs[5];
    Envi_x_unit = _SDFSpecs[6];         Envi_y_unit = _SDFSpecs[7];         Envi_z_unit = _SDFSpecs[8];
    Envi_x_length = _SDFSpecs[9];       Envi_y_length = _SDFSpecs[10];      Envi_z_length = _SDFSpecs[11];
    GridNo = (int)_SDFSpecs[12];
  }
  double SignedDistance(const Vector3 &Point) const
  {
    // This function is used to compute the distance from a 3D point to the environment terrain
    // The first job is to figure out the nearest neighbours of the Points

    double x_FloatIndex = (Point.x - Envi_x_min)/Envi_x_unit * 1.0;
    double y_FloatIndex = (Point.y - Envi_y_min)/Envi_y_unit * 1.0;
    double z_FloatIndex = (Point.z - Envi_z_min)/Envi_z_unit * 1.0;

    int x_leftindex = std::floor(x_FloatIndex);
    int y_leftindex = std::floor(y_FloatIndex);
    int z_leftindex = std::floor(z_FloatIndex);

    if(x_leftindex<0)
    {
      x_leftindex = 0;
    }
    else
    {
      if(x_leftindex>GridNo-2)
      {
        x_leftindex = GridNo-2;
      }
    }

    if(y_leftindex<0)
    {
      y_leftindex = 0;
    }
    else
    {
      if(y_leftindex>GridNo-2)
      {
        y_leftindex = GridNo-2;
      }
    }

    if(z_leftindex<0)
    {
      z_leftindex = 0;
    }
    else
    {
      if(z_leftindex>GridNo-2)
      {
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
  Vector3 SignedDistanceNormal(const Vector3 &Point) const
  {
    // This function is used to calculate the (1 x 3) Jacobian matrix given the current Position
    // This function is used to compute the distance from a 3D point to the environment terrain
    // The first job is to figure out the nearest neighbours of the Points

    double x_FloatIndex = (Point.x - Envi_x_min)/Envi_x_unit * 1.0;
    double y_FloatIndex = (Point.y - Envi_y_min)/Envi_y_unit * 1.0;
    double z_FloatIndex = (Point.z - Envi_z_min)/Envi_z_unit * 1.0;

    int x_leftindex = std::floor(x_FloatIndex);
    int y_leftindex = std::floor(y_FloatIndex);
    int z_leftindex = std::floor(z_FloatIndex);

    if(x_leftindex<0)
    {
      x_leftindex = 0;
    }
    else
    {
      if(x_leftindex>GridNo-2)
      {
        x_leftindex = GridNo-2;
      }
    }

    if(y_leftindex<0)
    {
      y_leftindex = 0;
    }
    else
    {
      if(y_leftindex>GridNo-2)
      {
        y_leftindex = GridNo-2;
      }
    }

    if(z_leftindex<0)
    {
      z_leftindex = 0;
    }
    else
    {
      if(z_leftindex>GridNo-2)
      {
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

struct RMPoint
{
  // This struct is used to save the information of a point in Reachability Map.
  RMPoint()
  {
  }
  RMPoint(const double & r, const Vector3 & Pos)
  {
    // Constructor
    Radius = r;
    Position = Pos;
    Direction = Position;
    double PositionLength = sqrt(Position.x * Position.x + Position.y * Position.y + Position.z * Position.z);
    Direction.x = Direction.x/PositionLength;
    Direction.y = Direction.y/PositionLength;
    Direction.z = Direction.z/PositionLength;
  }
  double Radius;
  Vector3 Position;
  Vector3 Direction;
};

struct ReachabilityMap
{
  // This struct is used to save the reachability map
  ReachabilityMap()
  {}
  ReachabilityMap(const std::map<int, std::vector<RMPoint>> & _RMLayers)
  {
    RMLayers = _RMLayers;
  };
  void ReachabilityMapPara(const double & _MaxRadius, const int & _LayerNumber, const int & _PointNumberOnInner, const double & _LayerDiff, const double & _MinRadius)
  {
    MaxRadius = _MaxRadius;
    LayerNumber = _LayerNumber;
    PointNumberOnInner = _PointNumberOnInner;
    LayerDiff = _LayerDiff;
    MinRadius = _MinRadius;
  }
  std::vector<Vector3> IdealReachablePointsFinder(const Robot & SimRobot, const int & LinkInfoIndex)
  {
    std::vector<Vector3> ReachablePoints;
    ReachablePoints.reserve(TotalPoint);
    double PivotalLinkIndex = EndEffectorPivotalIndex[LinkInfoIndex];
    Vector3 RefPoint, ZeroPos(0.0, 0.0, 0.0);
    SimRobot.GetWorldPosition(ZeroPos, PivotalLinkIndex, RefPoint);
    for (int i = 0; i < LayerNumber; i++)
    {
      std::vector<RMPoint> RMLayer_i = RMLayers[i];
      for (int j = 0; j < RMLayer_i.size(); j++)
      {
        Vector3 RMPointPos = RMLayer_i[j].Position + RefPoint;
        ReachablePoints.push_back(RMPointPos);
      }
    }
    return ReachablePoints;
  }
  std::vector<Vector3> ReachablePointsFinder(const Robot & SimRobot, const int & LinkInfoIndex, SignedDistanceFieldInfo & SDFInfo, const Vector3 & COMVel)
  {
    double Radius = EndEffectorRadius[LinkInfoIndex];
    double PivotalLinkIndex = EndEffectorPivotalIndex[LinkInfoIndex];

    Vector3 RefPoint, ZeroPos(0.0, 0.0, 0.0);
    SimRobot.GetWorldPosition(ZeroPos, PivotalLinkIndex, RefPoint);
    int ReachablePointNo = 0;
    return ReachablePointsGene(RefPoint, Radius, SDFInfo, ReachablePointNo, COMVel);
  }
  std::vector<Vector3> ReachablePointsFinder(const Robot & SimRobot, const int & LinkInfoIndex, SignedDistanceFieldInfo & SDFInfo)
  {
    double Radius = EndEffectorRadius[LinkInfoIndex];
    double PivotalLinkIndex = EndEffectorPivotalIndex[LinkInfoIndex];

    Vector3 RefPoint, ZeroPos(0.0, 0.0, 0.0);
    SimRobot.GetWorldPosition(ZeroPos, PivotalLinkIndex, RefPoint);
    Vector3 COMVel(0.0, 0.0, 0.0);
    int ReachablePointNo = 0;
    return ReachablePointsGene(RefPoint, Radius, SDFInfo, ReachablePointNo, COMVel);
  }
  std::vector<Vector3> ReachablePointsGene(const Vector3 & RefPoint, const double & Radius, SignedDistanceFieldInfo & SDFInfo, int & ReachablePointNo, const Vector3 & COMVel)
  {
    // Here MinRadius is used for comparison while Radius is the maximum allowed radius of robot's end effector.
    // This function is used to get presumably active ReachablePointsGene() from all sampled points.
    const double DisTol = 0.01;        // 1cm as a signed distance tolerance.
    std::vector<Vector3> ReachablePoints;
    ReachablePoints.reserve(TotalPoint);
    ReachablePointNo = 0;
    int LayerIndex = 0;
    double LayerRadius = MinRadius;
    if(Radius>MaxRadius)
    {
      // Then LayerNumber can all be used.
      LayerIndex = LayerNumber-1;
    }
    else
    {
      while (LayerRadius<Radius)
      {
        LayerRadius+=LayerDiff;
        LayerIndex++;
      }
      LayerRadius-=LayerDiff;
    }

    for (int i = 0; i < LayerIndex; i++)
    {
      std::vector<RMPoint> RMLayer_i = RMLayers[i];
      for (int j = 0; j < RMLayer_i.size(); j++)
      {
        Vector3 RMPointPos = RMLayer_i[j].Position + RefPoint;
        double CurrentDist = SDFInfo.SignedDistance(RMPointPos);
        if(CurrentDist*CurrentDist<DisTol*DisTol)
        {
          ReachablePoints.push_back(RMPointPos);
          ReachablePointNo++;
        }
      }
    }
    return ReachablePoints;
  }
  std::vector<Vector3> ContactFreePointsFinder(const double & radius, const std::vector<Vector3> & ReachablePoints,const std::vector<std::pair<Vector3, double>> & ContactFreeInfo)
  {
    // This function can only be called after ReachablePointsFinder() to reduce the extra point further.
    std::vector<Vector3> ContactFreePoints;
    ContactFreePoints.reserve(ReachablePoints.size());
    int ContactFreeNo = 0;
    for (int i = 0; i < ReachablePoints.size(); i++)
    {
      Vector3 ReachablePoint = ReachablePoints[i];
      bool ContactFreeFlag = true;
      int ContactFreeInfoIndex = 0;
      while (ContactFreeInfoIndex<ContactFreeInfo.size())
      {
        Vector3 RefPoint = ContactFreeInfo[ContactFreeInfoIndex].first;
        double Radius = ContactFreeInfo[ContactFreeInfoIndex].second;
        Vector3 PosDiff = ReachablePoint - RefPoint;
        double PosDiffDis = sqrt(PosDiff.x * PosDiff.x + PosDiff.y * PosDiff.y + PosDiff.z * PosDiff.z);
        if(PosDiffDis<=(Radius + radius))
        {
          ContactFreeFlag = false;
          break;
        }
        ContactFreeInfoIndex++;
      }
      if(ContactFreeFlag)
      {
        ContactFreePoints.push_back(ReachablePoints[i]);
        ContactFreeNo++;
      }
    }
    return ContactFreePoints;
  };
  std::map<int, std::vector<RMPoint>> RMLayers;       // Each layer contains several data points.
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

#endif
