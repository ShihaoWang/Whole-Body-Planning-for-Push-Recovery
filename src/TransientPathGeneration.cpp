#include "CommonHeader.h"
#include "NonlinearOptimizerInfo.h"

using namespace SplineLib;

static void SplinePiece1DObjGene(const double & sInit, const double & sGoal, const double & PosInit, const double & VelInit, const double & PosGoal, const double & VelGoal, double & a, double & b, double & c, double & d)
{
  // This function is used to generate the piecewiese cubic spline.
  double M00 = 3.0 * sInit * sInit;
  double M01 = 2.0 * sInit;
  double M02 = 1.0;

  double M10 = 3.0 * sGoal * sGoal;
  double M11 = 2.0 * sGoal;
  double M12 = 1.0;

  double M20 = (sInit * sInit * sInit - sGoal * sGoal * sGoal);
  double M21 = sInit * sInit - sGoal * sGoal;
  double M22 = sInit - sGoal;

  Vector3 Col1(M00, M10, M20);
  Vector3 Col2(M01, M11, M21);
  Vector3 Col3(M02, M12, M22);

  Matrix3 CubicMat(Col1, Col2, Col3);
  Vector3 CubicVec(VelInit, VelGoal, PosInit-PosGoal);

  Matrix3 CubicMatInv;
  CubicMat.getInverse(CubicMatInv);

  Vector3 abc = CubicMatInv * CubicVec;

  a = abc.x;
  b = abc.y;
  c = abc.z;
  d = PosInit - a * sInit * sInit * sInit - b * sInit * sInit - c * sInit;
  return;
}

static SplineInfo SplinePiece3DObjGene(const double & sInit, const double & sGoal, const Vector3 & PosInit, const Vector3 & NormalInit, const Vector3 & PosGoal, const Vector3 & NormalGoal)
{
  double ax, bx, cx, dx;
  SplinePiece1DObjGene(sInit, sGoal, PosInit.x, NormalInit.x, PosGoal.x, NormalGoal.x, ax, bx, cx, dx);

  double ay, by, cy, dy;
  SplinePiece1DObjGene(sInit, sGoal, PosInit.y, NormalInit.y, PosGoal.y, NormalGoal.y, ay, by, cy, dy);

  double az, bz, cz, dz;
  SplinePiece1DObjGene(sInit, sGoal, PosInit.z, NormalInit.z, PosGoal.z, NormalGoal.z, az, bz, cz, dz);

  Vector3 a(ax, ay, az);
  Vector3 b(bx, by, bz);
  Vector3 c(cx, cy, cz);
  Vector3 d(dx, dy, dz);

  SplineInfo SplineObj(sInit, sGoal, a, b, c, d);
  return SplineObj;
}

static std::vector<Vector3> BasePointsGene(const Vector3 & PosInit, const Vector3 & NormalInit, const Vector3 & PosGoal, const Vector3 & NormalGoal)
{
  // This function is used to generate the cubic spline for given robot's end effector path.
  const double scale = 0.25;
  SplineInfo BaseSpline = SplinePiece3DObjGene(0.0, 1.0, PosInit, scale * NormalInit, PosGoal, -scale * NormalGoal);
  const int segmentNo = 5;
  double sUnit = 1.0/(1.0 * segmentNo);
  std::vector<Vector3> BasePoints(segmentNo+1);
  for (int i = 0; i < segmentNo + 1; i++)
  {
    double s = 1.0 * i * sUnit;
    Vector3 BasePoint = BaseSpline.SplinePosVector(s);
    BasePoints[i] = BasePoint;
  }
  return BasePoints;
}

static double SelfCollisionDist(SelfLinkGeoInfo & SelfLinkGeoObj, const int & SwingLinkInfoIndex, const std::vector<Vector3> & PointVec){
  double Distol = 0.1;      // 10cm
  std::vector<double> DistVec;
  for (const Vector3 & Point: PointVec) {
    double  PointDist; Vector3 PointGrad;
    SelfLinkGeoObj.SelfCollisionDistNGrad(SwingLinkInfoIndex, Point, PointDist, PointGrad);
    DistVec.push_back(PointDist);
  }
  double SplineMin = *std::min_element(DistVec.begin(), DistVec.end());
  if(SplineMin>Distol) return SplineMin;
  else return min(Distol, min(DistVec.front(), DistVec.back()));
}

static std::vector<cSpline3> cSplineGene(const std::vector<Vector3> & Points, const int & SwingLinkInfoIndex, const double & SelfTol, const SelfLinkGeoInfo & SelfLinkGeoObj, ReachabilityMap & RMObject, int & PointIndex, Vector3 & ShiftPoint, bool & FeasibleFlag)
{
  // Algorithm stops when tolerance cannot be satisfied!
  Vec3f SplinePoints[Points.size()];
  for (int i = 0; i < Points.size(); i++)
  {
    Vec3f PointVec(Points[i].x, Points[i].y, Points[i].z);
    SplinePoints[i] = PointVec;
  }
  const int numPoints = sizeof(SplinePoints) / sizeof(SplinePoints[0]);
  cSpline3 splines[numPoints + 1];
  int numSplines = SplinesFromPoints(numPoints, SplinePoints, numPoints + 1, splines);
  std::vector<cSpline3> SplineObj;
  SplineObj.reserve(numSplines);

  std::vector<Vector3> ShiftPointVec;
  std::vector<double> ShiftPointDisVec;
  std::vector<int> SplineIndexVec;
  const int GridNo = 10;
  float sUnit = 1.0/(1.0 * GridNo);
  for (int i = 0; i < numSplines; i++)
  {
    for (int j = 0; j < GridNo-1; j++)
    {
      float s = 1.0 * j * sUnit;
      Vec3f ps = Position (splines[i], s);
      Vector3 SplinePoint(ps.x, ps.y, ps.z);
      double SplinePointDis = NonlinearOptimizerInfo::SDFInfo.SignedDistance(SplinePoint);
      if(SplinePointDis<0)
      {
        SplineIndexVec.push_back(i);
        ShiftPointVec.push_back(SplinePoint);
        ShiftPointDisVec.push_back(SplinePointDis);
      }
    }
  }

  switch (ShiftPointVec.size())
  {
    case 0:
    {
      FeasibleFlag = true;
      for (int i = 0; i < numSplines; i++)
      {
        SplineObj.push_back(splines[i]);
      }
    }
    break;
    default:
    {
      // Choose the point with the largest penetration.
      FeasibleFlag = false;
      int ShiftPointIndex = std::distance(ShiftPointDisVec.begin(), std::min_element(ShiftPointDisVec.begin(), ShiftPointDisVec.end()));
      PointIndex = SplineIndexVec[ShiftPointIndex];
      ShiftPoint = ShiftPointVec[ShiftPointIndex];
    }
    break;
  }
  return SplineObj;
}

static std::vector<Vector3> SpatialPointShifter(const std::vector<Vector3> & Points, const int & SwingLinkInfoIndex, const double & SelfTol, SelfLinkGeoInfo & SelfLinkGeoObj, ReachabilityMap & RMObject, bool & ShiftFeasFlag){
  // This function is used to shift Points according to SelfTol
  // Here we only allow this shift to be conducted N times.
  const int TotalShiftTime = 10;
  std::vector<Vector3> NewPoints;
  NewPoints.push_back(Points[0]);

  int PointCount = 1;
  while (PointCount<Points.size()-1)
  {
    bool SelfShitFlag, EnviShitFlag;
    Vector3 CurPt = Points[PointCount];
    // Self-Collision Distance
    double  CurSelfPtDist;
    Vector3 CurSelfPtGrad;
    SelfLinkGeoObj.SelfCollisionDistNGrad(SwingLinkInfoIndex, CurPt, CurSelfPtDist, CurSelfPtGrad);
    double CurEnvPtDist = NonlinearOptimizerInfo::SDFInfo.SignedDistance(CurPt);
    if((CurSelfPtDist<SelfTol)||(CurEnvPtDist<0))
    {
      // Then shift action is needed!
      Vector3 NewPt = CurPt;
      int ShiftTime = 0;
      double ShiftDistUnit = 0.025;     // 2.5cm for example
      while (ShiftTime<TotalShiftTime)
      {
        Vector3 SelfDirection(0.0, 0.0, 0.0);
        SelfShitFlag = false;
        if(CurSelfPtDist<SelfTol)
        {
          SelfDirection = CurSelfPtGrad;
          SelfShitFlag = true;
        }
        Vector3 EnviDirection(0.0, 0.0, 0.0);
        EnviShitFlag = false;
        if(CurEnvPtDist<0)
        {
          EnviDirection = NonlinearOptimizerInfo::SDFInfo.SignedDistanceNormal(NewPt);
          EnviShitFlag = true;
        }
        if((!SelfShitFlag)&&(!EnviShitFlag))
        {
          break;
        }
        Vector3 ShiftDirection = SelfDirection + EnviDirection;
        ShiftDirection.setNormalized(ShiftDirection);
        NewPt += ShiftDistUnit * ShiftDirection;
        SelfLinkGeoObj.SelfCollisionDistNGrad(SwingLinkInfoIndex, NewPt, CurSelfPtDist, CurSelfPtGrad);
        CurEnvPtDist = NonlinearOptimizerInfo::SDFInfo.SignedDistance(CurPt);
        ShiftTime++;
      }
      if((!SelfShitFlag)&&(!EnviShitFlag))
      {
        NewPoints.push_back(NewPt);
      }
      else
      {
        ShiftFeasFlag = false;
        return NewPoints;
      }
    }
    else
    {
      NewPoints.push_back(CurPt);
    }
    PointCount++;
  }
  NewPoints.push_back(Points[Points.size()-1]);
  ShiftFeasFlag = true;
  return NewPoints;
}

static Vector3 SinglePointShifter(const Vector3 & Point, const int & SwingLinkInfoIndex, const double & SelfTol, SelfLinkGeoInfo & SelfLinkGeoObj, ReachabilityMap & RMObject, bool & ShiftFeasFlag)
{
  double ShiftDistUnit = 0.025;     // 2.5cm for example
  Vector3 NewPoint = Point;
  const int TotalShiftTime = 10;
  int CurShiftTime = 0;

  bool SelfShitFlag, EnviShitFlag;
  double  CurSelfPtDist;
  Vector3 CurSelfPtGrad;
  SelfLinkGeoObj.SelfCollisionDistNGrad(SwingLinkInfoIndex, Point, CurSelfPtDist, CurSelfPtGrad);
  double CurEnvPtDist = NonlinearOptimizerInfo::SDFInfo.SignedDistance(Point);
  ShiftFeasFlag = false;
  while(CurShiftTime<TotalShiftTime)
  {
    Vector3 SelfDirection(0.0, 0.0, 0.0);
    SelfShitFlag = false;
    if(CurSelfPtDist<SelfTol)
    {
      SelfDirection = CurSelfPtGrad;
      SelfShitFlag = true;
    }
    Vector3 EnviDirection(0.0, 0.0, 0.0);
    EnviShitFlag = false;
    if(CurEnvPtDist<0)
    {
      EnviDirection = NonlinearOptimizerInfo::SDFInfo.SignedDistanceNormal(NewPoint);
      EnviShitFlag = true;
    }
    if((!SelfShitFlag)&&(!EnviShitFlag))
    {
      break;
    }
    Vector3 ShiftDirection = SelfDirection + EnviDirection;
    ShiftDirection.setNormalized(ShiftDirection);
    NewPoint += ShiftDistUnit * ShiftDirection;
    SelfLinkGeoObj.SelfCollisionDistNGrad(SwingLinkInfoIndex, NewPoint, CurSelfPtDist, CurSelfPtGrad);
    CurEnvPtDist = NonlinearOptimizerInfo::SDFInfo.SignedDistance(NewPoint);
    CurShiftTime++;
  }
  if((!SelfShitFlag)&&(!EnviShitFlag))
  {
    ShiftFeasFlag = true;
  }
  return NewPoint;
}

static std::vector<cSpline3> SplineObjGene(SelfLinkGeoInfo & SelfLinkGeoObj, ReachabilityMap & RMObject, const int & SwingLinkInfoIndex, SimPara & SimParaObj){
  // This function is used to generate a collision-free path!
  Vector3 PosInit     = SimParaObj.ContactInit;
  Vector3 NormalInit  = SimParaObj.DirectionInit;
  Vector3 PosGoal     = SimParaObj.ContactGoal;
  Vector3 NormalGoal  = SimParaObj.DirectionGoal;

  std::vector<cSpline3> SplineObj;
  std::vector<Vector3> Points = BasePointsGene(PosInit, NormalInit, PosGoal, NormalGoal);
  Vector3Writer(Points, "InitialTransitionPoints");
  double SelfTol = SelfCollisionDist(SelfLinkGeoObj, SwingLinkInfoIndex, Points);
  bool InitShiftFeasFlag;     // For the shift of initial pts.
  Points = SpatialPointShifter(Points, SwingLinkInfoIndex, SelfTol, SelfLinkGeoObj, RMObject, InitShiftFeasFlag);
  Vector3Writer(Points, "ShiftedTransitionPoints");
  bool FeasiFlag = false;
  SimParaObj.setTransPathFeasiFlag(FeasiFlag);
  if(!InitShiftFeasFlag) return SplineObj;
  // Then the task is to generate a path which is collision-free.
  const int TotalIter = 10;
  int CurrentIter = 0;
  while((FeasiFlag == false)&&(CurrentIter<=TotalIter)){
    Vector3 ShiftPoint;
    int ShiftPointIndex;
    SplineObj = cSplineGene(Points, SwingLinkInfoIndex, SelfTol, SelfLinkGeoObj, RMObject, ShiftPointIndex, ShiftPoint, FeasiFlag);
    if(!FeasiFlag){
      bool ShiftFlag;
      Vector3 NewShiftPoint = SinglePointShifter(ShiftPoint, SwingLinkInfoIndex, SelfTol, SelfLinkGeoObj, RMObject, ShiftFlag);
      if(!ShiftFlag)  return SplineObj;
      else {
        // Insert the NewShiftPoint back into Points vector
        std::vector<Vector3> NewPoints(Points.begin(), Points.begin() + ShiftPointIndex);
        NewPoints.push_back(NewShiftPoint);
        NewPoints.insert(NewPoints.end(), Points.begin() + ShiftPointIndex + 1, Points.end());
        Points = NewPoints;
      }
    }
    CurrentIter++;
  }
  SimParaObj.setTransPathFeasiFlag(FeasiFlag);
  return SplineObj;
}

static Vector3 InitDirectionGene(const Vector3 & PosInit, const Vector3 & PosGoal, const Vector3 & DirGoal){
  // This function calcualtes the initial spline direction such that the resulting path is a parabola.
  Vector3 Goal2Init = PosInit - PosGoal;
  Goal2Init.getNormalized(Goal2Init);
  double Proj = Goal2Init.dot(DirGoal);
  Vector3 Vert = DirGoal - Proj * Goal2Init;
  Vector3 DirInit = Vert - Proj * Goal2Init;
  return DirInit;
}

std::vector<cSpline3> TransientPathGene(const Robot & SimRobot, ReachabilityMap & RMObject, SelfLinkGeoInfo & SelfLinkGeoObj, SimPara & SimParaObj){
  // This function generates the transition path for robot's end effector.
  // The path direction is chosen such that initial path is a parabola.
  Vector3 DirGoal = NonlinearOptimizerInfo::SDFInfo.SignedDistanceNormal(SimParaObj.ContactGoal);
  Vector3 DirInit = InitDirectionGene(SimParaObj.ContactInit, SimParaObj.ContactGoal, DirGoal);
  DirInit.setNormalized(DirInit);
  SimParaObj.setDirectionInit(DirInit);
  SimParaObj.setDirectionGoal(DirGoal);
  int SwingLinkInfoIndex = SimParaObj.getSwingLinkInfoIndex();
  std::vector<cSpline3> SplineObj = SplineObjGene(SelfLinkGeoObj, RMObject, SwingLinkInfoIndex, SimParaObj);

  if(SimParaObj.getTransPathFeasiFlag()){
    const int SplineNumber = SplineObj.size();
    const int SplineGrid = 5;
    std::vector<Vector3> TransitionPoints(SplineNumber * SplineGrid + 1);

    double sUnit = 1.0/(1.0 * SplineGrid);
    int TransitionIndex = 0;
    for (int i = 0; i < SplineNumber; i++)
    {
      for (int j = 0; j < SplineGrid; j++)
      {
        double s = 1.0 * j * sUnit;
        Vec3f ps = Position (SplineObj[i], s);
        Vector3 SplinePoint(ps.x, ps.y, ps.z);
        TransitionPoints[TransitionIndex] = SplinePoint;
        TransitionIndex++;
      }
    }
    Vec3f ps = Position (SplineObj[SplineNumber-1], 1.0);
    Vector3 SplinePoint(ps.x, ps.y, ps.z);
    TransitionPoints[TransitionIndex] = SplinePoint;
    SimParaObj.DataRecorderObj.TransitionPoints = TransitionPoints;
    Vector3Writer(TransitionPoints, "TransitionPoints");
  }
  return SplineObj;
}
