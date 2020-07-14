#include <iostream>
#include <fstream>
#include <sstream>
#include "CommonHeader.h"
#include <omp.h>
#include "NonlinearOptimizerInfo.h"

/*
    This function computes the robot's whole-body configuration based on the inverted pendulum rotation motion around axis.
*/
// Calculates rotation matrix to euler angles
// The result is the same as MATLAB except the order
// of the euler angles ( x and z are swapped ).
static Vector3 RotMat2EulerAngles(Matrix3 & RotMat){
  float sy = sqrt(RotMat(0,0) * RotMat(0,0) +  RotMat(1,0) * RotMat(1,0) );

  bool singular = sy < 1e-6; // If

  float x, y, z;
  if (!singular)
  {
    x = atan2(RotMat(2,1) , RotMat(2,2));
    y = atan2(-RotMat(2,0), sy);
    z = atan2(RotMat(1,0), RotMat(0,0));
  }
  else
  {
    x = atan2(-RotMat(1,2), RotMat(1,1));
    y = atan2(-RotMat(2,0), sy);
    z = 0;
  }
  return Vector3(x, y, z);
}

static Vector3 RigidBodyRotation(const Vector3 & RigidBodyPoint, const double & RotAngle, const Vector3 & RotAxis, const Vector3 & AxisOri){
  // This function is used to calcualte the new RigidBodyPoint after the rotation around edge_a->edge_b axis for Angle degree.
  AngleAxisRotation RotMatrix(RotAngle, RotAxis);
  Vector3 RigidBodyPointNew;
  RotMatrix.transformPoint(RigidBodyPoint - AxisOri, RigidBodyPointNew);
  RigidBodyPointNew+=AxisOri;
  return RigidBodyPointNew;
}

static void StepIntegrator(InvertedPendulumInfo & InvertedPendulumObj, const Vector3 & RotAxis, const double & TimeStep){
  // This function is used to integrate robot's PIP dynamics.
  double L = InvertedPendulumObj.L;
  double g = InvertedPendulumObj.g;

  // Integration with the assumption that robot's acceleration remains to be constant during TimeStep.
  double Thetaddot = g/L * sin(InvertedPendulumObj.Theta);
  double ThetaOffset = InvertedPendulumObj.Thetadot * TimeStep + 0.5 * Thetaddot * TimeStep * TimeStep;
  double ThetadotOffset = Thetaddot * TimeStep;

  double ThetaNew = InvertedPendulumObj.Theta + ThetaOffset;
  double ThetadotNew = InvertedPendulumObj.Thetadot + ThetadotOffset;

  Vector3 COMPosNew = RigidBodyRotation(InvertedPendulumObj.COMPos, -ThetaOffset, RotAxis, InvertedPendulumObj.edge_a);
  Vector3 COMPosOnEdge = InvertedPendulumObj.edge_a + RotAxis.dot(COMPosNew - InvertedPendulumObj.edge_a) * RotAxis;
  Vector3 COMVelDir;
  COMVelDir.setNormalized(cross(COMPosNew - COMPosOnEdge, RotAxis));
  Vector3 COMVelNew = L * ThetadotNew * COMVelDir;

  InvertedPendulumObj.Theta     =  ThetaNew;
  InvertedPendulumObj.Thetadot  =  ThetadotNew;
  // InvertedPendulumObj.COMVel    =  (COMPosNew - InvertedPendulumObj.COMPos)/TimeStep;
  InvertedPendulumObj.COMVel    =  COMVelNew;
  InvertedPendulumObj.COMPos    =  COMPosNew;
  return;
}

static std::vector<double> GlobalFrameConfigUpdate(Robot & SimRobot, const double & ThetaOffset, const Vector3 & RotAxis, const Vector3 & AxisOri){
  // This function is used to update robot's configuration for global frame..
  // First part is frame's Euclidean position.
  Vector3 FramePos, FramePos1;
  SimRobot.GetWorldPosition(Vector3(0.0, 0.0, 0.0), 0, FramePos);   // This gets robot's position of global frame.
  SimRobot.GetWorldPosition(Vector3(0.0, 0.0, 0.0), 1, FramePos);   // This gets robot's position of global frame
  SimRobot.GetWorldPosition(Vector3(0.0, 0.0, 0.0), 2, FramePos);   // This gets robot's position of global frame.

  Vector3 FrameOri, FrameXaxis, FrameYaxis, FrameZaxis;
  SimRobot.GetWorldPosition(Vector3(0.0, 0.0, 0.0), 5, FrameOri);       // This gets robot's position of global frame.
  SimRobot.GetWorldPosition(Vector3(0.0, 0.0, 1.0), 5, FrameXaxis);
  SimRobot.GetWorldPosition(Vector3(0.0, 1.0, 0.0), 5, FrameYaxis);
  SimRobot.GetWorldPosition(Vector3(1.0, 0.0, 0.0), 5, FrameZaxis);

  Vector3 FramePosNew = RigidBodyRotation(FramePos, -ThetaOffset, RotAxis, AxisOri);

  Vector3 x_axis, y_axis, z_axis;
  x_axis = FrameXaxis - FrameOri;
  y_axis = FrameYaxis - FrameOri;
  z_axis = FrameZaxis - FrameOri;
  x_axis.setNormalized(x_axis);
  y_axis.setNormalized(y_axis);
  z_axis.setNormalized(z_axis);

  Matrix3 RotMat(z_axis, y_axis, x_axis);
  AngleAxisRotation RotMatrix(-ThetaOffset, RotAxis);
  Matrix3 NewRotMat;
  RotMatrix.getMatrix(NewRotMat);

  NewRotMat.mul(NewRotMat, RotMat);
  Vector3 EulerAngle = RotMat2EulerAngles(NewRotMat);    // Reverse order to update frame's yaw, pitch and roll.
  if(EulerAngle.x>M_PI)
  {
    EulerAngle.x-=2.0 * M_PI;
  }
  if(EulerAngle.x<-M_PI)
  {
    EulerAngle.x+=2.0 * M_PI;
  }
  double FrameConfig[] = {FramePosNew.x, FramePosNew.y, FramePosNew.z, EulerAngle.z , EulerAngle.y, EulerAngle.x};
  std:;vector<double> FrameConfigVec(FrameConfig, FrameConfig + 6);
  return FrameConfigVec;
}

Config WholeBodyDynamicsIntegrator(Robot & SimRobot, InvertedPendulumInfo & InvertedPendulumObj, const double & TimeDuration){
  Vector3 RotAxis = InvertedPendulumObj.edge_b - InvertedPendulumObj.edge_a;
  RotAxis.setNormalized(RotAxis);
  const int IntergrationStep = 11;
  int IntergrationIndex = 0;
  double TimeStep = TimeDuration/(1.0 * IntergrationStep - 1.0);
  double ThetaInit = InvertedPendulumObj.Theta;
  for (int IntergrationIndex = 0; IntergrationIndex < IntergrationStep; IntergrationIndex++){
    StepIntegrator(InvertedPendulumObj, RotAxis, TimeStep);
  }

  double ThetaOffset = InvertedPendulumObj.Theta - ThetaInit;
  std::vector<double> FrameConfig = GlobalFrameConfigUpdate(SimRobot, ThetaOffset, RotAxis, InvertedPendulumObj.edge_a);   // This would be 6 global coordinates.
  std::vector<double> UpdatedConfig = SimRobot.q;
  for (int i = 0; i < 6; i++){
    UpdatedConfig[i] = FrameConfig[i];
  }
  return Config(UpdatedConfig);
}

std::vector<double> CurrentBaseDeltaCal(const Robot & _SimRobot, const InvertedPendulumInfo & _InvertedPendulumObj, const double & TimeDuration){
  Robot SimRobot = _SimRobot;
  InvertedPendulumInfo InvertedPendulumObj = _InvertedPendulumObj;
  Vector3 RotAxis = InvertedPendulumObj.edge_b - InvertedPendulumObj.edge_a;
  RotAxis.setNormalized(RotAxis);
  const int IntergrationStep = 11;
  int IntergrationIndex = 0;
  double TimeStep = TimeDuration/(1.0 * IntergrationStep - 1.0);
  double ThetaInit = InvertedPendulumObj.Theta;
  for (int IntergrationIndex = 0; IntergrationIndex < IntergrationStep; IntergrationIndex++){
    StepIntegrator(InvertedPendulumObj, RotAxis, TimeStep);
  }

  std::vector<double> CurrentBaseDelta(6);

  double ThetaOffset = InvertedPendulumObj.Theta - ThetaInit;
  std::vector<double> FrameConfig = GlobalFrameConfigUpdate(SimRobot, ThetaOffset, RotAxis, InvertedPendulumObj.edge_a);   // This would be 6 global coordinates.
  std::vector<double> UpdatedConfig = SimRobot.q;
  for (int i = 0; i < 6; i++)
    CurrentBaseDelta[i] = FrameConfig[i] - _SimRobot.q[i];
  return CurrentBaseDelta;
}
