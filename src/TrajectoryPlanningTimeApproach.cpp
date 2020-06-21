#include "CommonHeader.h"
#include "NonlinearOptimizerInfo.h"
#include "gurobi_c++.h"

static Matrix getActiveTaskSpaceJacobian(const Robot & SimRobot, std::vector<int> SwingLinkChain, const int & SwingLinkInfoIndex){
  // Should be a 3 by 6 + SwingLinkChain.size() dimension matrix
  Matrix TaskSpaceJacobian;
  std::sort(SwingLinkChain.begin(), SwingLinkChain.end());
  SimRobot.GetPositionJacobian( NonlinearOptimizerInfo::RobotLinkInfo[SwingLinkInfoIndex].AvgLocalContact,
                                NonlinearOptimizerInfo::RobotLinkInfo[SwingLinkInfoIndex].LinkIndex, TaskSpaceJacobian);
  // However, this matrix contains redundant information.
  Matrix ActiveTaskSpaceJacobian;
  ActiveTaskSpaceJacobian.resize(3, 6 + SwingLinkChain.size());
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 6; j++) {
      ActiveTaskSpaceJacobian(i, j) = TaskSpaceJacobian(i, j);
    }
  }
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < SwingLinkChain.size(); j++) {
      ActiveTaskSpaceJacobian(i, j) = TaskSpaceJacobian(i, SwingLinkChain[j]);
    }
  }
  return ActiveTaskSpaceJacobian;
}

Config AccPhaseTaskSpaceMethod( Robot & SimRobotInner, const std::vector<double> & CurrentConfig, const std::vector<double> & CurrentVelocity,
                                const InvertedPendulumInfo & InvertedPendulumInner, SelfLinkGeoInfo & SelfLinkGeoObj,
                                EndEffectorPathInfo & EndEffectorPathObj, const std::vector<int> & SwingLinkChain, SimPara & SimParaObj){
  int SwingLinkInfoIndex = SimParaObj.getSwingLinkInfoIndex();
  double delta_t = SimParaObj.PhaseTimeStep;
  Matrix Jac = getActiveTaskSpaceJacobian(SimRobotInner, SwingLinkChain, SwingLinkInfoIndex);
  Vector3 CurrentPos;
  SimRobotInner.GetWorldPosition( NonlinearOptimizerInfo::RobotLinkInfo[SwingLinkInfoIndex].AvgLocalContact,
                                  NonlinearOptimizerInfo::RobotLinkInfo[SwingLinkInfoIndex].LinkIndex,
                                  CurrentPos);
  /*
    The problem is formulated to be a linear programming problem
                          max_{qdot, delta_s} delta_s
                    s.t   Jac(q) * qdot * delta_t = delta_s * d(q,s),
                                                        where d(q,s) = s' + k_{goal} * (GoalPos - CurrentPos) + k_{fit} * (s^* - CurrentPos)
                                                  |q_ref + qdot * delta_t|<=q bound
                                                  |qdot|<=qdot bound
                                                  |qdot - qdot_pre|<=qddot * delta_t
  */
  double k_goal = 1.0;
  double k_fit  = 1.0;
  double sPos = EndEffectorPathObj.Pos2s(CurrentPos);
  Vector3 splinePos, splineVel;
  EndEffectorPathObj.PosNTang(sPos, splinePos, splineVel);
  Vector3 GoalPos = SimParaObj.getContactGoal();

  try {
    GRBEnv env = GRBEnv();
    GRBModel model = GRBModel(env);

    // Create variables: delta_s and qdot
    std::vector<GRBVar> OptVariables;
    OptVariables.reserve(1 + SwingLinkChain.size());

    // delta_s
    std::string delta_s_var_name = "delta_s";
    GRBVar delta_s_var = model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS, delta_s_var_name);
    OptVariables.push_back(delta_s_var);
    int VarInd = 1;
    // qdot
    for (int i = 0; i < SwingLinkChain.size(); i++){
      std::string var_name = "qdot" + std::to_string(VarInd);
      double vel_bound = SimRobotInner.velMax[SwingLinkChain[i]];
      GRBVar qdot_var = model.addVar(-1.0 * vel_bound, vel_bound, 0.0, GRB_CONTINUOUS, var_name);
      OptVariables.push_back(qdot_var);
      VarInd += 1;
    }

    // Set objective
    GRBLinExpr obj = -OptVariables[0];        // Maximize delta_s <=> minimize -delta_s
    model.setObjective(obj);

    // Add constraint:
    // 0. J(q) * qdot * delta_t = delta_s * d(q,s)
    Vector3 d_q_s = splineVel + k_goal * (GoalPos - CurrentPos) + k_fit * (splinePos - CurrentPos);
    int ConsInd = 0;
    for (int i = 0; i < 3; i++){
      GRBLinExpr Jac_qdot_delta_t = 0;
      for (int j = 0; j < 6; j++)
        Jac_qdot_delta_t+=Jac(i,j) * CurrentVelocity[j];
      for (int j = 0; j < SwingLinkChain.size(); j++)
        Jac_qdot_delta_t+=Jac(i,6 + j) * OptVariables[1+j];

      GRBLinExpr delta_s_d_q_s = OptVariables[0] * d_q_s[i];
      std::string cons_name = "c" + std::to_string(ConsInd);
      model.addConstr(Jac_qdot_delta_t==delta_s_d_q_s, cons_name);
      ConsInd+=1;
    }

    // 1. Configuration bound
    for (int i = 0; i < SwingLinkChain.size(); i++) {
      GRBLinExpr UpdatedConfig = CurrentConfig[SwingLinkChain[i]] + OptVariables[1+i] * delta_t;
      std::string cons_name = "c" + std::to_string(ConsInd);
      double config_lb = SimRobotInner.qMin(SwingLinkChain[i]);
      double config_ub = SimRobotInner.qMax(SwingLinkChain[i]);
      model.addConstr(UpdatedConfig>=config_lb, cons_name);
      ConsInd+=1;
      cons_name = "c" + std::to_string(ConsInd);
      model.addConstr(UpdatedConfig<=config_ub, cons_name);
      ConsInd+=1;
    }

    // 2. Acceleration bound
    for (int i = 0; i < SwingLinkChain.size(); i++) {
      GRBLinExpr EstAcc = (OptVariables[1+i] - CurrentVelocity[SwingLinkChain[i]])/delta_t;
      std::string cons_name = "c" + std::to_string(ConsInd);
      double acc_bnd = SimRobotInner.accMax(SwingLinkChain[i]);
      model.addConstr(EstAcc>=-acc_bnd, cons_name);
      ConsInd+=1;
      cons_name = "c" + std::to_string(ConsInd);
      model.addConstr(EstAcc<=acc_bnd, cons_name);
      ConsInd+=1;
    }

    // Optimize model

    model.optimize();
    std::vector<double> OptValues(OptVariables.size());
    for (int i = 0; i < OptVariables.size(); i++){
      OptValues[i] = OptVariables[i].get(GRB_DoubleAttr_X);
      printf("OptVarible %d's Value: %f\n", i ,OptValues[i]);
    }

    cout << "Obj: " << model.get(GRB_DoubleAttr_ObjVal) << endl;

  } catch(GRBException e) {
    cout << "Error code = " << e.getErrorCode() << endl;
    cout << e.getMessage() << endl;
  } catch(...) {
    cout << "Exception during optimization" << endl;
  }
}
