#include "CommonHeader.h"
#include "NonlinearOptimizerInfo.h"
#include "gurobi_c++.h"

static double QCons_tol = 1e-3;

double TrajectoryPhaseTime( Robot & SimRobotInner, const std::vector<double> & CurrentConfig, const std::vector<double> & CurrentVelocity,
                            const std::vector<double> & NextConfig, std::vector<double> & NextVelocity,
                            const std::vector<int> & SwingLinkChain, SimPara & SimParaObj, const int & sIndex){
  int SwingLinkInfoIndex = SimParaObj.getSwingLinkInfoIndex();
  /*
    The problem is formulated to be a linear programming problem
                          max_{qdot, delta_s} -dela_t
                    s.t   q = q_ref + qdot_ref * delta_t +
                                                        where d(q,s) = s' + k_{goal} * (GoalPos - CurrentPos) + k_{fit} * (s^* - CurrentPos)
                                                  |q_ref + qdot * delta_t|<=q bound
                                                  |qdot|<=qdot bound
                                                  |qdot - qdot_pre|<=qddot * delta_t
  */
  try {
    GRBEnv env = GRBEnv();
    GRBModel model = GRBModel(env);

    // Create variables: delta_t, qdot and qddot
    std::vector<GRBVar> OptVariables;
    OptVariables.reserve(2 * SwingLinkChain.size() + 2);
    int VarInd = 0;
    // qdot
    for (int i = 0; i < SwingLinkChain.size(); i++){
      std::string var_name = "qdot" + std::to_string(VarInd);
      double vel_bound = SimRobotInner.velMax[SwingLinkChain[i]];
      GRBVar qdot_var = model.addVar(-1.0 * vel_bound, vel_bound, 0.0, GRB_CONTINUOUS, var_name);
      OptVariables.push_back(qdot_var);
      VarInd += 1;
    }
    // qddot
    for (int i = 0; i < SwingLinkChain.size(); i++){
      std::string var_name = "qddot" + std::to_string(VarInd);
      double acc_bound = SimRobotInner.accMax[SwingLinkChain[i]];
      GRBVar qddot_var = model.addVar(-1.0 * acc_bound, acc_bound, 0.0, GRB_CONTINUOUS, var_name);
      OptVariables.push_back(qddot_var);
      VarInd += 1;
    }
    // delta_t
    std::string delta_t_var_name = "delta_t";
    GRBVar delta_t_var = model.addVar(0.0, 100.0, 0.0, GRB_CONTINUOUS, delta_t_var_name);
    OptVariables.push_back(delta_t_var);
    // slack_s
    std::string slack_s_var_name = "slack_s";   // s = delta_t^2
    GRBVar slack_s_var = model.addVar(0.0, 1000.0, 0.0, GRB_CONTINUOUS, slack_s_var_name);
    OptVariables.push_back(slack_s_var);

    // Set objective
    GRBLinExpr obj = delta_t_var;        // Minimize delta_t
    model.setObjective(obj);

    // Add constraint:
    // 1. configuration bound
    int ConsInd = 0;
    for (int i = 0; i < SwingLinkChain.size(); i++) {
      GRBQuadExpr ConfigEquality = CurrentConfig[SwingLinkChain[i]] + CurrentVelocity[SwingLinkChain[i]] * delta_t_var +
                                    0.5 * OptVariables[SwingLinkChain.size() + i] * slack_s_var - NextConfig[SwingLinkChain[i]];
      std::string cons_name = "c" + std::to_string(ConsInd);
      model.addQConstr(ConfigEquality>=-QCons_tol, cons_name);
      ConsInd+=1;
      cons_name = "c" + std::to_string(ConsInd);
      model.addQConstr(ConfigEquality<=QCons_tol, cons_name);
      ConsInd+=1;
    }
    // 2. acceleration bound
    for (int i = 0; i < SwingLinkChain.size(); i++) {
      GRBLinExpr VelocityChange = OptVariables[i] - CurrentVelocity[SwingLinkChain[i]];
      GRBQuadExpr AccInte = OptVariables[SwingLinkChain.size() + i] * delta_t_var;
      std::string cons_name = "c" + std::to_string(ConsInd);
      model.addQConstr(VelocityChange==AccInte, cons_name);
      ConsInd+=1;
    }
    // 3. slack variable
    std::string cons_name = "s_delta_t_sq" + std::to_string(ConsInd);
    model.addQConstr(slack_s_var==delta_t_var * delta_t_var, cons_name);
    ConsInd+=1;

    // Optimize model
    model.optimize();
    std::vector<double> OptValues(OptVariables.size());
    for (int i = 0; i < OptVariables.size(); i++){
      OptValues[i] = OptVariables[i].get(GRB_DoubleAttr_X);
      printf("OptVarible %d's Value: %f\n", i ,OptValues[i]);
    }

    double delta_t = OptValues[OptValues.size()-2];
    double slack_s = OptValues[OptValues.size()-1];
    std::vector<double> OptConfiguration = CurrentConfig;
    for (int i = 0; i < SwingLinkChain.size(); i++) {
      OptConfiguration[SwingLinkChain[i]] += OptValues[i] * delta_t + 0.5 * OptValues[i + SwingLinkChain.size()] *delta_t * delta_t;
    }

    OptConfiguration = YPRShifter(OptConfiguration);

    std::string ConfigPath = "/home/motion/Desktop/Whole-Body-Planning-for-Push-Recovery/build/";
    std::string OptConfigFile = "UpdatedConfig" + std::to_string(sIndex) + ".config";
    RobotConfigWriter(OptConfiguration, ConfigPath, OptConfigFile);

    std::cout << "Obj: " << model.get(GRB_DoubleAttr_ObjVal) << endl;
  } catch(GRBException e) {
    std::cout << "Error code = " << e.getErrorCode() << endl;
    std::cout << e.getMessage() << endl;
  } catch(...) {
    std::cout << "Exception during optimization" << endl;
  }
}
