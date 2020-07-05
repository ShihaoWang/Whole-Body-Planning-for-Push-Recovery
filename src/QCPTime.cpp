#include "CommonHeader.h"
#include "NonlinearOptimizerInfo.h"
#include "gurobi_c++.h"

double AccPhaseTimePathMethodQCP(   const std::vector<double> & CurConfig,      const std::vector<double> & NextConfig,
                                    const std::vector<double> & CurVelocity,    std::vector<double> & NextVelocity,
                                    const std::vector<double> & VelocityBound,  const std::vector<double> & AccelerationBound,
                                    const std::vector<int> SwingLinkChain){
  /*
      Variables to be optimized:
      Joint's velocity, acceleration and time
  */

  try {
    GRBEnv env = GRBEnv();
    GRBModel model = GRBModel(env);

    // Create variables
    std::vector<GRBVar> OptVariables;
    OptVariables.reserve(SwingLinkChain.size() * 2 + 2);
    // qdot
    for (int i = 0; i < SwingLinkChain.size(); i++){
      std::string x_name = "qdot_" + std::to_string(i);
      double qdot_max = VelocityBound[SwingLinkChain[i]];
      GRBVar qdot_i = model.addVar(-1.0 * qdot_max, qdot_max, 0.0, GRB_CONTINUOUS, x_name);
      OptVariables.push_back(qdot_i);
    }
    // qddot
    for (int i = 0; i < SwingLinkChain.size(); i++){
      std::string x_name = "qddot_" + std::to_string(i);
      double qddot_max = AccelerationBound[SwingLinkChain[i]];
      GRBVar qdot_i = model.addVar(-1.0 * qddot_max, qddot_max, 0.0, GRB_CONTINUOUS, x_name);
      OptVariables.push_back(qdot_i);
    }
    std::string delta_t_name = "delta_t";
    GRBVar delta_t = model.addVar(0.0, 100.0, 0.0, GRB_CONTINUOUS, delta_t_name);
    OptVariables.push_back(delta_t);

    std::string s_name = "s";
    GRBVar s = model.addVar(0.0, 10000.0, 0.0, GRB_CONTINUOUS, s_name);
    OptVariables.push_back(s);

    // Set objective: s = delta_t^2
    GRBQuadExpr obj = OptVariables.back();
    model.setObjective(obj);

    // Set Constraint
    model.addQConstr(s == delta_t * delta_t, "slack");
    int ConsInd = 0;
    std::string cons_name;
    for (int i = 0; i < SwingLinkChain.size(); i++) {
      double delta_q = NextConfig[SwingLinkChain[i]] - CurConfig[SwingLinkChain[i]];
      cons_name = "integration" + std::to_string(ConsInd);
      model.addQConstr(delta_q == OptVariables[i] * delta_t + 0.5 * OptVariables[i+SwingLinkChain.size()] * s, cons_name);
      ConsInd++;
    }
    for (int i = 0; i < SwingLinkChain.size(); i++) {
      double qdot_max = VelocityBound[SwingLinkChain[i]];
      cons_name = "velocity" + std::to_string(ConsInd);
      model.addQConstr(OptVariables[i] + OptVariables[i+SwingLinkChain.size()] * delta_t<=qdot_max, cons_name);
      ConsInd++;
      cons_name = "velocity" + std::to_string(ConsInd);
      model.addQConstr(OptVariables[i] + OptVariables[i+SwingLinkChain.size()] * delta_t>=-qdot_max, cons_name);
      ConsInd++;
    }

    model.optimize();
    cout << "Objective Value: " << model.get(GRB_DoubleAttr_ObjVal)<< endl;
    std::vector<double> qdot_soln(OptVariables.size());
    for (int i = 0; i < OptVariables.size(); i++){
      qdot_soln[i] = OptVariables[i].get(GRB_DoubleAttr_X);
      std::cout<<qdot_soln[i]<<std::endl;
    }
    return sqrt(qdot_soln.back());

  } catch(GRBException e)
  {
    cout << "Error code = " << e.getErrorCode() << endl;
    cout << e.getMessage() << endl;
  } catch(...)
  {
    cout << "Exception during optimization" << endl;
  }
  return 0.0;
}
