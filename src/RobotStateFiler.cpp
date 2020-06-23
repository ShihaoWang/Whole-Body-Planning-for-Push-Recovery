// This function is used to load in the robot's specification
#include <iostream>
#include <fstream>
#include <sstream>
#include <Modeling/Robot.h>
#include "CommonHeader.h"

void RobotConfigLoader(Robot & SimRobot, const string & user_path, const string & file_name){
  string str_line, str_keyword;
  str_keyword = "\t";
  int flag = 0;

  string config_file_path = user_path + file_name;
  ifstream ConfigInfofile (config_file_path);

  std::vector<double> RobotConfig, RobotVelocity;
  if (ConfigInfofile.is_open())
  {
    while (getline (ConfigInfofile,str_line) )
    {
      size_t start_pos = str_line.find(str_keyword);        // Here gives out the position of the \t
      if (start_pos != string::npos)
      {
        string DOF_str ="";
        for (size_t j = 0; j < start_pos; j++)
        {
          DOF_str+= str_line[j];
        }
        const int DOF = stoi(DOF_str);
        RobotConfig.reserve(DOF);
        RobotVelocity.reserve(DOF);
        str_line.erase(0,start_pos+1);

        std::istringstream ss(str_line);
        string Config_i;
        while (ss >> Config_i)
        {
          RobotConfig.push_back(std::stod(Config_i));
          RobotVelocity.push_back(0);
        }
      }
      else
      {
        std::cout<<"Wrong! .config file cannot be found!"<<endl;
      }
    }
    ConfigInfofile.close();
    flag = 1;
  }
  else cout << "Unable to open file";

  // Update the SimRobot's status
  SimRobot.UpdateConfig(Config(RobotConfig));     // Here both the SimRobot.q and robot frames have already been updated.
  SimRobot.dq = RobotVelocity;
  return;
}

void RobotConfigWriter(const std::vector<double> & Config, const string &user_path, const string &config_file_name){
  std::ofstream ConfigInfoFile;
  std::string config_file_path = user_path + config_file_name;
  ConfigInfoFile.open (config_file_path);
  ConfigInfoFile<<std::to_string(Config.size())<<"\t";
  for (int i = 0; i < Config.size(); i++)
  {
    ConfigInfoFile << std::to_string(Config[i])<<" ";
  }
  ConfigInfoFile.close();
}

std::vector<int> TorsoLinkReader(const string & TorsoLinkFilePath){
  ifstream TorsoLinkFile (TorsoLinkFilePath);
  std::vector<int> TorsoLinkVec;
  int LinkIndex = -1;
  if (TorsoLinkFile.is_open())
  {
    string str_line;
    while (getline (TorsoLinkFile, str_line) )
    {
      int link_index = stoi(str_line);
      TorsoLinkVec.push_back(link_index);

    }
    TorsoLinkFile.close();
  }
  else std::cerr << "Unable to open file " <<TorsoLinkFilePath<<" does not exist!\n";

  if (TorsoLinkVec.size() == 0)
  {
    std::cerr<<"Robot Contact Status Info failed to be loaded!"<<"\n";
  }
  return TorsoLinkVec;
}

void StateTrajAppender(const char *stateTrajFile_Name, const double & Time_t, const std::vector<double> & Configuration){
  std::ofstream StateTrajWriter;
  StateTrajWriter.open(stateTrajFile_Name, std::ios_base::app);
  StateTrajWriter<<std::to_string(Time_t)<<"\t";
  StateTrajWriter<<std::to_string(Configuration.size())<<"\t";
  for (int i = 0; i < Configuration.size()-1; i++){
    StateTrajWriter<<std::to_string(Configuration[i])<<" ";
  }
  StateTrajWriter<<std::to_string(Configuration[Configuration.size()-1])<<"\n";
  StateTrajWriter.close();
}

void PushInfoFileAppender(const double & SimTime, const double & Fx_t, const double & Fy_t, const double & Fz_t, const string & SpecificPath){
  std::ofstream PushInfoFileWriter;
  string PushInfoFileStr = SpecificPath + "PushInfoFile.txt";
  const char *PushInfoFileStr_Name = PushInfoFileStr.c_str();
  PushInfoFileWriter.open(PushInfoFileStr_Name, std::ios_base::app);
  PushInfoFileWriter<<std::to_string(SimTime)<<" "<< std::to_string(Fx_t)<<" "<<std::to_string(Fy_t)<<" "<<std::to_string(Fz_t)<<"\n";
  PushInfoFileWriter.close();
  return;
}

void StateLogger(WorldSimulation & Sim, FailureStateInfo & FailureStateObj, LinearPath & CtrlStateTraj, LinearPath & PlanStateTraj, LinearPath & FailureStateTraj, std::vector<double> & qDes, const SimPara & SimParaObj){
  const char *FailureStateTrajStr_Name  = SimParaObj.FailureStateTrajStr.c_str();
  const char *CtrlStateTrajStr_Name     = SimParaObj.CtrlStateTrajStr.c_str();
  const char *PlanStateTrajStr_Name     = SimParaObj.PlanStateTrajStr.c_str();

  CtrlStateTraj.Append(Sim.time,    Sim.world->robots[0]->q);
  StateTrajAppender(CtrlStateTrajStr_Name, Sim.time, Sim.world->robots[0]->q);

  if(qDes.size()==0) qDes = PlanStateTraj.milestones[PlanStateTraj.times.size()-1];

  StateTrajAppender(PlanStateTrajStr_Name, Sim.time, qDes);
  PlanStateTraj.Append(Sim.time,    Config(qDes));

  if(!FailureStateObj.FailureStateFlag){
    FailureStateTraj.Append(Sim.time,    Sim.world->robots[0]->q);
    StateTrajAppender(FailureStateTrajStr_Name, Sim.time, Sim.world->robots[0]->q);
  }
  return;
}
