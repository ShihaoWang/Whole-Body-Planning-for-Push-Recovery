// This function is used to load in the robot's specification
#include <iostream>
#include <fstream>
#include <sstream>
#include "CommonHeader.h"

std::vector<LinkInfo> ContactInfoLoader(const string & ContactLinkFile, int & _ContactPointNo)
{
  string str_line, str_keyword;
  str_keyword = "Link";
  ifstream linkinfofile (ContactLinkFile);
  std::vector<LinkInfo> LinkInfoVec;
  LinkInfo contact_link_i;
  int LinkIndex = -1;
  if (linkinfofile.is_open())
  {
    while (getline (linkinfofile,str_line) )
    {
      if (str_line.find(str_keyword) != string::npos)
      {
        str_line.erase(str_line.begin(), str_line.begin()+4);
        int link_index = stoi(str_line);
        LinkInfo contact_link_point_i(link_index);
        contact_link_i = contact_link_point_i;
        LinkInfoVec.push_back(contact_link_i);
        LinkIndex = LinkIndex + 1;
      }
      else
      {
        // Each row will be converted into a Vector3 and push into LinkInfo
        istringstream ss(str_line);     // Here the string row can be separated by the comma
        string link_coordinate_i;
        vector <double> link_coordinates;
        while (ss >> link_coordinate_i)
        {
          link_coordinate_i.erase(link_coordinate_i.end()-1, link_coordinate_i.end());
          link_coordinates.push_back(stod(link_coordinate_i));
        }
        Vector3 link_coordinates_Vect3(link_coordinates[0], link_coordinates[1], link_coordinates[2]);
        LinkInfoVec[LinkIndex].AddLocalConact(link_coordinates_Vect3);
      }
    }
    for (int i = 0; i < LinkInfoVec.size(); i++)
    {
      LinkInfoVec[i].AvgContactUpdate();
    }
    linkinfofile.close();
  }
  else cout << "Unable to open file";

  int ContactPointNo = 0;

  for (int i = 0; i < LinkInfoVec.size(); i++)
  {
    for (int j = 0; j < LinkInfoVec[i].LocalContacts.size(); j++)
    {
      ContactPointNo = ContactPointNo + 1;
    }
  }
  _ContactPointNo = ContactPointNo;

  if (LinkInfoVec.size() == 0)
  {
    std::cerr<<"Robot Contact Link Info failed to be loaded!"<<"\n";
  }
  return LinkInfoVec;
}

std::vector<ContactStatusInfo> ContactStatusInfoLoader(const string & ContactStatusFile)
{
  string str_line, str_keyword;
  str_keyword = "Link";

  ifstream contactstatusinfofile (ContactStatusFile);
  std::vector<ContactStatusInfo> ContactStatusInfoVec;
  ContactStatusInfo contact_link_i;
  int LinkIndex = -1;
  if (contactstatusinfofile.is_open())
  {
    while (getline (contactstatusinfofile,str_line) )
    {
      if (str_line.find(str_keyword) != string::npos)
      {
        str_line.erase (str_line.begin(), str_line.begin()+4);
        int link_index = stoi(str_line);
        ContactStatusInfo contact_link_point_i(link_index);
        contact_link_i = contact_link_point_i;
        ContactStatusInfoVec.push_back(contact_link_i);
        LinkIndex = LinkIndex + 1;
      }
      else
      {
        ContactStatusInfoVec[LinkIndex].AddLocalConactStatus(std::stoi(str_line));
      }
    }
    contactstatusinfofile.close();
  }
  else std::cerr << "\nUnable to open file";

  if (ContactStatusInfoVec.size() == 0)
  {
    std::cerr<<"\nRobot Contact Status Info failed to be loaded!"<<"\n";
  }
  return ContactStatusInfoVec;
}
