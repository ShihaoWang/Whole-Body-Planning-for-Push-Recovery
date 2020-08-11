import sys, os, time
from klampt import *
from klampt import vis
from klampt.vis.glrobotprogram import GLSimulationPlugin
from klampt.model.trajectory import Trajectory, RobotTrajectory
from scipy.interpolate import interp1d
import ipdb
import copy
from scipy.spatial import ConvexHull
import draw_hull
from OpenGL.GL import *
import math
import numpy as np

CurCase = "flat_1Contact"
ExpNo = 0

class MyGLPlugin(vis.GLPluginInterface):
    def __init__(self, world):
        vis.GLPluginInterface.__init__(self)
        self.world = world
        self.quit = False
        self.starp = False

    def mousefunc(self, button, state, x, y):
        print("mouse",button,state,x,y)
        if button==2:
            if state==0:
                print("Click list...",[o.getName() for o in self.click_world(x,y)])
            return True
        return False

    def motionfunc(self, x, y, dx, dy):
        return False

    def keyboardfunc(self, c, x, y):
        print("Pressed", c)
        return True

    def click_world(self, x, y):
        """Helper: returns a list of world objects sorted in order of
        increasing distance."""
        #get the viewport ray
        (s, d) = self.click_ray(x, y)

        #run the collision tests
        collided = []
        for g in self.collider.geomList:
            (hit, pt) = g[1].rayCast(s, d)
            if hit:
                dist = vectorops.dot(vectorops.sub(pt, s), d)
                collided.append((dist,g[0]))
        return [g[1] for g in sorted(collided)]

def my_draw_hull(h):
    glEnable(GL_LIGHTING)
    glEnable(GL_BLEND)
    glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA)
    glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE,[1.0,0.25,0.5,0.5])
    draw_hull.draw_hull(h)

def StringList2NumberList(str_list, string_name):
    # This function is used to convert the string list to a certain type of list
    if string_name =="float":
        res_list = [float(i) for i in str_list]
    else:
        res_list = [int(i) for i in str_list]
    return res_list

def ConfigurationLoaderfn(Config_Name):
    # This function is only used to load in the initial configuraiton
    # The initial file will be in the .config format
    with open(Config_Name,'r') as robot_angle_file:
        robotstate_angle_i = robot_angle_file.readlines()
    config_temp = [x.replace('\t',' ') for x in robotstate_angle_i]
    config_temp = [x.replace('\n','') for x in config_temp]
    config_temp = [float(i) for i in config_temp[0].split()]

    DOF = int(config_temp[0])
    # Config_Init = np.array(config_temp[1:])
    Config_Init = config_temp[1:]
    return DOF, Config_Init

def ReadTxtfn(file_name):
    Empty_List = []
    with open(file_name) as Txt_File:
        Txt_File_Str = Txt_File.read().splitlines()
        for Txt_File_Str_i in Txt_File_Str:
            Txt_File_Str_i.replace("\'","")
            Empty_List.append(float(Txt_File_Str_i))
    return Empty_List

def StateLoaderfn(*args):
    if len(args) == 2:
        # In this case, the robot is only given the configuration file
        config_file_path = args[1] + args[0]
        DOF, Config_Init = ConfigurationLoaderfn(config_file_path)
        # Then the Velocty_Init is set to be a zero value list
        Velocity_Init = []
        for i in range(0,DOF):
            Velocity_Init.append(0)
    else:
        if len(args) == 3:
            config_file_path = args[2] + args[0]
            velocity_file_path = args[2] + args[1]
            Config_Init = ReadTxtfn(config_file_path)
            Velocity_Init = ReadTxtfn(velocity_file_path)
            DOF = len(Config_Init)
        else:
            raise RuntimeError("Input name should be either one config file or two txt files!")
    return DOF, Config_Init, Velocity_Init

def ContactLinkReader(File_Path_Name):
    # This function is used to read-in the formation of the certain link and its associated contact points
    # The output of this function is a dictionary of several keys with multiple list values
    # File_Name = "./User_File/Contact_Link.txt"
    ContactLinkDictionary = dict()
    # The format of this function should be an integet with a list of contact points
    Link_Number_i = -1
    with open(File_Path_Name) as Txt_File:
        Txt_File_Str = Txt_File.read().splitlines()
        Dictionary_Value_Add_Flag = 0
        for Txt_File_Str_i in Txt_File_Str:
            if "Link" in Txt_File_Str_i:
                # This indicates a contact link number
                # Then the next few points would be the local coordinates of the contact extremities
                Txt_File_Str_i = Txt_File_Str_i.translate(None, 'Link')     # This step is to get the link number out
                Link_Number_i = int(Txt_File_Str_i)
                ContactLinkDictionary[Link_Number_i] = []
                continue
            if Link_Number_i>=0:
                # Here the Txt_File_Str_i is a string of 'a, b, c' format
                Txt_File_Str_i = Txt_File_Str_i.split(",")
                del Txt_File_Str_i[-1]
                Txt_File_Flt_i = StringList2NumberList(Txt_File_Str_i,"float")
                ContactLinkDictionary[Link_Number_i].append(Txt_File_Flt_i)
    return ContactLinkDictionary

def ContactStatusReader(File_Name, Path_Name):
    # This function is used to read-in the formation of the certain link and its associated contact points
    # The output of this function is a dictionary of several keys with multiple list values
    # File_Name = "./User_File/Contact_Link.txt"
    Contact_Status_Dictionary = dict()
    # The format of this function should be an integet with a list of contact points
    Link_Number_i = -1
    File_Path_Name = Path_Name + File_Name
    with open(File_Path_Name) as Txt_File:
        Txt_File_Str = Txt_File.read().splitlines()
        Dictionary_Value_Add_Flag = 0
        for Txt_File_Str_i in Txt_File_Str:
            if "Link" in Txt_File_Str_i:
                # This indicates a contact link number
                # Then the next few points would be the local coordinates of the contact extremities
                Txt_File_Str_i = Txt_File_Str_i.translate(None, 'Link')     # This step is to get the link number out
                Link_Number_i = int(Txt_File_Str_i)
                Contact_Status_Dictionary[Link_Number_i] = []
                continue
            if Link_Number_i>=0:
                Contact_Status_Dictionary[Link_Number_i].append(int(Txt_File_Str_i))
    return Contact_Status_Dictionary

def ConvexEdgeReader(File_Name, Path_Name):
    # This function is used to load in the Convex Edge for visualization
    # The format of this function should be an inteer with a list of contact points
    File_Path_Name = Path_Name + File_Name
    Convex_Edge_List = []
    with open(File_Path_Name) as Txt_File:
        Txt_File_Str = Txt_File.read().splitlines()
        Dictionary_Value_Add_Flag = 0
        for Txt_File_Str_i in Txt_File_Str:
            Txt_File_Str_i = Txt_File_Str_i.slit(" ")
            Txt_File_Flt_i = StringList2NumberList(Txt_File_Str_i,"float")
            Convex_Edge_List.append(Txt_File_Flt_i)
    return Convex_Edge_List

def PIPTrajReader(file_path):

    PIPList = []
    EdgeAList = []
    EdgeBList = []
    EdgeCOMList = []
    EdgexList = []
    EdgeyList = []
    EdgezList = []
    EdgeVertexList = []

    EdgeA_path = file_path + "/EdgeATraj.txt"
    with open(EdgeA_path,'r') as EdgeA_path_file:
        EdgeA_tot = EdgeA_path_file.readlines()
        for i in range(0, len(EdgeA_tot)):
            EdgeAList_i = []
            EdgeA_tot_i = EdgeA_tot[i].split(" ")
            EdgeAString = StringList2NumberList(EdgeA_tot_i[0:-1], "float")
            for j in range(0, len(EdgeAString)/3):
                EdgeAList_i.append(EdgeAString[3*j:3*j+3])
            EdgeAList.append(EdgeAList_i)

    EdgeB_path = file_path + "/EdgeBTraj.txt"
    with open(EdgeB_path,'r') as EdgeB_path_file:
        EdgeB_tot = EdgeB_path_file.readlines()
        for i in range(0, len(EdgeB_tot)):
            EdgeBList_i = []
            EdgeB_tot_i = EdgeB_tot[i].split(" ")
            EdgeBString = StringList2NumberList(EdgeB_tot_i[0:-1], "float")
            for j in range(0, len(EdgeBString)/3):
                EdgeBList_i.append(EdgeBString[3*j:3*j+3])
            EdgeBList.append(EdgeBList_i)

    EdgeCOM_path = file_path + "/EdgeCOMTraj.txt"
    with open(EdgeCOM_path,'r') as EdgeCOM_path_file:
        EdgeCOM_tot = EdgeCOM_path_file.readlines()
        for i in range(0, len(EdgeCOM_tot)):
            EdgeCOMList_i = []
            EdgeCOM_tot_i = EdgeCOM_tot[i].split(" ")
            EdgeCOMString = StringList2NumberList(EdgeCOM_tot_i[0:-1], "float")
            for j in range(0, len(EdgeCOMString)/3):
                EdgeCOMList_i.append(EdgeCOMString[3*j:3*j+3])
            EdgeCOMList.append(EdgeCOMList_i)

    Edgex_path = file_path + "/EdgexTraj.txt"
    with open(Edgex_path,'r') as Edgex_path_file:
        Edgex_tot = Edgex_path_file.readlines()
        for i in range(0, len(Edgex_tot)):
            EdgexList_i = []
            Edgex_tot_i = Edgex_tot[i].split(" ")
            EdgexString = StringList2NumberList(Edgex_tot_i[0:-1], "float")
            for j in range(0, len(EdgexString)/3):
                EdgexList_i.append(EdgexString[3*j:3*j+3])
            EdgexList.append(EdgexList_i)

    Edgey_path = file_path + "/EdgeyTraj.txt"
    with open(Edgey_path,'r') as Edgey_path_file:
        Edgey_tot = Edgey_path_file.readlines()
        for i in range(0, len(Edgey_tot)):
            EdgeyList_i = []
            Edgey_tot_i = Edgey_tot[i].split(" ")
            EdgeyString = StringList2NumberList(Edgey_tot_i[0:-1], "float")
            for j in range(0, len(EdgeyString)/3):
                EdgeyList_i.append(EdgeyString[3*j:3*j+3])
            EdgeyList.append(EdgeyList_i)

    Edgez_path = file_path + "/EdgezTraj.txt"
    with open(Edgez_path,'r') as Edgez_path_file:
        Edgez_tot = Edgez_path_file.readlines()
        for i in range(0, len(Edgez_tot)):
            EdgezList_i = []
            Edgez_tot_i = Edgez_tot[i].split(" ")
            EdgezString = StringList2NumberList(Edgez_tot_i[0:-1], "float")
            for j in range(0, len(EdgezString)/3):
                EdgezList_i.append(EdgezString[3*j:3*j+3])
            EdgezList.append(EdgezList_i)

    EdgeVertex_path = file_path + "/EdgeVertexTraj.txt"
    with open(EdgeVertex_path,'r') as EdgeVertex_path_file:
        EdgeVertex_tot = EdgeVertex_path_file.readlines()
        for i in range(0, len(EdgeVertex_tot)):
            EdgeVertexList_i = []
            EdgeVertex_tot_i = EdgeVertex_tot[i].split(" ")
            EdgezVertextring = StringList2NumberList(EdgeVertex_tot_i[0:-1], "float")
            for j in range(0, len(EdgezVertextring)/3):
                EdgeVertexList_i.append(EdgezVertextring[3*j:3*j+3])
            EdgeVertexList.append(EdgeVertexList_i)

    PIPList.append(EdgeAList)
    PIPList.append(EdgeBList)
    PIPList.append(EdgeCOMList)
    PIPList.append(EdgexList)
    PIPList.append(EdgeyList)
    PIPList.append(EdgezList)
    PIPList.append(EdgeVertexList)
    return PIPList

def ConvexEdgesPlot(SimRobot, convex_edges_list, vis):
    Convex_Edges_Number = len(convex_edges_list)/2
    COMPos = SimRobot.getCom()
    for i in range(0, Convex_Edges_Number):
        EdgeA = convex_edges_list[2*i]
        EdgeB = convex_edges_list[2*i + 1]
        # Three edges to be added: A->B, A -> COM, B-> COM
        Edge_Index = str(i)
        vis.add("Edge:" + Edge_Index, Trajectory([0, 1], [EdgeA, EdgeB]))
        vis.hideLabel("Edge:" + Edge_Index, True)
        vis.setAttribute("Edge:" + Edge_Index,'width', 5.0)

def PIPVisualizer(i, EdgeA, EdgeB, EdgeCOM, Edgex, Edgey, Edgez, COMPos, vis):
    scale = 0.25
    Edge_Index = str(i)
    vis.add("PIPEdge:" + Edge_Index, Trajectory([0, 1], [EdgeA, EdgeB]))
    vis.hideLabel("PIPEdge:" + Edge_Index, True)
    vis.setAttribute("PIPEdge:" + Edge_Index,'width', 7.5)

    vis.add("PIPEdgeCOM:" + Edge_Index, Trajectory([0, 1], [COMPos, EdgeCOM]))
    vis.hideLabel("PIPEdgeCOM:" + Edge_Index, True)
    vis.setAttribute("PIPEdgeCOM:" + Edge_Index,'width', 7.5)
    vis.setColor("PIPEdgeCOM:" + Edge_Index, 65.0/255.0, 199.0/255.0, 244.0/255.0, 1.0)

    # Local Coordinates
    Edgex_i = [ 0.0, 0.0, 0.0]
    Edgex_i[0] = EdgeCOM[0] + scale * Edgex[0]
    Edgex_i[1] = EdgeCOM[1] + scale * Edgex[1]
    Edgex_i[2] = EdgeCOM[2] + scale * Edgex[2]
    vis.add("PIPEdgex:" + Edge_Index, Trajectory([0, 1], [EdgeCOM, Edgex_i]))
    vis.hideLabel("PIPEdgex:" + Edge_Index, True)
    vis.setAttribute("PIPEdgex:" + Edge_Index,'width', 7.5)
    vis.setColor("PIPEdgex:" + Edge_Index, 1.0, 0.0, 0.0, 1.0)

    Edgey_i = [ 0.0, 0.0, 0.0]
    Edgey_i[0] = EdgeCOM[0] + scale * Edgey[0]
    Edgey_i[1] = EdgeCOM[1] + scale * Edgey[1]
    Edgey_i[2] = EdgeCOM[2] + scale * Edgey[2]
    vis.add("PIPEdgey:" + Edge_Index, Trajectory([0, 1], [EdgeCOM, Edgey_i]))
    vis.hideLabel("PIPEdgey:" + Edge_Index, True)
    vis.setAttribute("PIPEdgey:" + Edge_Index,'width', 7.5)
    vis.setColor("PIPEdgey:" + Edge_Index, 155.0/255.0, 244.0/255.0, 66.0/255.0, 1.0)

    Edgez_i = [ 0.0, 0.0, 0.0]
    Edgez_i[0] = EdgeCOM[0] + scale * Edgez[0]
    Edgez_i[1] = EdgeCOM[1] + scale * Edgez[1]
    Edgez_i[2] = EdgeCOM[2] + scale * Edgez[2]
    # print Edgez_i
    vis.add("PIPEdgez:" + Edge_Index, Trajectory([0, 1], [EdgeCOM, Edgez_i]))
    vis.hideLabel("PIPEdgez:" + Edge_Index, True)
    vis.setAttribute("PIPEdgez:" + Edge_Index,'width', 7.5)
    vis.setColor("PIPEdgez:" + Edge_Index, 68.0/255.0, 65.0/255.0, 244.0/255.0, 1.0)

def PIPCleaner(i, vis):
    # This function is used to remove the previous PIP to make sure that no residual PIPs exist.
    Edge_Index = str(i)
    vis.hide("PIPEdge:" + Edge_Index, True)
    vis.hide("PIPEdgeCOM:" + Edge_Index, True)
    vis.hide("PIPEdgex:" + Edge_Index, True)
    vis.hide("PIPEdgey:" + Edge_Index, True)
    vis.hide("PIPEdgez:" + Edge_Index, True)

def RobotCOMPlot(SimRobot, vis):
    COMPos_start = SimRobot.getCom()
    COMPos_end = COMPos_start[:]
    COMPos_end[2] = COMPos_end[2] - 7.50
    vis.add("COM", Trajectory([0, 1], [COMPos_start, COMPos_end]))
    vis.hideLabel("COM",True)
    vis.setColor("COM", 0.0, 204.0/255.0, 0.0, 1.0)
    vis.setAttribute("COM",'width', 5.0)

def ContactDataGene(ReachableContacts_data):
    RowNo, ColumnNo = ReachableContacts_data.shape
    RowStart = 0
    RowEnd = RowNo

    point = []
    for i in range(RowStart, RowEnd):
        point_start = [0.0, 0.0, 0.0]
        ReachableContact_i = ReachableContacts_data[i]
        point_start[0] = ReachableContact_i[0]
        point_start[1] = ReachableContact_i[1]
        point_start[2] = ReachableContact_i[2]
        point.append(point_start)
    return ContactDataGene

def ContactDataUnplot(vis, StepNo, LimbNo, ReachableContacts_data):
    RowNo, ColumnNo = ReachableContacts_data.shape
    RowStart = 0
    RowEnd = RowNo
    for i in range(RowStart, RowEnd):
        ContactName = "Stage" + str(StepNo) + "LinkNo" + str(LimbNo) + "Point:" + str(i)
        vis.hide(ContactName, True)

def ContactDataPlot(vis, StepNo, LimbNo, ReachableContacts_data):
    RowNo, ColumnNo = ReachableContacts_data.shape
    RowStart = 0
    RowEnd = RowNo

    for i in range(RowStart, RowEnd):
        point_start = [0.0, 0.0, 0.0]
        ReachableContact_i = ReachableContacts_data[i]
        point_start[0] = ReachableContact_i[0]
        point_start[1] = ReachableContact_i[1]
        point_start[2] = ReachableContact_i[2]
        ContactName = "Stage" + str(StepNo) +"LinkNo" + str(LimbNo) + "Point:" + str(i)
        vis.add(ContactName, point_start)
        vis.hideLabel(ContactName, True)
        vis.setColor(ContactName, 65.0/255.0, 199.0/255.0, 244.0/255.0, 1.0)

def TransitionDataUnplot(vis, StepNo, LimbNo):
    TransName = "Stage" + str(StepNo) + "LinkNo" + str(LimbNo) + "Path"
    vis.hide(TransName, True)

def TransitionDataPlot(vis, StepNo, LimbNo, ReachableContacts_data):
    RowNo, ColumnNo = ReachableContacts_data.shape
    RowStart = 0
    RowEnd = RowNo
    Traj = []
    for i in range(RowStart, RowEnd):
        point_i = [0.0, 0.0, 0.0]
        ReachableContact_i = ReachableContacts_data[i]
        point_i[0] = ReachableContact_i[0]
        point_i[1] = ReachableContact_i[1]
        point_i[2] = ReachableContact_i[2]
        Traj.append(point_i)
    TransName = "Stage" + str(StepNo) + "LinkNo" + str(LimbNo) + "Path"
    vis.add(TransName, Trajectory([0, 1], Traj))
    vis.hideLabel(TransName, True)
    vis.setColor(TransName, 255.0/255.0, 255.0/255.0, 51.0/255.0, 1.0)

def ContactDataLoader(Specificpath, StepNo, LimbNo, IdealReachableContact):
    IdealReachableContacts = Specificpath+ "/" + str(StepNo) + "_" + str(LimbNo) + "_" + IdealReachableContact + ".bin"
    f_IdealReachableContacts = open(IdealReachableContacts, 'rb')
    IdealReachableContacts_data = np.fromfile(f_IdealReachableContacts, dtype=np.double)
    IdealReachableContacts_data = IdealReachableContacts_data.reshape((IdealReachableContacts_data.size/3, 3))
    return IdealReachableContacts_data

def ImpulseInfoReader(file_path):
    PushInfoFilePath = file_path + "/PushInfoFile.txt"
    with open(PushInfoFilePath,'r') as PushInfoFilePathFile:
        PushInfoTot = PushInfoFilePathFile.readlines()
        # Actually we only need two information: Time and Impulsive Force
        StartPushInfoTot = PushInfoTot[0].split(" ")
        StartTime = float(StartPushInfoTot[0])
        EndPushInfoTot = PushInfoTot[len(PushInfoTot)-1].split(" ")
        EndTime = float(EndPushInfoTot[0])
        ImpulseX = float(EndPushInfoTot[1])
        ImpulseY = float(EndPushInfoTot[2])
        ImpulseZ = float(EndPushInfoTot[3].replace('\n',''))
    return StartTime, EndTime, [ImpulseX, ImpulseY, ImpulseZ]

def ImpulsePlot(vis, SimRobot, SimulationTime, ImpulseObj):
    StartTime = ImpulseObj[0]
    EndTime = ImpulseObj[1]
    ImpulForce = ImpulseObj[2]
    DurationTime = EndTime - StartTime

    OriPoint = SimRobot.link(19).getWorldPosition([0.0, 0.0, 0.0])
    ImpulseScale = 0.001 * (SimulationTime - StartTime)/DurationTime
    EndPoint = OriPoint[:]
    EndPoint[0]+=ImpulForce[0] * ImpulseScale
    EndPoint[1]+=ImpulForce[1] * ImpulseScale
    EndPoint[2]+=ImpulForce[2] * ImpulseScale

    if SimulationTime<StartTime:
        pass
    elif SimulationTime>EndTime:
        vis.hide("Impulse", True)
    else:
        vis.add("Impulse", Trajectory([0, 1], [OriPoint, EndPoint]))
        vis.hideLabel("Impulse",True)
        vis.setColor("Impulse", 235.0/255.0, 52.0/255.0, 52.0/255.0, 1.0)
        vis.setAttribute("Impulse",'width', 5.0)

def PlanningInfoReader(file_path):
    PlanningInfoFilePath = file_path + "/PlanningInfoFile.txt"
    PlanningTimeList = []
    PlanningResList = []
    with open(PlanningInfoFilePath,'r') as PlanningInfoFile:
        PlanningInfoTot = PlanningInfoFile.readlines()
        # Actually we only need two information: Time and Impulsive Force
        for PlanningInfo_i in PlanningInfoTot:
            PlanningInfo_i_Str = PlanningInfo_i.split(" ")
            PlanningTime_i = float(PlanningInfo_i_Str[2].replace('\n',''))
            PlanningTimeList.append(PlanningTime_i)
            Step = int(PlanningInfo_i_Str[0])
            Limb = int(PlanningInfo_i_Str[1])
            PlanningResList.append([Step, Limb])
    return PlanningTimeList, PlanningResList

def PlanningInterval(SimulationTime, PlanningTimeList, TimeStep):
    if (SimulationTime>max(PlanningTimeList)):
        return len(PlanningTimeList)-1
    else:
        for i in range(0, len(PlanningTimeList)):
            if (SimulationTime<(PlanningTimeList[i]+TimeStep)):
                if(i>0):
                    return i-1
                else:
                    return i
def ContactDataRefine(ContactPts, ContactWeights_array):
    # This function is used to address the Klampt visualization problem
    ContactPointSize = ContactWeights_array.size/3
    if ContactPointSize>100:
        ContactWeightLists = []
        for i in range(ContactPointSize):
            ContactWeight_i = ContactWeights_array[i]
            ContactWeight_i_value = ContactWeight_i[0]**2 + ContactWeight_i[1]**2 + ContactWeight_i[2]**2
            ContactWeightLists.append(ContactWeight_i_value)
        sorted_indices = sorted(range(len(ContactWeightLists)), key=lambda k: ContactWeightLists[k])
        sorted_indices.reverse()
        ContactPts_ = []
        ContactWeights_ = []
        for i in range(0, 20):
            sorted_index = sorted_indices[i]
            ContactPt_i = ContactPts[sorted_index]
            ContactWeight_i = ContactWeights_array[sorted_index]
            ContactPts_.append(ContactPt_i.tolist())
            ContactWeights_.append(ContactWeight_i.tolist())
        return np.array(ContactPts_), np.array(ContactWeights_)
    else:
        return ContactPts, ContactWeights_array

def WeightedContactDataPlot(vis, StepNo, LimbNo, OptimalContact_data, OptimalContactWeights_data):
    scale = 1.0
    for i in range(OptimalContact_data.size/3):
        point_start = [0.0, 0.0, 0.0]
        ReachableContact_i = OptimalContact_data[i]
        point_start[0] = ReachableContact_i[0]
        point_start[1] = ReachableContact_i[1]
        point_start[2] = ReachableContact_i[2]

        point_end = [0.0, 0.0, 0.0]
        ReachableContactWeight_i = OptimalContactWeights_data[i]
        point_end[0] = point_start[0] + scale * ReachableContactWeight_i[0]
        point_end[1] = point_start[1] + scale * ReachableContactWeight_i[1]
        point_end[2] = point_start[2] + scale * ReachableContactWeight_i[2]

        ContactName = "Stage" + str(StepNo) +"LinkNo" + str(LimbNo) + "Point:" + str(i)
        print i
        vis.add(ContactName, Trajectory([0, 1], [point_start, point_end]))
        vis.hideLabel(ContactName, True)
        vis.setColor(ContactName, 0.0, 204.0/255.0, 0.0, 1.0)
        vis.setAttribute(ContactName, 'width', 5.0)

def WeightedContactDataUnPlot(vis, StepNo, LimbNo, OptimalContact_data):
    for i in range(OptimalContact_data.size/3):
        ContactName = "Stage" + str(StepNo) +"LinkNo" + str(LimbNo) + "Point:" + str(i)
        vis.hide(ContactName, True)

def EndEffectorTrajPlot(vis, SimulationTime, PlanningObj, SpecificPath, PlotDuration, TimeStep):
    PlanningTimeList = PlanningObj[0]
    PlanningPairList = PlanningObj[1]
    if SimulationTime<min(PlanningTimeList):
        pass
    elif SimulationTime>(max(PlanningTimeList) + PlotDuration):
        StepLimbPair = PlanningPairList[len(PlanningPairList)-1]
        StepNo = StepLimbPair[0]
        LimbNo = StepLimbPair[1]
        for i in range(0, LimbNo+1):
            CandidateContacts_data = ContactDataLoader(SpecificPath, StepNo, i, "CandidateContacts")
            CandidateContactWeights_data = ContactDataLoader(SpecificPath, StepNo, i, "CandidateContactWeights")
            CandidateContacts_data, CandidateContactWeights_data = ContactDataRefine(CandidateContacts_data, CandidateContactWeights_data)

            WeightedContactDataPlot(vis, StepNo, LimbNo, CandidateContacts_data, CandidateContactWeights_data)

            PathWaypoints_data = ContactDataLoader(SpecificPath, StepNo, i, "PathWaypoints")
            TransitionDataPlot(vis, StepNo, LimbNo, PathWaypoints_data)

            # WeightedContactDataUnPlot(vis, CandidateContacts_data)
            # TransitionDataUnplot(vis, StepNo, i)
    else:
        PlanIndex = PlanningInterval(SimulationTime, PlanningTimeList, TimeStep)
        StepLimbPair = PlanningPairList[PlanIndex]
        StepNo = StepLimbPair[0]
        LimbNo = StepLimbPair[1]
        for i in range(0, LimbNo+1):
            CandidateContacts_data = ContactDataLoader(SpecificPath, StepNo, i, "CandidateContacts")
            CandidateContactWeights_data = ContactDataLoader(SpecificPath, StepNo, i, "CandidateContactWeights")
            CandidateContacts_data, CandidateContactWeights_data = ContactDataRefine(CandidateContacts_data, CandidateContactWeights_data)

            WeightedContactDataPlot(vis, StepNo, i, CandidateContacts_data, CandidateContactWeights_data)

            PathWaypoints_data = ContactDataLoader(SpecificPath, StepNo, i, "PathWaypoints")
            TransitionDataPlot(vis, StepNo, i, PathWaypoints_data)
        if PlanIndex>0:
            for i in range(0, PlanIndex):
                StepLimbPair = PlanningPairList[i]
                StepNo = StepLimbPair[0]
                LimbNo = StepLimbPair[1]
                for j in range(0, LimbNo+1):
                    CandidateContacts_data = ContactDataLoader(SpecificPath, StepNo, j, "CandidateContacts")
                    CandidateContactWeights_data = ContactDataLoader(SpecificPath, StepNo, j, "CandidateContactWeights")
                    CandidateContacts_data, CandidateContactWeights_data = ContactDataRefine(CandidateContacts_data, CandidateContactWeights_data)

                    WeightedContactDataUnPlot(vis, StepNo, j, CandidateContacts_data)
                    TransitionDataUnplot(vis, StepNo, j)

def ExperimentVisualizer(world, ContactLinkDictionary, ExpTraj, PIPInfoList, ImpulseObj, PlanningObj, SpecificPath, Para):

    ExpViewer = MyGLPlugin(world)
    vis.pushPlugin(ExpViewer)
    vis.add("world", world)
    vis.show()

    SimRobot = world.robot(0)

    FailureStateTraj = ExpTraj[0]
    CtrlStateTraj = ExpTraj[1]
    PlanStateTraj = ExpTraj[2]

    EdgeAList = PIPInfoList[0]
    EdgeBList = PIPInfoList[1]
    EdgeCOMList = PIPInfoList[2]
    EdgexList = PIPInfoList[3]
    EdgeyList = PIPInfoList[4]
    EdgezList = PIPInfoList[5]
    EdgeVertexList = PIPInfoList[6]

    PlanningTimeList = PlanningObj[0]
    PlanningResList = PlanningObj[1]

    TimeStep = CtrlStateTraj.times[1] - CtrlStateTraj.times[0]
    PlotDuration = 10 * TimeStep
    StateTrajLength = len(CtrlStateTraj.times)
    PIPTrajLength = len(EdgeAList)              # Here PIPTraj could be less than the length of state

    StateTraj = []

    StateType = Para[0]
    VisMode = Para[1]

    if StateType == "F":
        print "Failure State Traj!"
        StateTraj = FailureStateTraj.milestones
    elif StateType == "C":
        print "Controlled State Traj!"
        StateTraj = CtrlStateTraj.milestones
    else:
        print "Planned State Traj!"
        StateTraj = PlanStateTraj.milestones

    SimulationTime = 0.0
    while vis.shown():
        # This is the main plot program
        for i in range(0, StateTrajLength):
            SimulationTime = 1.0 * i * TimeStep
            vis.lock()
            Config = StateTraj[i]
            SimRobot.setConfig(Config)
            ImpulsePlot(vis, SimRobot, SimulationTime, ImpulseObj)

            if VisMode == "PIP":
                COMPos = SimRobot.getCom()
                if (i >= (StateTrajLength - PIPTrajLength)):
                    EdgeIndex = i + PIPTrajLength - StateTrajLength
                    EdgeAList_i = EdgeAList[EdgeIndex]
                    EdgeBList_i = EdgeBList[EdgeIndex]
                    EdgeCOMList_i = EdgeCOMList[EdgeIndex]
                    EdgexList_i = EdgexList[EdgeIndex]
                    EdgeyList_i = EdgeyList[EdgeIndex]
                    EdgezList_i = EdgezList[EdgeIndex]
                    EdgeVertexList_i = EdgeVertexList[EdgeIndex]
                    for j in range(0, len(EdgeAList_i)):
                        EdgeA = EdgeAList_i[j]
                        EdgeB = EdgeBList_i[j]
                        EdgeCOM = EdgeCOMList_i[j]
                        Edgex = EdgexList_i[j]
                        Edgey = EdgeyList_i[j]
                        Edgez = EdgezList_i[j]
                        PIPVisualizer(j, EdgeA, EdgeB, EdgeCOM, Edgex, Edgey, Edgez, COMPos, vis)
            elif VisMode == "Poly":
                FeasiFlag = 1        # for j in range(0, len(EdgeAList_i)):
                if (i >= (StateTrajLength - PIPTrajLength)):
                    EdgeIndex = i + PIPTrajLength - StateTrajLength
                    EdgeVertexList_i = EdgeVertexList[EdgeIndex]
                    try:
                        h = ConvexHull(EdgeVertexList_i)
                    except:
                        FeasiFlag = 0
                    if FeasiFlag == 1:
                        h = ConvexHull(EdgeVertexList_i)
                        hrender = draw_hull.PrettyHullRenderer(h)
                        vis.add("ContactPolytope", h)
                        vis.setDrawFunc("ContactPolytope", my_draw_hull)
                    else:
                        print "Input Contact Polytope Infeasible!"
            else:
                EndEffectorTrajPlot(vis, SimulationTime, PlanningObj, SpecificPath, PlotDuration, TimeStep)
            vis.unlock()
            time.sleep(1.0 * TimeStep)
            if VisMode == "PIP" and i >= (StateTrajLength - PIPTrajLength):
                EdgeIndex = i + PIPTrajLength - StateTrajLength
                for i in range(0, len(EdgeAList[EdgeIndex])):
                    PIPCleaner(i, vis)

        if VisMode == "Traj":
            PlanningPairList = PlanningObj[1]
            StepLimbPair = PlanningPairList[len(PlanningPairList)-1]
            StepNo = StepLimbPair[0]
            LimbNo = StepLimbPair[1]
            for i in range(0, LimbNo+1):
                CandidateContacts_data = ContactDataLoader(SpecificPath, StepNo, i, "CandidateContacts")
                CandidateContactWeights_data = ContactDataLoader(SpecificPath, StepNo, i, "CandidateContactWeights")

                CandidateContacts_data, CandidateContactWeights_data = ContactDataRefine(CandidateContacts_data, CandidateContactWeights_data)

                WeightedContactDataUnPlot(vis, StepNo, i, CandidateContacts_data)
                TransitionDataUnplot(vis, StepNo, i)
        if VisMode == "Poly":
            vis.hide("ContactPolytope", True)

def main():
    ExpNo = 0
    StateType = "P"
    VisMode = "Traj"
    if(len(sys.argv[1:])==1):
        ExpNo = int(sys.argv[1:][0])
    elif(len(sys.argv[1:])==2):
        ExpNo = int(sys.argv[1:][0])
        StateType = sys.argv[1:][1]
    else:
        ExpNo = int(sys.argv[1:][0])
        StateType = sys.argv[1:][1]
        VisChoice = str(sys.argv[1:][2])
        if VisChoice == "T" or VisChoice == "Traj":
            VisMode = "Traj"
        if VisChoice == "Poly" or VisChoice == "Polytope":
            VisMode = "Poly"
        if VisChoice == "PIP":
            VisMode = "PIP"

    world = WorldModel()                    	# WorldModel is a pre-defined class
    curDir = os.getcwd()
    CurCasePath = curDir[0:-4] + "-Data/result/" + CurCase

    XML_path = CurCasePath + "/Environment.xml"
    result = world.readFile(XML_path)         	# Here result is a boolean variable indicating the result of this loading operation
    if not result:
        raise RuntimeError("Unable to load model " + XML_path)
    ContactLinkDictionary = ContactLinkReader(curDir + "/../user/ContactLink.txt")
    PlanStateTraj = Trajectory(world.robot(0))
    CtrlStateTraj = Trajectory(world.robot(0))
    FailureStateTraj = Trajectory(world.robot(0))
    SpecificPath = CurCasePath + "/" + str(ExpNo)
    PlanStateTraj.load(SpecificPath + "/PlanStateTraj.path")
    CtrlStateTraj.load(SpecificPath+ "/CtrlStateTraj.path")
    FailureStateTraj.load(SpecificPath + "/FailureStateTraj.path")
    ExpTraj = [FailureStateTraj, CtrlStateTraj, PlanStateTraj]
    PIPInfoList = PIPTrajReader(SpecificPath)
    StartTime, EndTime, ImpulForce = ImpulseInfoReader(SpecificPath)
    ImpulseObj = [StartTime, EndTime, ImpulForce]
    PlanningTimeList, PlanningResList = PlanningInfoReader(SpecificPath)
    PlanningObj = [PlanningTimeList, PlanningResList]
    ExperimentVisualizer(world, ContactLinkDictionary, ExpTraj, PIPInfoList, ImpulseObj, PlanningObj, SpecificPath, [StateType, VisMode])

if __name__ == "__main__":
    # StateType = "F"      # "F", "C", "P"
    main()
