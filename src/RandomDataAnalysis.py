import sys, os, time
import math
import numpy as np

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

if __name__ == "__main__":
    main()
