from collections import Counter
import glob
import os
import re
import sys 
import time
from itertools import product
import multiprocessing as mp
import itertools
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 
import seaborn as sns 
import networkx as nx
import jitterplot
from palettable.colorbrewer.qualitative import Set1_9
import process_data as PD
import graphs


def get_traj_df():
    TRAJ_DIR = "/home/version1/trajectory_lengths/" # ADD YOUR OWN DATA DIRECTORY
    files = sorted(glob.glob(TRAJ_DIR + 'trajectory_lengths_N*_*_genome*.txt'))
    Traj = pd.read_csv(files[0],header = None)
    Traj.columns = ["N","genomeID","initID","sig","asym",'adj','graph_nodes',"indep","intrinsic_indep","numsteps"]
    Traj['filenum'] = [0]*len(Traj)
    for i in np.linspace(0,len(files)-2,len(files)-1):
        Traj1 = pd.read_csv(files[int(i)+1],header = None)
        Traj1.columns = ["N","genomeID","initID","sig","asym",'adj','graph_nodes',"indep","intrinsic_indep","numsteps"]
        Traj1['filenum'] = [i+1]*len(Traj1)
        Traj = pd.concat([Traj,Traj1],ignore_index=True)
        print(i)

    Traj = Traj.drop_duplicates()
    Traj = Traj.sort_values(['filenum','N','genomeID','initID','sig','asym','adj'],ascending=[True,True,True,True,True,True,True])
    Traj = Traj.reset_index(drop=True)
    return Traj
