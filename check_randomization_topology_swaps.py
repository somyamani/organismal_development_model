import argparse
import glob
import os
import re
import sys 
import time
import feather

from itertools import product
import multiprocessing as mp
import numpy as np
from networkx.algorithms import approximation as approx
import networkx as nx
import pandas as pd


# load original data and randomized data. Out of the graphs that were randomized, how many graphs changed topology?

def load_original_data(n,ver):
    DATA_DIR = "/home/version%d/data_linmaps/"%(ver) # ADD YOUR OWN DATA DIRECTORY
    files3 = sorted([x for x in glob.glob(DATA_DIR + 'N3/*.pickle') if '_random' not in x])[:n]
    files4 = sorted([x for x in glob.glob(DATA_DIR + 'N4/*.pickle') if '_random' not in x])[:n]
    files5 = sorted([x for x in glob.glob(DATA_DIR + 'N5/*.pickle') if '_random' not in x])[:n]
    files6 = sorted([x for x in glob.glob(DATA_DIR + 'N6/*.pickle') if '_random' not in x])[:n]
    files7 = sorted([x for x in glob.glob(DATA_DIR + 'N7/*.pickle') if '_random' not in x])[:n]
    columns_load = ['num_components','unicell','scc','chain','tree','dag','cyclic']

    df_big3 = pd.concat([pd.read_pickle(f) for f in files3], ignore_index=True)
    df_big3 = df_big3.loc[:,columns_load]
    l3 = len(df_big3)
    print('N3 loaded')

    df_big4 = pd.concat([pd.read_pickle(f) for f in files4], ignore_index=True)
    df_big4 = df_big4.loc[:,columns_load]
    l4 = len(df_big4)
    print('N4 loaded')

    df_big34 =pd.concat([df_big3,df_big4], ignore_index=True)
    df_big3 = None
    df_big4 = None

    df_big5 = pd.concat([pd.read_pickle(f) for f in files5], ignore_index=True)
    df_big5 = df_big5.loc[:,columns_load]
    l5 = len(df_big5)
    print('N5 loaded')

    df_big345 = pd.concat([df_big34,df_big5],ignore_index=True)
    df_big5 = None
    df_big34 = None
    
    df_big6 = pd.concat([pd.read_pickle(f) for f in files6], ignore_index=True)
    df_big6 = df_big6.loc[:,columns_load]
    l6 = len(df_big6)
    print('N6 loaded')
    df = pd.concat([df_big345,df_big6],ignore_index=True)
    df_big6 = None
    df_big345 = None
    df['N'] = ['N3']*l3 + ['N4']*l4 + ['N5']*l5 + ['N6']*l6
    return df

def load_randomized_data(n,ver):
    DATA_DIR = "/home/somya/Somya/version%d/data_linmaps/"%(ver)
    files3 = sorted(glob.glob(DATA_DIR + 'N3/*_random2.pickle'))[:n]
    files4 = sorted(glob.glob(DATA_DIR + 'N4/*_random2.pickle'))[:n]
    files5 = sorted(glob.glob(DATA_DIR + 'N5/*_random2.pickle'))[:n]
    files6 = sorted(glob.glob(DATA_DIR + 'N6/*_random2.pickle'))[:n]
    files7 = sorted(glob.glob(DATA_DIR + 'N7/*_random2.pickle'))[:n]
    columns_load = ['num_components','unicell','scc','chain','tree','dag','cyclic','num_changes']

    df_big3 = pd.concat([pd.read_pickle(f) for f in files3], ignore_index=True)
    df_big3 = df_big3.loc[:,columns_load]
    l3 = len(df_big3)
    print('N3 loaded')

    df_big4 = pd.concat([pd.read_pickle(f) for f in files4], ignore_index=True)
    df_big4 = df_big4.loc[:,columns_load]
    l4 = len(df_big4)
    print('N4 loaded')

    df_big34 =pd.concat([df_big3,df_big4], ignore_index=True)
    df_big3 = None
    df_big4 = None

    df_big5 = pd.concat([pd.read_pickle(f) for f in files5], ignore_index=True)
    df_big5 = df_big5.loc[:,columns_load]
    l5 = len(df_big5)
    print('N5 loaded')
    df_big345 = pd.concat([df_big34,df_big5],ignore_index=True)
    df_big5 = None
    df_big34 = None

    df_big6 = pd.concat([pd.read_pickle(f) for f in files6], ignore_index=True)
    df_big6 = df_big6.loc[:,columns_load]
    l6 = len(df_big6)
    print('N6 loaded')
    df = pd.concat([df_big345,df_big6],ignore_index=True)
    df_big6 = None
    df_big345 = None
    df['N'] = ['N3']*l3 + ['N4']*l4 + ['N5']*l5 + ['N6']*l6
    return df

def get_conversions(df_ori,df_rand):
    df_to_acyclic =  df_rand.loc[((df_ori.chain + df_ori.tree + df_ori.dag) == 0)&(df_rand.num_changes > 0), ['chain','tree','dag']]
    df_to_cyclic =  df_rand.loc[((df_ori.cyclic + df_ori.scc) == 0)&(df_rand.num_changes > 0), ['scc','cyclic']]
    f_acyclic = sum(df_to_acyclic.chain + df_to_acyclic.tree + df_to_acyclic.dag)/len(df_to_acyclic)
    f_cyclic = sum(df_to_cyclic.cyclic + df_to_cyclic.scc)/len(df_to_cyclic)
    return f_acyclic, f_cyclic, len(df_to_acyclic), len(df_to_cyclic)

def get_graphtype_conversion_fractions(df_ori, df_rand):
    # unicell, scc, cyclic, chain, tree, dag
    gtypes = ['unicell', 'scc','chain','tree','dag','cyclic']
    conv_mat = np.zeros((6,6))
    for i in range(len(gtypes)):
        for j in range(len(gtypes)):
            X = gtypes[i]
            Y = gtypes[j]
            a = sum(df_rand.loc[(df_ori[X] == 1), 'num_changes'] > 0)
            conv_mat[i,j] = 100 * sum(df_rand.loc[(df_ori[X] == 1)&(df_ori.num_components == 1)&(df_rand.num_changes > 0),Y])/a
    return conv_mat


if __name__ == "__main__":
    df_ori = load_original_data(10)
    df_rand = load_randomized_data(10)
   #f1,f2, n1,n2 = get_conversions(df_ori,df_rand)
   #print(f1*100,n1)
   #print(f2*100,n2)
    conv_mat = get_graphtype_conversion_fractions(df_ori, df_rand)
    sns.heatmap(conv_mat)
     


 
