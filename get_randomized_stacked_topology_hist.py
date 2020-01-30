mport glob
import os
import re
import sys 
import time

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 
import seaborn as sns 

import process_data as PD
import graphs

def concat_data(n,rand_true):
    DATA_DIR = "/home/somya/Somya/version1/data_linmaps/"
    if rand_true = 0:
        files3 = [x for x in sorted(glob.glob(DATA_DIR + 'N3/*.pickle')) if 'random' not in x][:n]
        files4 = [x for x in sorted(glob.glob(DATA_DIR + 'N4/*.pickle')) if 'random' not in x][:n]
        files5 = [x for x in sorted(glob.glob(DATA_DIR + 'N5/*.pickle')) if 'random' not in x][:n]
        files6 = [x for x in sorted(glob.glob(DATA_DIR + 'N6/*.pickle')) if 'random' not in x][:n]
        files7 = [x for x in sorted(glob.glob(DATA_DIR + 'N7/*.pickle')) if 'random' not in x][:n]
    else:
        files3 = sorted(glob.glob(DATA_DIR + 'N3/*.pickle'))[:n]
        files4 = sorted(glob.glob(DATA_DIR + 'N4/*.pickle'))[:n]
        files5 = sorted(glob.glob(DATA_DIR + 'N5/*.pickle'))[:n]
        files6 = sorted(glob.glob(DATA_DIR + 'N6/*.pickle'))[:n]
        files7 = sorted(glob.glob(DATA_DIR + 'N7/*.pickle'))[:n]

    
    columns_load = ['n_gene','n_nodes','n_edges','num_components','f_asy','f_sig','f_adj','unicell','scc','chain','tree','dag','cyclic']
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
    df_big3456 = pd.concat([df_big345,df_big6],ignore_index=True)
    df_big6 = None
    df_big345 = None

    df_big7 = pd.concat([pd.read_pickle(f) for f in files7], ignore_index=True)
    df_big7 = df_big7.loc[:,columns_load]
    l7 = len(df_big7)
    print('N7 loaded')
    df = pd.concat([df_big3456,df_big7],ignore_index=True)
    df_big7 = None
    df['N'] = [3]*l3 + [4]*l4 + [5]*l5 + [6]*l6 +[7]*l7

    return df

def stacked_burger_plot(df):

# get stacked bar plot for graph topology abundance
    order = ['u', 's', 'C', 'c', 'd', 't']
    top = ['unicell','scc','cyc','chain','dag','tree']
    topology = np.array(['']*len(df))
    for i,o in enumerate(top):
        topology[np.where(df[o]==1)[0]]=order[i]# only takes in a single alphabet
    df['topology'] = topology

    n_counts = [Counter(df.loc[df.N==n, 'topology']) for n in range(3,8)]

    X = range(3,8)
    Y = np.array([[n_counts[i][o] for i in range(5)] for o in order], dtype=float)
    for i in range(5):
        Y[:,i] = Y[:,i] / Y[:,i].sum()
    width = 0.75

    col = list(np.array(Set1_9.mpl_colors)[[6, 3, 2, 0, 1]]) + ['k']
    fig, ax = plt.subplots()
    bottom = np.zeros(5)
    for i in range(len(order)):
        ax.bar(X, Y[i], width, bottom=bottom, color=col[i], edgecolor='white', linewidth=0.5)
        bottom += Y[i]
    plt.xticks([3,4,5,6,7],('3','4','5','6','7'))
    plt.xlim((2.5,7.5))

