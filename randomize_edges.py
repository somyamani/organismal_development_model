import argparse
import glob
import os
import re
import sys 
import time

from itertools import product
import multiprocessing as mp
import numpy as np
from networkx.algorithms import approximation as approx
import networkx as nx
import pandas as pd
from random import randint
import process_data


N_PROC = 55
CHUNK=100

N_BATCH = 2 

DATA_DIR = "/home/version1/data_linmaps/"
# INPUT YOUR OWN DATA DIRECTORY ABOVE

parser = argparse.ArgumentParser()
parser.add_argument('-f','--fName', help='Filename to be processed - full path required', dest='fName', default='')
ARGS = parser.parse_args()


def load_single_file(fName):
  #  cols = ['n_gene','n_fixed', 'n_osc', 'basin', 'ss', 'f_sig', 'm_sig', 'f_adj', 'm_adj', 'f_asy', 'm_asy', 'edges']
    ext = os.path.splitext(fName)[1]
    return pd.read_pickle(fName)

def int_arr_to_str(int_arr, delim=';'):
    return delim.join([str(x) for x in int_arr])

def randomize_edges(edges):
    edges1 = edges
    connections = np.array([int(x) for x in edges.split(';')], dtype=int)
    if connections.size <= 2:
        return edges , False
    edges = connections.reshape(int(connections.size/2), 2)
    str_edges1 = np.array([';'.join([str(x) for x in e]) for e in edges])
    set_edges1 = set(str_edges1)
    count = 0 
    while True:
        count += 1
        rand_idx = np.arange(edges.shape[0])[np.argsort(np.random.rand(edges.shape[0]))]
        edges[:,1] = edges[rand_idx,1]

        str_edges = np.array([';'.join([str(x) for x in e]) for e in edges])
        set_edges = set(str_edges)
        
        if len(set_edges) == edges.shape[0]:
            return ';'.join([str(x) for x in edges.ravel()]), set_edges1 != set_edges
        if count == 1000:
            return edges1 , False # if randomization does not work in 1000 goes, return the original graph
#   else:
#       count = np.array([len(re.findall(s, one_str)) for s in str_edges])
#       multi_idx = np.where(count>1)


def process_database(fName):

    timeS = time.time()
    # Load csv file
    df = load_single_file(fName)
    print(f"Time taken to load: {(time.time()-timeS)/60.} mins")
    # Remove results that did not converge
    df = df.drop(index=df.loc[df.edges.str.len()==1].index).reset_index(drop=True)
    df2 = pd.DataFrame()
    
    #set up POOL
    pool = mp.Pool(N_PROC)
    # randomise edges
    out = np.array(list(pool.map(randomize_edges, df.edges,CHUNK)))
    df2['edges'] = out[:,0]
    
    print(f"Time taken to get random graphs: {(time.time()-timeS)/60.} mins")
    ### Close POOL
    pool.close()
    # Get random graph properties
    df2 = process_data.update_dataframe_with_graph_properties(df2)
    print(f"Time taken to get random graph properties: {(time.time()-timeS)/60.} mins")
    # Save in high density format
    new_name = os.path.splitext(fName)[0] + '_random.pickle'
    #df.to_hdf(new_name, 'df', table=True, mode='w', complib='blosc')
    df2.to_pickle(new_name)
    print(f"Time taken to write: {(time.time()-timeS)/60.} mins")
    return 0    

if __name__ == "__main__":

    timeS = time.time()

    process_database(ARGS.fName)

    print(f"Time to run...  {(time.time()-timeS)/60.} minutes")

