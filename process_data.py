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

DATA_DIR = "/home/somya/data_linmaps/"
OUT_DIR = "/home/somya/data_linmaps/compiled_data/"

# INPUT YOUR OWN DATA DIRECTORY AND OUTPUT DIRECTORY


N_PROC=55
CHUNK=100

N_BATCH = 2

CELL_TYPES = {n: np.array([[int(x) for x in list(np.binary_repr(i, width=n))] for i in range(2**n)]) for n in range(3,9)}


parser = argparse.ArgumentParser()
parser.add_argument('-f','--fName', help='Filename to be processed - full path required', dest='fName', default='')
ARGS = parser.parse_args()

def load_hdf_and_amalgamate(files):
    return pd.concat([pd.read_hdf(f) for f in files], ignore_index=True)

def amalgamate_csv_files_to_dataframe():
    files = sorted(glob.glob(DATA_DIR + 'N3/*.txt'))
    cols = ['n_gene','n_fixed', 'n_osc', 'basin', 'ss', 'f_sig', 'm_sig', 'f_adj', 'm_adj', 'f_asy', 'm_asy', 'edges']
    return pd.concat([pd.read_csv(f, names=cols) for f in files], ignore_index=True)

def int_arr_to_str(int_arr, delim=';'):
    return delim.join([str(x) for x in int_arr])

def add_column(df, X, Y, some_function):
#   df['new_parameter'] = df['edges']
#   df['new_parameter'] = df.edges.apply(some_function)
    df[Y] = df.edges.apply(some_function)
    return df

def add_column_parallel(df, X, Y, some_function):
    pool = mp.Pool(N_PROC)
    df[Y] = pool.map(some_function, df[X], CHUNK)
    pool.close()
    return df

def fn_n_nodes(x):
    # find number of nodes in input graph (edge-list)
    return len(set(x.split(';')))

def fn_n_edges(x):
    # find number of edges in input graph (edge-list)
    return (len(re.findall(';', x))+1)/2

def get_graph_properties(edges):
    # Set up graph
    connections = np.array([int(x) for x in edges.split(';')])
    
    nodes = sorted(list(set(connections)))
    # Calculate Properties
    properties = []
    timings = {}

    if connections[0] > 0:
        edges = connections.reshape(int(connections.size/2), 2)
        timeS = time.time()        

        # directed graph
        G = nx.DiGraph()
        G.add_edges_from(edges)

        # undirected graph
        U = nx.Graph()
        U.add_edges_from(edges)
        # graph generated    
        

        # property 1: number of components
        num_comp = nx.number_connected_components(U)
        properties.append(num_comp)
    
        # property 2: number of strongly connected components
        num_strong_comp = nx.number_strongly_connected_components(G)
        properties.append(num_strong_comp)
        
        # property 3: average in/out degree
        indeg = [];
        outdeg = [];
        indeg_ls = list(G.in_degree())
        outdeg_ls = list(G.out_degree())
        
        for x in np.arange(len(nodes)):        
            indeg.append(indeg_ls[x][1])
            outdeg.append(outdeg_ls[x][1])
        av_deg = np.mean(indeg)
        properties.append(av_deg)
        
        # property 4: link density
        linkden = connections.size/(len(nodes) * len(nodes))
        properties.append(linkden)
        
        # property 5: number of self loops
        numloop = list(G.selfloop_edges())
        numloop = len(numloop)
        properties.append(numloop) 
#       # property 6: number of simple cycles (excluding self loops)
#       numcyc = list(nx.simple_cycles(G))
#       numcyc = len(numcyc) - numloop
#       properties.append(numcyc)

#       timings.update({'p6':time.time()-timeS})
#       print('p6')
#       print(timings['p6'])
#       timeS = time.time()        

        # find all components
        components = list(nx.connected_components(U));
        
        ischain = [None]*len(components)
        istree = [None]*len(components)
        isdag = [None]*len(components)
        unicel = [None]*len(components)
        isscc = [None]*len(components)
        iscyc = [None]*len(components)
        iseul = [None]*len(components)
        indeg_by_comp = []
        outdeg_by_comp = []
        node_conn = [0]*len(components)
        av_clust = [0.]*len(components)
        assort = [0.]*len(components)
        indeg_cen_av = [0.]*len(components)
        indeg_cen_max = [0.]*len(components)
        indeg_cen_min = [0.]*len(components)
        outdeg_cen_av = [0.]*len(components)
        outdeg_cen_max = [0.]*len(components)
        outdeg_cen_min = [0.]*len(components)
        bet_cen_av = [0.]*len(components)
        bet_cen_max = [0.]*len(components)
        bet_cen_min = [0.]*len(components)
        eig_cen_av = [0.]*len(components)
        eig_cen_max = [0.]*len(components)
        eig_cen_min = [0.]*len(components)
        triangles_av = [0.]*len(components)
        triangles_max = [0.]*len(components)
        triangles_min = [0.]*len(components)
        squares_av = [0.]*len(components)
        squares_max = [0.]*len(components)
        squares_min = [0.]*len(components)
        transitivity = [0.]*len(components)
        rc = [0.]*len(components)
        loopnumber = [0]*len(components)


        for compnum in np.arange(len(components)):
            # property 6: ischain?(remove self-loops and then test this property)
            # want: how many chains does the graph contain.. look at each component, not the whole graph in one go.
            # most graphs are single components.
            G1 = G.subgraph(list(components[compnum]))
            Gnoself = G1.copy()
            Gnoself.remove_edges_from(Gnoself.selfloop_edges())
            Unoself = nx.Graph()
            Unoself. add_edges_from(Gnoself.edges)
            
            # if all in and out degrees are 1, graph is a chain..do not include in trees
            indeg2 = []; 
            outdeg2 = []; 
            indeg_ls2 = list(Gnoself.in_degree())
            outdeg_ls2 = list(Gnoself.out_degree())
            # nx gives indeg and outdeg as tuples (nodename, in/out deg). which is why i need the for loop below
            for x in np.arange(len(G1.nodes())):
                indeg2.append(indeg_ls2[x][1])
                outdeg2.append(outdeg_ls2[x][1])
            indeg_by_comp.append(int_arr_to_str(indeg2, delim=';'))
            outdeg_by_comp.append(int_arr_to_str(outdeg2, delim=';'))

            indeg2 = np.array(indeg2)
            outdeg2 = np.array(outdeg2)
            in_min_out = indeg2 - outdeg2
            ischain[compnum] = int((np.sum(in_min_out) == 0) & (np.sum(np.abs(in_min_out)) == 2) & (np.all(indeg2<=1)) & (np.all(outdeg2<=1)))
            # property 7: istree(remove chains first)
            istree[compnum] = int((nx.is_tree(Gnoself) - ischain[compnum]) > 0)
            # property 8: isdag(only looking at DAGs other than trees and chains)
            isdag[compnum] = int((int(nx.is_directed_acyclic_graph(Gnoself)) - istree[compnum] - ischain[compnum]) > 0)
            if isdag[compnum] > 0:
                loopnumber[compnum] = len(list(Gnoself.edges)) - (len(list(Gnoself.nodes)) - 1)
            # property 9: single celled
            unicel[compnum] = int(len(Gnoself.nodes) == 1)
            istree[compnum] = int(istree[compnum]) - int(unicel[compnum]) # nx counts single node with no self-edge as a tree   
            # property 10: isscc (excluding unicellular)
            num_strong_comp2 = nx.number_strongly_connected_components(Gnoself)
            isscc[compnum] = int(num_strong_comp2 == 1)
            isscc[compnum] = int((isscc[compnum] - unicel[compnum]) > 0)
            # property 11: iscyc(cyclic graphs other than those with a single scc and single celled graphs) 
            iscyc[compnum] = int((isdag[compnum] + istree[compnum]+ ischain[compnum] + isscc[compnum] + unicel[compnum]) == 0)        
            # property 12: is eulerian
            iseul[compnum] = int(nx.is_eulerian(Gnoself))
            # property 13: node connectivity
            node_conn[compnum] = approx.node_connectivity(Gnoself)
            # property 14: clustering coefficient
            av_clust[compnum] = nx.average_clustering(Gnoself)
            # property 15: assortativity(pearson's coefficient)
            try:
                assort[compnum] = nx.degree_pearson_correlation_coefficient(Gnoself)#####################check
            except:
                assort[compnum] = 0.0
            # property 16,17,18: in degree centrality (average, maximum and minimum)
            indeg_cen = []
            dict1 = nx.in_degree_centrality(Gnoself)
            for a1 in dict1:
                indeg_cen.append(dict1[a1])            
            indeg_cen_av[compnum] = np.average(indeg_cen)
            indeg_cen_max[compnum] = max(indeg_cen)
            indeg_cen_min[compnum] = min(indeg_cen)
            # property 19,20,21: out degree centrality (average, maximum, minimum)
            outdeg_cen = []
            dict1 = nx.out_degree_centrality(Gnoself)
            for a1 in dict1:
                outdeg_cen.append(dict1[a1])             
            outdeg_cen_av[compnum] = np.average(outdeg_cen)
            outdeg_cen_max[compnum] = max(outdeg_cen)
            outdeg_cen_min[compnum] = min(outdeg_cen)
            # property 22,23,24: betweenness centrality (average,maximum, minimum)
            bet_cen = []
            dict1 = nx.betweenness_centrality(Gnoself)
            for a1 in dict1:
                bet_cen.append(dict1[a1])             
            bet_cen_av[compnum] = np.average(bet_cen)
            bet_cen_max[compnum] = max(bet_cen)
            bet_cen_min[compnum] = min(bet_cen)
            # property 25,26,27: eigen vector centrality (average,maximum, minimum)
            eig_cen = []
            try:
                dict1 = nx.eigenvector_centrality(Gnoself)
                for a1 in dict1:
                    eig_cen.append(dict1[a1])             
                eig_cen_av[compnum] = np.average(eig_cen)
                eig_cen_max[compnum] = max(eig_cen)
                eig_cen_min[compnum] = min(eig_cen)
            except nx.PowerIterationFailedConvergence:
                pass
            # property 28,29,30: number of triangles for each node (average,maximum, minimum)
            triangles = []
            dict1 = nx.triangles(Unoself)
            for a1 in dict1:
                triangles.append(dict1[a1])             
            if len(triangles):
                triangles_av[compnum] = np.average(triangles)
                triangles_max[compnum] = max(triangles)
                triangles_min[compnum] = min(triangles)
            # property 31: transitivity (fraction of all possible triangles present in the graph)
            transitivity[compnum] = nx.transitivity(Gnoself)
            # property 32,33,34: square clustering for each node(fraction of all possible squares present at a node)
            squares = []
            dict1 = nx.square_clustering(Gnoself)
            for a1 in dict1:
                squares.append(dict1[a1])             
            if len(squares):
                squares_av[compnum] = np.average(squares)
                squares_max[compnum] = max(squares)
                squares_min[compnum] = min(squares)
            # propery 35: rich club coefficient
            if len(list(Unoself.nodes())) > 3:
                rc[compnum] = 0.0
#               rc[compnum] = nx.rich_club_coefficient(Unoself).values()# only works if graph has 4 or more edges
            # property 36 and 37: number of source and target nodes

        iseul = sum(iseul)
        iscyc = sum(iscyc)
        isscc = sum(isscc)
        unicel = sum(unicel)
        isdag = sum(isdag)
        istree = sum(istree)
        ischain = sum(ischain)
        indeg_by_comp = ';'.join([str(x) for x in indeg_by_comp])
        outdeg_by_comp = ';'.join([str(x) for x in outdeg_by_comp])
        node_conn = ';'.join([str(x) for x in node_conn])# node connectivity for each component
        avav_clust = np.average(av_clust)# average clustering coefficient over all components
        av_clust = ';'.join([str(round(x,2)) for x in av_clust])# average clustering coefficients for each component
        av_assort = np.average(assort)# average assortativity over all components
        assort = ';'.join([str(round(x,2)) for x in assort])# assortativity for each component
        indeg_cen_avav = np.average(indeg_cen_av)# average indeg centrality over all components
        indeg_cen_av =';'.join([str(round(x,2)) for x in indeg_cen_av])# average indeg centrality for each component
        indeg_cen_maxmax = max(indeg_cen_max)# maximum indeg centrality across all components
        indeg_cen_max = ';'.join([str(round(x,2)) for x in indeg_cen_max])# maximum indeg centrality for each component
        indeg_cen_minmin = min(indeg_cen_min)# minimum indeg centrality across all components
        indeg_cen_min = ';'.join([str(round(x,2)) for x in indeg_cen_min])# minimum indeg centrality for each component 

        outdeg_cen_avav = np.average(outdeg_cen_av)
        outdeg_cen_av =';'.join([str(round(x,2)) for x in outdeg_cen_av])
        outdeg_cen_maxmax = max(outdeg_cen_max)
        outdeg_cen_max = ';'.join([str(round(x,2)) for x in outdeg_cen_max])
        outdeg_cen_minmin = min(outdeg_cen_min)
        outdeg_cen_min = ';'.join([str(round(x,2)) for x in outdeg_cen_min]) 
        bet_cen_avav = np.average(bet_cen_av)
        bet_cen_av =';'.join([str(round(x,2)) for x in bet_cen_av])
        bet_cen_maxmax = max(bet_cen_max)
        bet_cen_max = ';'.join([str(round(x,2)) for x in bet_cen_max])
        bet_cen_minmin = min(bet_cen_min)
        bet_cen_min = ';'.join([str(round(x,2)) for x in bet_cen_min]) 
        eig_cen_avav = np.average(eig_cen_av)
        eig_cen_av =';'.join([str(round(x,2)) for x in eig_cen_av])
        eig_cen_maxmax = max(eig_cen_max)
        eig_cen_max = ';'.join([str(round(x,2)) for x in eig_cen_max])
        eig_cen_minmin = min(eig_cen_min)
        eig_cen_min = ';'.join([str(round(x,2)) for x in eig_cen_min]) 
        triangles_avav = np.average(triangles_av)
        triangles_av =';'.join([str(x) for x in triangles_av])
        triangles_maxmax = max(triangles_max)
        triangles_max = ';'.join([str(x) for x in triangles_max])
        triangles_minmin = min(triangles_min)
        triangles_min = ';'.join([str(x) for x in triangles_min])
        transitivity_av = np.average(transitivity)
        transitivity_max = max(transitivity)
        transitivity_min = min(transitivity)
        transitivity = ';'.join([str(x) for x in transitivity]) 
        squares_avav = np.average(squares_av)
        squares_maxmax = max(squares_max)
        squares_minmin = min(squares_min)
        squares_av =';'.join([str(x) for x in squares_av])
        squares_max = ';'.join([str(x) for x in squares_max])
        squares_min = ';'.join([str(x) for x in squares_min])
        rc_av = np.average(rc)
        rc_max = max(rc)
        rc_min = min(rc) 
        rc = ';'.join([str(x) for x in rc])
        ln = [loopnumber[x] for x in np.nonzero(loopnumber)[0]]
        if any(ln):        
            loopnumber_av = np.average(ln)
        else:
            loopnumber_av = 0.0
        loopnumber = ';'.join([str(x) for x in loopnumber])

        # check.. sum of iscyc, isscc, unicel, dag,tree, chain should be the total number of components        
        if num_comp != (iscyc+isscc+unicel+isdag+istree+ischain):
            print('Number of components is wrong!!!!!!')
            print(num_comp)
            print([iscyc,isscc,unicel,isdag,istree,ischain])
            sys.exit()

        properties.append(indeg_by_comp)# string
        properties.append(outdeg_by_comp)#string
        properties.append(ischain)#int
        properties.append(istree)#int
        properties.append(isdag)#int
        properties.append(unicel)#int
        properties.append(isscc)#int
        properties.append(iscyc)#int
        properties.append(iseul)#int
        properties.append(loopnumber_av)#float
        properties.append(loopnumber)#string
        properties.append(node_conn)#string
        properties.append(avav_clust)#float
        properties.append(av_clust)#string
        properties.append(av_assort)#float
        properties.append(assort)#string
        properties.append(indeg_cen_avav)#float
        properties.append(indeg_cen_av)#string
        properties.append(indeg_cen_maxmax)#float
        properties.append(indeg_cen_max)#string
        properties.append(indeg_cen_minmin)#float
        properties.append(indeg_cen_min)#string
        properties.append(outdeg_cen_avav)#float
        properties.append(outdeg_cen_av)#string
        properties.append(outdeg_cen_maxmax)#float
        properties.append(outdeg_cen_max)#string
        properties.append(outdeg_cen_minmin)#float
        properties.append(outdeg_cen_min)#string
        properties.append(bet_cen_avav)#float
        properties.append(bet_cen_av)#string
        properties.append(bet_cen_maxmax)#float
        properties.append(bet_cen_max)#string
        properties.append(bet_cen_minmin)#float
        properties.append(bet_cen_min)#string
        properties.append(eig_cen_avav)#float
        properties.append(eig_cen_av)#string
        properties.append(eig_cen_maxmax)#float
        properties.append(eig_cen_max)#string
        properties.append(eig_cen_minmin)#float
        properties.append(eig_cen_min)#string
        properties.append(triangles_avav)#float
        properties.append(triangles_av)#string
        properties.append(triangles_maxmax)#float
        properties.append(triangles_max)#string
        properties.append(triangles_minmin)#float
        properties.append(triangles_min)#string
        properties.append(transitivity_av)# float
        properties.append(transitivity_max)#float
        properties.append(transitivity_min)#float
        properties.append(transitivity)#string
        properties.append(squares_avav)#float
        properties.append(squares_av)#string
        properties.append(squares_maxmax)#float
        properties.append(squares_max)#string
        properties.append(squares_minmin)#float
        properties.append(squares_min)#string
        properties.append(rc_av)# float
        properties.append(rc_max)#float
        properties.append(rc_min)#float
        properties.append(rc)#string



        
        # append more properties.....
				    # property 14: 
        
        # property x: in-degree sequence
        #indeg = # list(G.in_degree())[iterate over number of nodes][1]                                                                                                                   
        # property y: out-degree sequence
        #outdeg = # list(G.in_degree())[iterate over number of nodes][1]                                                                                                                       
        #.....
    else:
        properties = [0]*2 + [0.]*2 + [0] + ['']*2 + [0]*7 + [0.] + ['']*2 + [0.,'']*17 + [0.]*3 + [''] + [0.,'']*3 + [0.,0.,0.,'']
        
    # return list of properties
    return properties
    
def update_dataframe_with_graph_properties(df, random=False):
    pool = mp.Pool(N_PROC)
    names = ['num_components', 'num_strong_components','mean_degree','link_density','num_selfloop','indeg', 'outdeg', 'chain','tree','dag','unicell','scc','cyclic', 'is_eulerian','av_loopnum','loopnum','node_connectivity','clustering_av','clustering','assort_av','assort','indeg_cent_avav','indeg_cent_av','indeg_cent_maxmax','indeg_cent_max','indeg_cent_minmin','indeg_cent_min','outdeg_cent_avav','outdeg_cent_av','outdeg_cent_maxmax','outdeg_cent_max','outdeg_cent_minmin','outdeg_cent_min','bet_cent_avav','bet_cent_av','bet_cent_maxmax','bet_cent_max','bet_cent_minmin','bet_cent_min','eig_cent_avav','eig_cent_av','eig_cent_maxmax','eig_cent_max','eig_cent_minmin','eig_cent_min','triangles_avav','triangles_av','triangles_maxmax','triangles_max','triangles_minmin','triangles_min','transitivity_av','transitivity_max','transitivity_min','transitivity','squares_avav','squares_av','squares_maxmax','squares_max','squares_minmin','squares_min','rc_av','rc_max','rc_min','rc']
    if random:
            names = [f"ran_{s}" for s in names]
            properties = np.array(pool.map(get_graph_properties, df.ran_edges, CHUNK))
    else:
            properties = np.array(pool.map(get_graph_properties, df.edges, CHUNK))
 
    pool.close()

    types = [int]*2 + [float]*2 + [int] + [str]*2 + [int]*7 + [float] + [str]*2 + [float, str]*17 + [float]*3 + [str] + [float,str]*3 + [float]*3 + [str]

    for i in range(len(names)):
        df[names[i]] = properties[:,i].astype(types[i])

    return df

def load_dataframe():
    if os.path.exists(os.path.join(OUT_DIR, 'all_data.h5')):
        df = pd.read_hdf(os.path.join(OUT_DIR, 'all_data.h5'))
    else:
        df = amalgamate_csv_files_to_dataframe()
        df.to_hdf(os.path.join(OUT_DIR, 'all_data.h5'), 'df', table=True, mode='a')
    return df


      

def load_single_file(fName):
    cols = ['n_gene','n_fixed', 'n_osc', 'basin', 'ss', 'f_sig', 'm_sig', 'f_adj', 'm_adj', 'f_asy', 'm_asy', 'edges']
    ext = os.path.splitext(fName)[1]
    if ext == '.txt':
        return pd.read_csv(fName, names=cols)
    elif ext == '.csv':
        return pd.read_csv(fName)
    elif ext == '.hd5':
        return pd.read_hdf(fName)
    elif ext == '.pickle':
        return pd.read_pickle(fName)


def fix_column_types(df):
    for col in df.columns: 
        if df[col].dtype == object: 
            print('Function not working yet...')
            continue
            try: 
                if '.' in df.loc[0,col]:
                    df[col] = df[col].astype(float) 
                else:
                    i1 = int(df.loc[0,col]) 
                    df[col] = df[col].astype(int) 
            except: 
                pass 
    return df

def get_max_of_string_of_ints(st, delim=';'):
    return max([int(x) for x in st.split(delim)])

def process_database(fName):
#   accepted_directories = [f"n{n}" for n in range(3,6)]
#   for root, dirs, files in os.walk(data_dir):
#       if sum([1 for acc in accepted_directories if acc in root]):
#           for fName in files:
#               if os.path.splitext(fName)[1] != '.txt':
#                   continue

                timeS = time.time()

                # Load csv file
#               df = load_single_file(os.path.join(root, fName))
                df = load_single_file(fName)
                print(f"Time taken to load: {(time.time()-timeS)/60.} mins")
                # Remove results that did not converge
                df = df.drop(index=df.loc[df.edges.str.len()==1].index).reset_index(drop=True)
                df = add_column_parallel(df, 'edges', 'n_nodes', fn_n_nodes)
                df = add_column_parallel(df, 'edges', 'n_edges', fn_n_edges)
                print(f"Time taken to add_columns: {(time.time()-timeS)/60.} mins")
                # Get graph properties
                df = update_dataframe_with_graph_properties(df)
                print(f"Time taken to get_graph_props: {(time.time()-timeS)/60.} mins")

                ### Set up POOL
                pool = mp.Pool(N_PROC)


                ### Close POOL
                pool.close()

                
                # Save in pickle format
                new_name = os.path.splitext(fName)[0] + '.pickle'
                df.to_pickle(new_name)

                os.remove(fName)

    

if __name__ == "__main__":

    timeS = time.time()

    process_database(ARGS.fName)

    print(f"Time to run...  {(time.time()-timeS)/60.} minutes")





