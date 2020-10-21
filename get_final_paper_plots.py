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

N_PROC=55
CHUNK=100


def get_graph_properties(edges):
    # Set up graph
    connections = np.array([int(x) for x in edges.split(';')])
    num_edges = 0
    nodes = sorted(list(set(connections)))
    # Calculate Properties
    properties = []

    if connections[0] > 0:
        edges = connections.reshape(int(connections.size/2), 2)
        num_edges = int(connections.size/2)
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

        # find all components
        components = list(nx.connected_components(U));

        ischain = [None]*len(components)
        istree = [None]*len(components)
        isdag = [None]*len(components)
        unicel = [None]*len(components)
        isscc = [None]*len(components)
        iscyc = [None]*len(components)
        for compnum in np.arange(len(components)):
            # property 6: ischain?(remove self-loops and then test this property)
            # want: how many chains does the graph contain.. look at each component, not the whole graph in one go.
            # most graphs are single components.
            G1 = G.subgraph(list(components[compnum]))
            Gnoself = G1.copy()
            Gnoself.remove_edges_from(nx.selfloop_edges(Gnoself))
            Unoself = nx.Graph()
            Unoself. add_edges_from(Gnoself.edges)
            # if all in and out degrees are 1, graph is a chain..do not include in trees
            indeg2 = [];
            outdeg2 = [];
            indeg_ls2 = list(Gnoself.in_degree())
            outdeg_ls2 = list(Gnoself.out_degree())
            for x in np.arange(len(G1.nodes())):
                indeg2.append(indeg_ls2[x][1])
                outdeg2.append(outdeg_ls2[x][1])
            indeg2 = np.array(indeg2)
            outdeg2 = np.array(outdeg2)
            in_min_out = indeg2 - outdeg2
            ischain[compnum] = int((np.sum(in_min_out) == 0) & (np.sum(np.abs(in_min_out)) == 2) & (np.all(indeg2<=1)) & (np.all(outdeg2<=1)))
            # property 7: istree(remove chains first)
            istree[compnum] = int((nx.is_tree(Gnoself) - ischain[compnum]) > 0)
            # property 8: isdag(only looking at DAGs other than trees and chains)
            isdag[compnum] = int((int(nx.is_directed_acyclic_graph(Gnoself)) - istree[compnum] - ischain[compnum]) > 0)
            # property 9: single celled
            unicel[compnum] = int(len(Gnoself.nodes) == 1)
            istree[compnum] = int(istree[compnum]) - int(unicel[compnum]) # nx counts single node with no self-edge as a tree   
            # property 10: isscc (excluding unicellular)
            num_strong_comp2 = nx.number_strongly_connected_components(Gnoself)
            isscc[compnum] = int(num_strong_comp2 == 1)
            isscc[compnum] = int((isscc[compnum] - unicel[compnum]) > 0)
            # property 11: iscyc(cyclic graphs other than those with a single scc and single celled graphs) 
            iscyc[compnum] = int((isdag[compnum] + istree[compnum]+ ischain[compnum] + isscc[compnum] + unicel[compnum]) == 0)
        iscyc = sum(iscyc)
        isscc = sum(isscc)
        unicel = sum(unicel)
        isdag = sum(isdag)
        istree = sum(istree)
        ischain = sum(ischain)
    return iscyc, isscc, unicel, isdag, istree, ischain, num_comp, num_edges, indeg2,outdeg2

def reset_idx_Ind(inputs):
    [filenum,N,genomeID,initID,sig,adj,asym] = inputs
    idx_ind = Ind.index[(Ind.filenum == filenum)&(Ind.N == N)&(Ind.genomeID == genomeID)&(Ind.initID == initID)&(Ind.sig == sig)&(Ind.asym == asym)&(Ind.adj == adj)].tolist()
    return idx_ind


def topology_regeneration_plots():

    REP_DIR = "/home/somya/Somya/version1/reproducibility_data/"
    files3 = sorted(glob.glob(REP_DIR + 'reproducibility_N*_*_genome*.txt'))#[:1]
    IND_DIR = "/home/somya/Somya/version1/independance_data/"
    files4 = sorted(glob.glob(IND_DIR + 'cell_categories_N*_*_genome*.txt'))#[:1]
    INT_DIR = "/home/somya/Somya/version1/intrinsic_independance/"
    files5 = sorted(glob.glob(INT_DIR + 'intrinsically_independant_cells_N*_*_genome*.txt'))
    PLU_DIR = "/home/somya/Somya/version1/reproducibility_data/pluripotent_cell_identity/"    
    files6 = sorted(glob.glob(PLU_DIR + 'pluripotent_N*_*_genome*.txt'))
    BAS_DIR = "/home/somya/Somya/version1/basin_sizes/"
    files7 = sorted(glob.glob(BAS_DIR + 'basin_sizes_N*_*_genome*.txt'))

    Rep = pd.read_csv(files3[0],header = None)
    Rep.columns = ["N","genomeID","initID","sig","asym","adj","graph_nodes","frac_rep_graph","nongraph_nodes","frac_rep_nongraph",'graph']
    Rep['filenum'] = [0]*len(Rep) 

    Ind = pd.read_csv(files4[0],header = None)
    Ind.columns = ["N","genomeID","initID","sig","asym","adj","ab","Ab","ABCDE","ABcDE",'ABcdE','ABCDe','ABcDe','ABcde','aBCDE','aBcDE','aBcdE','aBCDe','aBcDe','aBcde']
    Ind['filenum'] = [0]*len(Ind)

    Int = pd.read_csv(files5[0],header = None)
    Int.columns = ["N","genomeID","initID","sig","asym","adj","CDE","CDe",'cDE','cDe','cdE','cde']
    Int['filenum'] = [0]*len(Int)
    
    Plu = pd.read_csv(files6[0],header = None)
    Plu.columns = ["N","genomeID","initID","sig","asym","adj","frac_rep_graph",'traj_lengths','pluripotent_cells',"nongraph_nodes","frac_rep_nongraph",'graph']
    Plu['filenum'] = [0]*len(Plu)

    Bas = pd.read_csv(files7[0],header = None)
    Bas.columns = ["N","genomeID","initID","sig","asym","adj","basins"]
    #Bas['filenum'] = [0]*len(Bas)
    Rep['basin_sizes'] = list(set(Bas['basins']))*len(Rep)


     
    for i in np.linspace(0,len(files3)-2,len(files3)-1):
        Rep1 = pd.read_csv(files3[int(i)+1],header = None)
        Rep1.columns = ["N","genomeID","initID","sig","asym","adj","graph_nodes","frac_rep_graph","nongraph_nodes","frac_rep_nongraph",'graph']
        Rep1['filenum'] = [i+1]*len(Rep1)
        
        Bas1 = pd.read_csv(files7[int(i)+1],header = None)
        Bas1.columns = ["N","genomeID","initID","sig","asym","adj","basins"]
        Rep1['basin_sizes'] = list(set(Bas1['basins']))*len(Rep1)
        
        Rep = pd.concat([Rep,Rep1], ignore_index=True)

        Ind1 = pd.read_csv(files4[int(i)+1],header = None)
        Ind1.columns = ["N","genomeID","initID","sig","asym","adj","ab","Ab","ABCDE","ABcDE",'ABcdE','ABCDe','ABcDe','ABcde','aBCDE','aBcDE','aBcdE','aBCDe','aBcDe','aBcde']
        Ind1['filenum'] = [i+1]*len(Ind1)
        Ind = pd.concat([Ind,Ind1], ignore_index=True)

        Int1 = pd.read_csv(files5[int(i)+1],header = None)
        Int1.columns = ["N","genomeID","initID","sig","asym","adj","CDE","CDe",'cDE','cDe','cdE','cde']
        Int1['filenum'] = [i+1]*len(Int1)
        Int = pd.concat([Int,Int1], ignore_index=True)

        Plu1 = pd.read_csv(files6[int(i)+1],header = None)
        Plu1.columns = ["N","genomeID","initID","sig","asym","adj","frac_rep_graph",'traj_lengths','pluripotent_cells',"nongraph_nodes","frac_rep_nongraph",'graph']
        Plu1['filenum'] = [i+1]*len(Plu1)
        Plu = pd.concat([Plu,Plu1], ignore_index=True)
        

        print(i+1)

    Rep = Rep.drop_duplicates()
    Ind = Ind.drop_duplicates()
    Int = Int.drop_duplicates()
    Plu = Plu.drop_duplicates()
    
    
    Rep = Rep.sort_values(['filenum','N','genomeID','initID','sig','asym','adj'],ascending=[True,True,True,True,True,True,True])
    Ind = Ind.sort_values(['filenum','N','genomeID','initID','sig','asym','adj'],ascending=[True,True,True,True,True,True,True])
    Int = Int.sort_values(['filenum','N','genomeID','initID','sig','asym','adj'],ascending=[True,True,True,True,True,True,True])
    Plu = Plu.sort_values(['filenum','N','genomeID','initID','sig','asym','adj'],ascending=[True,True,True,True,True,True,True])    
   #Bas = Bas.sort_values(['filenum','N','genomeID','initID','sig','asym','adj'],ascending=[True,True,True,True,True,True,True])    

    Rep = Rep.reset_index(drop=True)
    Ind = Ind.reset_index(drop=True)
    Int = Int.reset_index(drop=True)
    Plu = Plu.reset_index(drop=True)
   #Bas = Bas.reset_index(drop=True)

    # get graph properties: number of components, graph type (unicell, scc, cyclic, dag, tree, chain)
    scc = [0]*len(Rep)
    chain = [0]*len(Rep)
    tree = [0]*len(Rep)
    dag = [0]*len(Rep)
    cyc = [0]*len(Rep)
    unicel = [0]*len(Rep)
    numcomp = [0]*len(Rep)
    edges = [0]*len(Rep)
    indeg2 = [0]*len(Rep)
    outdeg2 = [0]*len(Rep)
    with mp.Pool(N_PROC) as pool:
#       results = np.array(pool.map(get_graph_properties,Rep.graph,CHUNK))
#       cyc,scc, unicel, dag, tree, chain, numcomp, edges = results.T

        results = pool.map(get_graph_properties,Rep.graph,CHUNK)
        for i in range(len(Rep)):
            cyc[i], scc[i], unicel[i], dag[i], tree[i], chain[i], numcomp[i], edges[i], indeg2[i], outdeg2[i] = results[i]
    Rep['scc'] = scc
    Rep['tree'] = tree
    Rep['cyc'] = cyc
    Rep['dag'] = dag
    Rep['chain'] = chain
    Rep['unicell'] = unicel
    Rep['num_components'] = numcomp
    Rep['num_edges'] = edges
    # remove rows corresponding to num_components>1 from both dat and rep
    a1 = np.where(Rep['num_components'] > 1)
    Rep = Rep.drop(Rep.index[list(a1)])
    Ind = Ind.drop(Ind.index[list(a1)])
    Int = Int.drop(Int.index[list(a1)])
    Plu = Plu.drop(Plu.index[list(a1)])

    indeg3 = [indeg2[i] for i in list(np.where(np.array(numcomp)==1)[0])]
    outdeg3 = [outdeg2[i] for i in list(np.where(np.array(numcomp)==1)[0])]
    Rep['pluripotent_cells'] = Plu['pluripotent_cells']
    Rep['traj_lengths'] = Plu['traj_lengths']
    return Rep, Ind, Int, indeg3, outdeg3

def get_list_from_string(st,delim=';'):
    a = np.array([int(x) for x in st.split(delim)])
    return a

def get_expanded_list(x):
    x1 = len(x)
    return [x1]*x1

def numbasin_basinsize_scatter(Rep):
    with mp.Pool(N_PROC) as pool:
        bas_size_list = np.array(pool.map(get_list_from_string, Rep.basin_sizes, CHUNK))
        x1 = np.array(pool.map( get_expanded_list, bas_size_list, CHUNK))
    print('got number of basins, basin_size lists')
    numbasins = np.concatenate(x1).ravel()
    bas_size_list = np.concatenate(bas_size_list).ravel()
    return numbasins, bas_size_list

def get_basinsize_list_nwise(N,Rep):
    #returns a list of all basin sizes pooled from across all systems with a given N
    bas_n = Rep.basin_sizes[Rep.N==N].reset_index(drop=True)
    with mp.Pool(N_PROC) as pool:
        bas_n_list = np.array(pool.map(get_list_from_string,bas_n,CHUNK))
    bas_av_list = [np.average(x) for x in bas_n_list]
    bas_n_list = np.concatenate(bas_n_list).ravel()
    return bas_n_list, bas_av_list

def get_Nwise_basinsize_hist(Rep):
    bas = [get_basinsize_list_nwise(i,Rep) for i in [3,4,5,6,7]]
    basN = [bas[i][0] for i in range(len(bas))]
    bas_av = [bas[i][1] for i in range(len(bas))]

    print('got basN, bas_av') 
    m1 = max([max(basN[i]) for i in range(len(basN))])
    a1 = [0]*len(basN); a2 = [0]*len(basN)
    for i1 in range(len(basN)):
        a1[i1],a2[i1] = np.histogram(basN[i1], bins=[x+1 for x in range(m1+1)])
        a1[i1] = [x/sum(a1[i1]) for x in a1[i1]]
    print('got histogram N')    

    m2 = max([max(bas_av[i]) for i in range(len(bas_av))])
    b1 = [0]*len(bas_av); b2 = [0]*len(bas_av)
    for i1 in range(len(bas_av)):
        b1[i1],b2[i1] = np.histogram(bas_av[i1], bins=np.linspace(1,m2,10))#minimum basin size is 1
        b1[i1] = [x/sum(b1[i1]) for x in b1[i1]]
    print('got histogram av')

    #plot histogram
    color=['xkcd:medium green','xkcd:denim blue','xkcd:gold','xkcd:pale red','grey']
    plt.figure();plt.bar(a2[0][0:-1]-0.3, a1[0], width=0.1, color=color[0])
    plt.bar(a2[1][0:-1]-0.15, a1[1], width=0.1, color=color[1])
    plt.bar(a2[2][0:-1], a1[2], width=0.1, color=color[2])
    plt.bar(a2[3][0:-1]+0.15, a1[3], width=0.1, color=color[3])
    plt.bar(a2[4][0:-1]+0.3, a1[4], width=0.1, color=color[4])
    plt.legend(['N=3','4','5','6','7'])

    plt.figure();plt.bar(b2[0][0:-1]-0.06, b1[0], width=0.02, color=color[0])
    plt.bar(b2[1][0:-1]-0.03, b1[1], width=0.02, color=color[1])
    plt.bar(b2[2][0:-1], b1[2], width=0.02, color=color[2])
    plt.bar(b2[3][0:-1]+0.03, b1[3], width=0.02, color=color[3])
    plt.bar(b2[4][0:-1]+0.06, b1[4], width=0.02, color=color[4])
    plt.legend(['N=3','4','5','6','7'])
    return a1, a2, b1, b2


def get_basinsize_list_acyclic(Rep):
    #returns a list of all basin sizes pooled from across all systems with a given N

    bas_tree = Rep.basin_sizes[Rep.tree==1].reset_index(drop=True)
    bas_dag = Rep.basin_sizes[Rep.dag==1].reset_index(drop=True)
    bas_chain = Rep.basin_sizes[Rep.chain==1].reset_index(drop=True)
    bas_cyclic = Rep.basin_sizes[(Rep.cyc)==1].reset_index(drop=True)
    bas_scc = Rep.basin_sizes[Rep.scc==1].reset_index(drop=True)
    bas_unicellular = Rep.basin_sizes[(Rep.cyc+Rep.scc+Rep.tree+Rep.dag+Rep.chain)==0].reset_index(drop=True)
    
    with mp.Pool(N_PROC) as pool:
        bas_tree_list = np.array(pool.map(get_list_from_string,bas_tree,CHUNK))
        bas_dag_list = np.array(pool.map(get_list_from_string,bas_dag,CHUNK))
        bas_chain_list = np.array(pool.map(get_list_from_string,bas_chain,CHUNK))
        bas_cyclic_list = np.array(pool.map(get_list_from_string,bas_cyclic,CHUNK))
        bas_scc_list = np.array(pool.map(get_list_from_string,bas_scc,CHUNK))
        bas_unicellular_list = np.array(pool.map(get_list_from_string,bas_unicellular,CHUNK))

    bas_av_tree_list = np.array([np.average(x) for x in bas_tree_list])
    bas_tree_list = np.concatenate(bas_tree_list).ravel()

    bas_av_dag_list = np.array([np.average(x) for x in bas_dag_list])
    bas_dag_list = np.concatenate(bas_dag_list).ravel()

    bas_av_chain_list = np.array([np.average(x) for x in bas_chain_list])
    bas_chain_list = np.concatenate(bas_chain_list).ravel()

    bas_av_cyclic_list = np.array([np.average(x) for x in bas_cyclic_list])
    bas_cyclic_list = np.concatenate(bas_cyclic_list).ravel()

    bas_av_scc_list = np.array([np.average(x) for x in bas_scc_list])
    bas_scc_list = np.concatenate(bas_scc_list).ravel()

    bas_av_unicellular_list = np.array([np.average(x) for x in bas_unicellular_list])
    bas_unicellular_list = np.concatenate(bas_unicellular_list).ravel()
    return bas_unicellular_list, bas_cyclic_list, bas_chain_list, bas_scc_list, bas_dag_list, bas_tree_list, bas_av_unicellular_list, bas_av_cyclic_list, bas_av_chain_list, bas_av_scc_list, bas_av_dag_list, bas_av_tree_list 


def get_acyclic_basinsize_hist(Rep):
    bas_unicellular_list, bas_cyclic_list, bas_chain_list, bas_scc_list, bas_dag_list, bas_tree_list, bas_av_unicellular_list, bas_av_cyclic_list, bas_av_chain_list, bas_av_scc_list, bas_av_dag_list, bas_av_tree_list= get_basinsize_list_acyclic(Rep)
    basN = [bas_unicellular_list, bas_scc_list, bas_cyclic_list, bas_chain_list, bas_tree_list, bas_dag_list]
    bas_av = [bas_av_unicellular_list, bas_av_scc_list, bas_av_cyclic_list, bas_av_chain_list, bas_av_tree_list, bas_av_dag_list]
    print('got basN, bas_av')
    
    basN = [bas_unicellular_list, np.concatenate((bas_scc_list, bas_cyclic_list)), bas_chain_list, bas_tree_list, bas_dag_list] 
    bas_av = [bas_av_unicellular_list, np.concatenate((bas_av_scc_list, bas_av_cyclic_list)), bas_av_chain_list, bas_av_tree_list, bas_av_dag_list]
    m1 = max([max(basN[i]) for i in range(len(basN))])
    a1 = [0]*len(basN); a2 = [0]*len(basN)
    for i1 in range(len(basN)):
        a1[i1],a2[i1] = np.histogram(basN[i1], bins=[x+1 for x in range(m1+1)])
        a1[i1] = [x/sum(a1[i1]) for x in a1[i1]]
    print('got histogram N')    

    m2 = max([max(bas_av[i]) for i in range(len(bas_av))])
    b1 = [0]*len(bas_av); b2 = [0]*len(bas_av)
    
    for i1 in range(len(bas_av)):
        b1[i1],b2[i1] = np.histogram(bas_av[i1], bins=np.linspace(1,m2,10))# minimum basin size is 1
        b1[i1] = [x/sum(b1[i1]) for x in b1[i1]]
    print('got histogram av')

    # acyclic vs cyclic/unicellular
    c1 = [0]*2
    c1[0] = [b1[0][i] + b1[1][i] for i in range(len(b1[0]))]# cyclic/scc/unicellular
    c1[0] = [x/sum(c1[0]) for x in c1[0]]
    c1[1] = [b1[2][i] + b1[3][i] + b1[4][i] for i in range(len(b1[2]))]#chain/tree/dag
    c1[1] = [x/sum(c1[1]) for x in c1[1]]

    #plot histogram
    color=['xkcd:gunmetal','xkcd:light grey','xkcd:pale red','k','xkcd:sky blue']
    
   # color = list(np.array(Set1_9.mpl_colors)[[6, 3, 2, 0, 1]]) + ['k']
    plt.figure();#plt.bar(a2[0][0:-1]-0.3, a1[0], width=0.1, color=color[0])
    plt.bar(a2[0][0:-1]-0.15, a1[0], width=0.1, color=color[0])
    plt.bar(a2[1][0:-1], a1[0], width=0.1, color=color[1])
    plt.bar(a2[2][0:-1]+0.15, a1[2], width=0.1, color=color[2])
    plt.bar(a2[3][0:-1]+0.3, a1[3], width=0.1, color=color[3])
    plt.bar(a2[4][0:-1]+0.35, a1[4], width=0.1, color=color[4])
    plt.legend(['unicellular','cyclic','chain','tree','DAG'])
    
    plt.figure();#plt.bar(b2[0][0:-1]-0.0625, b1[0], width=0.015, color=color[0])
    plt.bar(b2[0][0:-1]-0.0375, b1[0], width=0.015, color=color[0])
    plt.bar(b2[1][0:-1]-0.0125, b1[1], width=0.015, color=color[1])
    plt.bar(b2[2][0:-1]+0.0125, b1[2], width=0.015, color=color[2])
    plt.bar(b2[3][0:-1]+0.0375, b1[3], width=0.015, color=color[3])
    plt.bar(b2[4][0:-1]+0.0625, b1[4], width=0.015, color=color[4])
    plt.legend(['unicellular','cyclic','chain','tree','DAG'])
    
    plt.plot(b2[0][0:-1],c1[0], color = 'xkcd:light grey')
    plt.plot(b2[0][0:-1],c1[1], color = 'k')

    return a1, a2, b1, b2

def get_true_root_nodes(graph):
    # all other nodes in the graph are reachable from true root nodes
    edges = np.array([int(x) for x in graph.split(';')])# edge list for the homeostatic graph
    edges = edges.reshape(int(edges.size/2), 2)
    G = nx.DiGraph()
    G.add_edges_from(edges)
    G.remove_edges_from(nx.selfloop_edges(G))
    nodes = list(set(G.nodes))
    true_root_nodes = np.zeros((len(nodes)))
    for i in range(len(nodes)):
        x1 = np.array([nx.has_path(G,source = nodes[i], target = nodes[j]) for j in range(len(nodes)) if j!=i])
        true_root_nodes[i] = not np.any(x1==0)
    return true_root_nodes

def trajlength_vs_graph_pathlengths(Rep):
    Rep_acyc = Rep[(Rep.chain>0)|(Rep.tree>0)|(Rep.dag>0)].reset_index(drop=True)
    
    Traj_len = np.array(())# regeneration trajectory lengths form all plu_root_nodes 
    Graph_pathlen = np.array(())# max path length from root node to any other node in lineage graph
    True_root_nodes = np.array(())# number of root nodes in graph
    Plu_nodes_num = np.array(())# number of pluripotent nodes in graph
    Plu_root_nodes = np.array(())# number of pluripotent nodes that are root nodes    

    for i in range(len(Rep_acyc)):
        t1 = np.array([int(x)-1 for x in Rep_acyc.traj_lengths[i].split(';')])# trajectory lengths to steady state for all cells in organism
        p1 = np.array([int(x) for x in Rep_acyc.pluripotent_cells[i].split(';')])# identity of all pluripotent cells
        t1 = t1[p1>0]# trajectory lengths from pluripotent cells to the same organism
        
        edges = np.array([int(x) for x in Rep_acyc.graph[i].split(';')])# edge list for the homeostatic graph
        edges = edges.reshape(int(edges.size/2), 2)
        G = nx.DiGraph()
        G.add_edges_from(edges)
        G.remove_edges_from(nx.selfloop_edges(G))
        nodes = set(G.nodes)
    
        plu_nodes = np.array(list(nodes))[p1>0]
        Plu_nodes_num = np.append(Plu_nodes_num,len(plu_nodes))

        receivers = set([x[1] for x in list(G.edges)])
        donors = set([x[0] for x in list(G.edges)])
        roots = list(nodes - receivers)
        leaf = list(nodes - donors)
        true_root_nodes = get_true_root_nodes(Rep_acyc.graph[i])# boolean array
        true_root_nodes = np.array(list(nodes))[true_root_nodes>0]
        True_root_nodes = np.append(True_root_nodes,len(true_root_nodes))

        plu_root_nodes = np.array(list(set(plu_nodes).intersection( set(true_root_nodes))))# all the pluripotent root nodes
        Plu_root_nodes = np.append(Plu_root_nodes,len(plu_root_nodes))
        if any(plu_root_nodes):
            traj_plu_root = [t1[x] for x in range(len(t1)) if plu_nodes[x] in set(plu_root_nodes) ]# traj length from all pluripotent root nodes
            Traj_len = np.append(Traj_len,traj_plu_root)

            max_graph_path = np.zeros((len(plu_root_nodes)))
            for j in range(len(plu_root_nodes)):
                max_graph_path[j] = max([nx.shortest_path_length(G, source=plu_root_nodes[j], target=leaf[k]) for k in range(len(leaf))])
            Graph_pathlen = np.append(Graph_pathlen,max_graph_path)           
        print(i)
    return Traj_len, Graph_pathlen, Plu_nodes_num, True_root_nodes, Plu_root_nodes

def get_regeneration_traj_plots(Rep):
# plot1: rootnodes vs pluripotent nodes

    x1_chain, x2_chain, x3_chain, x4_chain, x5_chain = trajlength_vs_graph_pathlengths(Rep[Rep.chain>0])
    x1_tree, x2_tree, x3_tree, x4_tree, x5_tree = trajlength_vs_graph_pathlengths(Rep[Rep.tree>0])
    x1_dag, x2_dag, x3_dag, x4_dag, x5_dag = trajlength_vs_graph_pathlengths(Rep[Rep.dag>0])

   # cell_cat = [not reg- not root, root | reg- not root, root]

    graph_top = ['chain','tree','dag'] 
    X = range(1,4) 
    Y = np.zeros((3,4)) 

    #chains
    Y[0,0] = np.sum(Rep.graph_nodes[Rep.chain>0] - x4_chain - x3_chain + x5_chain)
    Y[0,1] = np.sum(x4_chain - x5_chain)
    Y[0,2] = np.sum(x3_chain - x5_chain)
    Y[0,3] = np.sum(x5_chain)

    #trees
    Y[1,0] = np.sum(Rep.graph_nodes[Rep.tree>0] - x4_tree - x3_tree + x5_tree)
    Y[1,1] = np.sum(x4_tree - x5_tree)
    Y[1,2] = np.sum(x3_tree - x5_tree)
    Y[1,3] = np.sum(x5_tree)

    #dags
    Y[2,0] = np.sum(Rep.graph_nodes[Rep.dag>0] - x4_dag - x3_dag + x5_dag)
    Y[2,1] = np.sum(x4_dag - x5_dag)
    Y[2,2] = np.sum(x3_dag - x5_dag)
    Y[2,3] = np.sum(x5_dag)

    for i in range(3): 
       #print(Y[:,i].sum())
        Y[i] = Y[i] / Y[i].sum() 

    width = 0.75    
    col = ['xkcd:medium grey','xkcd:light grey','xkcd:brown','xkcd:puce']
    fig, ax = plt.subplots() 
    bottom = np.zeros(3) 
    for i in range(4): 
        ax.bar(X, Y[:,i], width, bottom=bottom, color=col[i], edgecolor='white', linewidth=0.5) 
        bottom += Y[:,i] 
    plt.xticks([1,2,3],('chain','tree','DAG')) 
    plt.xlim((0.5,3.5)) 


# plot2: regeneration trajectory length vs graph path length

#plt.figure();jitterplot.jitter_normal(x2_chain,x1_chain,max_noise=0.1,c='r',s=1)
#plt.figure();jitterplot.jitter_normal(x2_dag,x1_dag,max_noise=0.1,c='k',s=1)
#plt.figure();jitterplot.jitter_normal(x2_tree,x1_tree,max_noise=0.15,c='b',s=1)
    plt.figure();
    x1 = x1_chain-x2_chain
    x2 = x1_tree-x2_tree
    x3 = x1_dag-x2_dag
    #bins = np.linspace(0,50,11)
    
    f1 = plt.hist([x1,x2,x3],log=True,density=True,color=['red','blue','black'])
    #plt.xticks(bins)
    #ylim = np.sum(f1[0][0])
    #plt.yticks(np.linspace(0,ylim,3),('0','0.5','1'))
    return 0

    
def count_links_chains_in_acyclic_graph(edges):
    edges = np.array([int(x) for x in edges.split(';')])
    edges = edges.reshape(int(edges.size/2), 2) 
    G = nx.DiGraph()
    G.add_edges_from(edges)
    G.remove_edges_from(nx.selfloop_edges(G))
    nodes = set(G.nodes)
    receivers = set([x[1] for x in list(G.edges)])
    donors = set([x[0] for x in list(G.edges)])
    roots = list(nodes - receivers) 
    leaves = list(nodes - donors)
    x = [len(list(nx.all_simple_paths(G,source=i,target=j)))for i in roots for j in leaves]
    num_chains_in_graph = sum(x) - 1 #since chains have 0 branches
    num_links = len(G.edges) - len(nodes) + 1;
    return num_links, num_chains_in_graph


def acyclize_cyclic_graph(edges):
    edges = np.array([int(x) for x in edges.split(';')])
    edges = edges.reshape(int(edges.size/2),2)
    G = nx.DiGraph()
    G.add_edges_from(edges)
    G.remove_edges_from(nx.selfloop_edges(G))
    edges2 = np.array(list(G.edges))
    # merge cyclic edges of the graph
    # find all sccs in the graph
    scc_list = list(nx.strongly_connected_components(G))# each element of the list is a set of nodes in that scc (singleton set = single node not part of any cycle)
    node_scc_map = np.zeros((max(G.nodes)+1))
    edge_donors = np.zeros((len(edges2))); edge_receivers = np.zeros((len(edges2)))
    for i in range(len(scc_list)):
        node_scc_map[list(scc_list[i])] = i
    for i in range(len(edges2)):
        edge_donors[i] = node_scc_map[edges2[i][0]]
        edge_receivers[i] = node_scc_map[edges2[i][1]]
    edges_new = np.array([edge_donors,edge_receivers]).astype(int).T
    # make a new graph with all sccs as nodes, and edges connecting any node in scc-i with any node in scc-j as a an edge between new nodes i and j
    G2 = nx.DiGraph()
    G2.add_edges_from(edges_new)
    G2.remove_edges_from(nx.selfloop_edges(G2))
    #print(np.array(G2.edges))
    edges_new = np.array(list(G2.edges))
    num_links = 0; num_chains = -1
    if edges_new.size>0:
        # count links and chains in the new graph
        edges_new_str = ';'.join([str(x) for x in edges_new.reshape(len(edges_new)*2)])
        num_links, num_chains = count_links_chains_in_acyclic_graph(edges_new_str)
    # return (total number of sccs)/ (total number of nodes) in the original graph as a measure of 'acyclicity'
    return len(scc_list)/len(G.nodes), len(edges_new)/len(edges2), num_links, num_chains 


def get_branches_links_scatter(Rep, acyc):    
    x2 = (Rep.graph_nodes*Rep.frac_rep_graph + Rep.nongraph_nodes*Rep.frac_rep_nongraph)/(Rep.graph_nodes + Rep.nongraph_nodes)
    Rep['reg_cap'] = Rep['frac_rep_graph']/x2#regeneration capacity
    if acyc:
        acyc_graphs = Rep.graph[(Rep.tree>0)|(Rep.chain>0)|(Rep.dag>0)].reset_index(drop=True)
        reg_cap = Rep.reg_cap[(Rep.tree>0)|(Rep.chain>0)|(Rep.dag>0)].reset_index(drop=True)
        num_links = []
        num_chains = []
        node_frac = [1]*len(acyc_graphs)# the fraction of nodes left in the graph after acyclization. Does not affect acyclic graphs
        edge_frac = [1]*len(acyc_graphs)# fraction of edges left after acyclization
        for i in acyc_graphs:
            a1, a2 = count_links_chains_in_acyclic_graph(i)
            num_links.append(a1) 
            num_chains.append(a2)
    else:
        cyc_graphs = Rep.graph[(Rep.cyc>0)|(Rep.scc>0)].reset_index(drop=True)
        reg_cap = Rep.reg_cap[(Rep.cyc>0)|(Rep.scc>0)].reset_index(drop=True)
        num_links = []
        num_chains = []
        node_frac = []
        edge_frac = []
        for i in cyc_graphs:
            print(i)
            a1, a2, a3, a4 = acyclize_cyclic_graph(i)
            num_links.append(a3)
            num_chains.append(a4)
            node_frac.append(a1)
            edge_frac.append(a2)
        
    num_chains = np.array(num_chains);num_links = np.array(num_links); node_frac = np.array(node_frac); edge_frac = np.array(edge_frac)
    # zoomed in scatter plots
    # if reg==0:
       # c='xkcd:dark grey'
    # else:
    #    c = ['xkcd:dark grey' if x<=1 else 'xkcd:puce' for x in acyc_reg[(num_chains<=20)&(num_links<=20)]]
    c='xkcd:medium grey'
    plt.figure();jitterplot.jitter_normal(num_chains[(num_chains<=20)&(num_links<=20)&(edge_frac>0)],num_links[(num_chains<=20)&(num_links<=20)&(edge_frac>0)],c=c, max_noise=0.1);

    # zoomed out scatter plots:
    plt.figure();jitterplot.jitter_normal(num_chains[edge_frac>0],num_links[edge_frac>0],c='xkcd:medium grey', max_noise=0.1);

    return num_links, num_chains, np.array(reg_cap), node_frac, edge_frac


def frac_relaxed_trees(num_links, num_chains, edge_frac, num_all_graphs):
    x = np.linspace(0,1,11)# threshold for nl/nb to be considered as trees
    y = np.linspace(0,1,11)# threshold for edge_frac to be considered acyclic)

    x1 = np.array([(num_links[num_chains>0]/num_chains[num_chains>0])<=x[i] for i in range(11)])
    y1 = np.array([edge_frac[num_chains>0] >= y[i] for i in range(11)])
    frac_trees = np.zeros((11,11))# fraction of trees across all graphs as a function of edge-frac and nl/nb thresholds
    for i in range(11):
        for j in range(11):
            frac_trees[i,j] = sum(x1[i]&y1[j])/num_all_graphs
            print([i,j])
    plt.imshow(frac_trees, cmap='hot', vmin=0.01, vmax=0.3)
    plt.xlabel('edge_frac');plt.ylabel('nl/nb')
    plt.colorbar()
    return frac_trees

 

def triangle_fraction_plot(Rep):
    f1 = np.array(list(set(Rep['filenum'])))
    X = np.array([[0,0.0,0.0]]*len(f1))
    for i in range(len(f1)):
        Rep1 = Rep.loc[Rep['filenum'] == f1[i]]
        X[i][0] = list(Rep1['N'])[0]
        x1 = Rep1['frac_rep_graph']
        x3 = Rep1['frac_rep_nongraph']
        x2 = (Rep1.graph_nodes*Rep1.frac_rep_graph + Rep1.nongraph_nodes*Rep1.frac_rep_nongraph)/(Rep1.graph_nodes + Rep1.nongraph_nodes)
        X[i][1] = sum((x1-x2)>0)/len(x1)
        X[i][2] = sum((x3-x2)>0)/len(x3)
    color=['xkcd:medium green','xkcd:denim blue','xkcd:gold','xkcd:pale red','grey']
    plt.figure();sns.boxplot(x=X[:,0], y=X[:,2],palette = color);plt.ylim((-0.1,1.1))#lower triangle fraction
    plt.figure();sns.boxplot(x=X[:,0], y=X[:,1],palette = color);plt.ylim((-0.1,1.1))#upper triangle fraction
    return X
  
def get_regeneration_param_lineplots(Rep):
    x2 = (Rep.graph_nodes*Rep.frac_rep_graph + Rep.nongraph_nodes*Rep.frac_rep_nongraph)/(Rep.graph_nodes + Rep.nongraph_nodes) 
    # sig
    plt.figure()
    sns.lineplot(y=x2,x=Rep['sig'])
    a = sns.lineplot(y='frac_rep_graph',x='sig',data=Rep)
    a.set(ylim=(0,1))
    #asym
    plt.figure()
    b = sns.lineplot(y=x2,x=Rep['asym'])
    b = sns.lineplot(y='frac_rep_graph',x='asym',data=Rep)
    b.set(ylim=(0,1))

    #adj
    plt.figure()
    c = sns.lineplot(y=x2,x=Rep['adj'])
    c = sns.lineplot(y='frac_rep_graph',x='adj',data=Rep)
    c.set(ylim=(0,1))

    return 0 

def get_jitter_single_genome(Rep):
    Rep = Rep.loc[Rep['N']==7] 
    unicel =(Rep.scc + Rep.cyc + Rep.chain + Rep.dag + Rep.tree) == 0
    Rep['uni'] = unicel
    f1 = np.array(list(set(Rep['filenum'])))
    for i in range(len(f1)):
        Rep1 = Rep.loc[Rep['filenum'] == f1[i]]
        x2 = (Rep1.graph_nodes*Rep1.frac_rep_graph + Rep1.nongraph_nodes*Rep1.frac_rep_nongraph)/(Rep1.graph_nodes + Rep1.nongraph_nodes)
       # plt.figure()
       # jitterplot.jitter(x = x2,y = Rep1.frac_rep_graph, c = 'k')
       # plt.xlim((-0.1,1.1));plt.ylim((-0.1,1.1))
        # coloured according to topology
        plt.figure()
        top = ['uni','scc','cyc','chain','dag']
        col = list(np.array(Set1_9.mpl_colors)[[6, 3, 2, 0, 1]]) + ['k']
        stdev_x = .03*(max(x2)-min(x2))
        stdev_y = .03*(max(Rep1.frac_rep_graph)-min(Rep1.frac_rep_graph))
        
        for j2,j in enumerate(top):
            tx2 = x2[Rep1[j] == 1]
            tf = Rep1.loc[Rep1[j]==1,'frac_rep_graph']
            tx2 = tx2 + np.random.randn(len(tx2)) * stdev_x
            tf = tf + np.random.randn(len(tf)) * stdev_y
            plt.scatter(tx2, tf, s=0.07, c=col[j2], marker='o', cmap=None, norm=None, vmin=None, vmax=None, alpha=None, linewidths=None)
       # plt.legend(('unicellular','scc','cyclic','chain','DAG'))
        tx2 = x2[Rep1['tree'] == 1]
        tf = Rep1.loc[Rep1['tree']==1,'frac_rep_graph']
        jitterplot.jitter(x = tx2,y = tf, c =col[-1])
        plt.xlim((-0.1,1.1));plt.ylim((-0.1,1.1))
    return 0

def get_acyclic_jitter_single_genome(Rep):
    Rep = Rep.loc[Rep['N']==7] 
    unicel =(Rep.scc + Rep.cyc + Rep.chain + Rep.dag + Rep.tree) == 0
    Rep['uni'] = unicel
    Rep['cyclic'] = (Rep.scc+Rep.cyc) > 0 
    f1 = np.array(list(set(Rep['filenum'])))
    for i in [22]:#range(len(f1)):
        Rep1 = Rep.loc[Rep['filenum'] == f1[i]]
        x2 = (Rep1.graph_nodes*Rep1.frac_rep_graph + Rep1.nongraph_nodes*Rep1.frac_rep_nongraph)/(Rep1.graph_nodes + Rep1.nongraph_nodes)
       # plt.figure()
       # jitterplot.jitter(x = x2,y = Rep1.frac_rep_graph, c = 'k')
       # plt.xlim((-0.1,1.1));plt.ylim((-0.1,1.1))
        # coloured according to topology
        plt.figure()
        top = ['uni','cyclic','chain','dag']
       # col = [np.array([[0.9,0.9,0.9]]),np.array([[0.9,0.9,0.9]]),np.array([[1,0,0]]),np.array([[0,0,0]])]
        col = ['xkcd:gunmetal', 'xkcd:medium grey', 'xkcd:pale red', 'xkcd:sky blue']
        sizes = [0.01, 0.01, 0.3,0.3]
        stdev_x = .03*(max(x2)-min(x2))
        stdev_y = .03*(max(Rep1.frac_rep_graph)-min(Rep1.frac_rep_graph))
        x_fillbet = np.arange(-0.1,1.1,0.1)
        y1_fillbet = np.arange(-0.1,1.1,0.1)-0.1
        y2_fillbet = np.arange(-0.1,1.1,0.1)+0.1   
        plt.fill_between(x_fillbet, y1_fillbet, y2_fillbet, alpha = 0.2, color = [0.2,0.2,0]);plt.axis([-0.1,1.1,-0.1,1.1]) 
        for j2,j in enumerate(top):
            tx2 = x2[Rep1[j] == 1]
            tf = Rep1.loc[Rep1[j]==1,'frac_rep_graph']
            tx2 = tx2 + np.random.randn(len(tx2)) * stdev_x
            tf = tf + np.random.randn(len(tf)) * stdev_y
            plt.scatter(tx2, tf, s=sizes[j2], c=col[j2], marker='o', cmap=None, norm=None, vmin=None, vmax=None, alpha=None, linewidths=None)
       # plt.legend(('unicellular','scc','cyclic','chain','DAG'))
        tx2 = x2[Rep1['tree'] == 1]
        tf = Rep1.loc[Rep1['tree']==1,'frac_rep_graph']
        jitterplot.jitter(x = tx2, y = tf, s=0.7, c='k' )
        plt.xlim((-0.1,1.1));plt.ylim((-0.1,1.1))
    return 0

def get_graph(edge_list):
# get nx graph from edge-list
    edges =  np.array([int(x) for x in edge_list.split(';')])
    edges = edges.reshape(int(edges.size/2), 2)
    G = nx.DiGraph()
    G.add_edges_from(edges)
    return G

def get_acyclic_jitter_single_genome_unique_top(Rep):
    Rep = Rep.loc[Rep['N']==7] 
    unicel =(Rep.scc + Rep.cyc + Rep.chain + Rep.dag + Rep.tree) == 0
    Rep['uni'] = unicel
    Rep['cyclic'] = (Rep.scc+Rep.cyc) > 0 
    f1 = np.array(list(set(Rep['filenum'])))
    for i in [22]:#range(len(f1)):
        Rep1 = Rep.loc[Rep['filenum'] == f1[i]]
        x2 = (Rep1.graph_nodes*Rep1.frac_rep_graph + Rep1.nongraph_nodes*Rep1.frac_rep_nongraph)/(Rep1.graph_nodes + Rep1.nongraph_nodes)
       # plt.figure()
       # jitterplot.jitter(x = x2,y = Rep1.frac_rep_graph, c = 'k')
       # plt.xlim((-0.1,1.1));plt.ylim((-0.1,1.1))
        # coloured according to topology
        plt.figure()
        top = ['uni','cyclic','chain','dag', 'tree']
       # col = [np.array([[0.9,0.9,0.9]]),np.array([[0.9,0.9,0.9]]),np.array([[1,0,0]]),np.array([[0,0,0]])]
        col = ['xkcd:gunmetal', 'xkcd:medium grey', 'xkcd:pale red', 'xkcd:sky blue', 'k']
        sizes = [0.1, 0.1, 3, 3, 3]
        stdev_x = .01*(max(x2)-min(x2))
        stdev_y = .01*(max(Rep1.frac_rep_graph)-min(Rep1.frac_rep_graph))
        
        x_fillbet = np.arange(-0.1,1.1,0.1)
        y1_fillbet = np.arange(-0.1,1.1,0.1)-0.1
        y2_fillbet = np.arange(-0.1,1.1,0.1)+0.1   
        plt.fill_between(x_fillbet, y1_fillbet, y2_fillbet, alpha = 0.2, color = [0.2,0.2,0]);plt.axis([-0.1,1.1,-0.1,1.1]) 
        for j2,j in enumerate(top):
            Rep2 = Rep1.loc[Rep[j]==1].reset_index(drop=True)
            tx2 = x2[Rep1[j] == 1].reset_index(drop=True)
            tf = Rep2.frac_rep_graph
            # find all unique combinations of tx2 and tf
            fg_fp = np.concatenate(([tx2],[tf]),axis=0).T
            uniq_fg_fp = np.unique(fg_fp,axis=0)
            fg = []
            fp = []
            # find all graphs at each (tx2,tf) set
            for k in range(np.shape(uniq_fg_fp)[0]):
                print(k,sum((tx2 == uniq_fg_fp[k][0])&(tf == uniq_fg_fp[k][1])))
               # print(tx2[:10],tf[:10],uniq_fg_fp[k])
                G_all = Rep2.graph[(tx2 == uniq_fg_fp[k][0])&(tf == uniq_fg_fp[k][1])].reset_index(drop=True)# all graphs that have this particular (fg,fp)
                G_uniq = [get_graph(G_all[0])]
                
                for k2 in range(len(G_all)):
                   # print(k2,G_all[k2])
                    G2 = get_graph(G_all[k2])
                    iso=0
                    for k3 in range(len(G_uniq)):
                        G1 = G_uniq[k3]
                        iso +=  nx.is_isomorphic(G1,G2) == 0
                    if iso == len(G_uniq):
                        G_uniq = G_uniq + [G2]
                fg = fg + [uniq_fg_fp[k][0]]*len(G_uniq)
                fp = fp + [uniq_fg_fp[k][1]]*len(G_uniq)
            # reduce each graph set by removing any graphs isomorphic to one already present in the set   
            fg = fg + np.random.randn(len(fg)) * stdev_x
            fp = fp + np.random.randn(len(fp)) * stdev_y
            plt.scatter(fg, fp, s=sizes[j2], c=col[j2], marker='o', cmap=None, norm=None, vmin=None, vmax=None, alpha=None, linewidths=None)
       # plt.legend(('unicellular','scc','cyclic','chain','DAG'))
       # tx2 = x2[Rep1['tree'] == 1]
       # tf = Rep1.loc[Rep1['tree']==1,'frac_rep_graph']
       # jitterplot.jitter(x = tx2, y = tf, s=0.7, c='k' )
        plt.xlim((-0.1,1.1));plt.ylim((-0.1,1.1))
    return 0

def get_graph_nodes_hist(Rep):
    gn = Rep.graph_nodes
    x1 = gn[np.where(Rep.N==3)[0]]
    x2 = gn[np.where(Rep.N==4)[0]]
    x3 = gn[np.where(Rep.N==5)[0]]
    x4 = gn[np.where(Rep.N==6)[0]]
    x5 = gn[np.where(Rep.N==7)[0]]
    bins = np.linspace(0,50,11)
    
    f1 = plt.hist([x1,x2,x3,x4,x5],bins,density=True,color=['xkcd:medium green','xkcd:denim blue','xkcd:gold','xkcd:pale red','black'],label=['N=3','4','5','6','7'])
    plt.legend(loc='upper right')
    plt.xticks(bins)
    ylim = np.sum(f1[0][0])
    plt.yticks(np.linspace(0,ylim,3),('0','0.5','1'))
    return 0
 
def get_graph_nodes_hist_acyclic(Rep):
    gn = Rep.graph_nodes# all graphs
    x1 = gn[(Rep.cyc==0)&(Rep.scc==0)&(Rep.tree==0)&(Rep.chain==0)&(Rep.dag==0)]# unicellular graphs
    x2 = gn[(Rep.cyc>0)|(Rep.scc>0)]# cyclic graphs
    x3 = gn[Rep.chain>0]# chain
    x4 = gn[Rep.tree>0]# tree
    x5 = gn[Rep.dag>0]# dag
    bins = np.linspace(0,50,11)
    
    f1 = plt.hist([x1,x2,x3,x4,x5],bins,density=True,color=['xkcd:dark grey','xkcd:medium grey','xkcd:bright red','black', 'blue'],label=['unicellular','cyclic','chain','tree','DAG'])
    plt.legend(loc='upper right')
    plt.xticks(bins)
    ylim = np.sum(f1[0][0])
    plt.yticks(np.linspace(0,ylim,3),('0','0.5','1'))

    return 0

def get_graph_nodes_bar_acyclic(Rep):
    gn = Rep.graph_nodes
    x1 = gn[(Rep.cyc>0)|(Rep.scc>0)|(Rep.tree>0)|(Rep.chain>0)|(Rep.dag>0)]# all graphs except unicellular
    x2 = gn[(Rep.tree>0)|(Rep.chain>0)|(Rep.dag>0)]# all graphs except unicellular + cyclic    
    x3 = gn[(Rep.tree>0)|(Rep.dag>0)]#trees and dags
    x4 = gn[Rep.dag>0]# dags
    barheights = np.zeros((5,18))
    barwindows = np.linspace(0,90,19)
    
    barheights[0] = [len(gn[(gn>=barwindows[i])&(gn<barwindows[i+1])]) for i in range(len(barwindows)-1)] 
    barheights[1] = [len(x1[(x1>=barwindows[i])&(x1<barwindows[i+1])]) for i in range(len(barwindows)-1)]
    barheights[2] = [len(x2[(x2>=barwindows[i])&(x2<barwindows[i+1])]) for i in range(len(barwindows)-1)]
    barheights[3] = [len(x3[(x3>=barwindows[i])&(x3<barwindows[i+1])]) for i in range(len(barwindows)-1)]
    barheights[4] = [len(x4[(x4>=barwindows[i])&(x4<barwindows[i+1])]) for i in range(len(barwindows)-1)]
    color =['xkcd:gunmetal','xkcd:medium grey','xkcd:pale red','black', 'xkcd:sky blue'];label=['unicellular','cyclic','chain','tree','DAG']
    plt.figure(); 
    for i in range(5):
        plt.bar(np.linspace(0,85,18), barheights[i], color = color[i], log=True, align='edge', width=4.5)#, edgecolor='white')
    plt.legend(label)
    return 0

def get_graph_nodes_hist_topologywise(Rep,top,bins):
    gn = Rep.graph_nodes[Rep[top]==1]
    x1 = gn[np.where(Rep.N==3)[0]]
    x2 = gn[np.where(Rep.N==4)[0]]
    x3 = gn[np.where(Rep.N==5)[0]]
    x4 = gn[np.where(Rep.N==6)[0]]
    x5 = gn[np.where(Rep.N==7)[0]]
    n = max(gn)
    
    
    f1 = plt.hist([x1,x2,x3,x4,x5],bins,density=True,color=['xkcd:medium green','xkcd:denim blue','xkcd:gold','xkcd:pale red','black'],label=['N=3','4','5','6','7'])
    plt.legend(loc='upper right')
    plt.xticks(bins)
    ylim = np.sum(f1[0][0])
    plt.yticks(np.linspace(0,ylim,3),('0','0.5','1'))
    return ylim

def get_node_param_lineplots(Rep):

    color=['xkcd:medium green','xkcd:denim blue','xkcd:gold','xkcd:pale red','black']
    sns.lineplot(x='asym',y='graph_nodes',data=Rep ,hue='N', palette=color)
    plt.figure();sns.lineplot(x='sig',y='graph_nodes',data=Rep ,hue='N', palette=color)
    plt.figure();sns.lineplot(x='adj',y='graph_nodes',data=Rep ,hue='N', palette=color)
    plt.figure();f1 = sns.lineplot(x='graph_nodes',y='num_edges',data=Rep ,hue='N',hue_order=[7,6,5,4,3], palette=color[::-1])#[::-1] reverses the list
    f1.set(yscale='log')
    return 0

def get_node_edge_scatter_jitter_plots(Rep):
    #scatter plot:
    plt.figure()
    plt.scatter(x = Rep.graph_nodes[Rep.N==7],y = np.log10(Rep.num_edges[Rep.N==7]),s=0.1,c ='k',alpha=0.1)
    plt.scatter(x = Rep.graph_nodes[Rep.N==6],y = np.log10(Rep.num_edges[Rep.N==6]),s=0.1,c ='xkcd:pale red',alpha=0.1)
    plt.scatter(x = Rep.graph_nodes[Rep.N==5],y = np.log10(Rep.num_edges[Rep.N==5]),s=0.1,c ='xkcd:gold',alpha=0.1)
    plt.scatter(x = Rep.graph_nodes[Rep.N==4],y = np.log10(Rep.num_edges[Rep.N==4]),s=0.1,c ='xkcd:denim blue',alpha=0.1)
    plt.scatter(x = Rep.graph_nodes[Rep.N==3],y = np.log10(Rep.num_edges[Rep.N==3]),s=0.1,c ='xkcd:medium green',alpha=0.1)

    #jitter plot:
    plt.figure()
    jitterplot.jitter(x = Rep.graph_nodes[Rep.N==7],y = Rep.num_edges[Rep.N==7],s=0.1,c ='k')
    jitterplot.jitter(x = Rep.graph_nodes[Rep.N==6],y = Rep.num_edges[Rep.N==6],s=0.1,c ='xkcd:pale red')
    jitterplot.jitter(x = Rep.graph_nodes[Rep.N==5],y = Rep.num_edges[Rep.N==5],s=0.1,c ='xkcd:gold')
    jitterplot.jitter(x = Rep.graph_nodes[Rep.N==4],y = Rep.num_edges[Rep.N==4],s=0.1,c ='xkcd:denim blue')
    jitterplot.jitter(x = Rep.graph_nodes[Rep.N==3],y = Rep.num_edges[Rep.N==3],s=0.1,c ='xkcd:medium green')

def get_topology_regeneration_swarms(Rep):
    unicel =(Rep.scc + Rep.cyc + Rep.chain + Rep.dag + Rep.tree) == 0
    Rep['unicell'] = unicel
    f1 = np.array(list(set(Rep['filenum'])))
    t1 = ['unicell','scc','cyc','chain','dag','tree']
    X = np.array([[0.0,0.0]]*(6*len(f1)))
    T1 = t1*len(f1)
    count = 0
    for i in range(len(f1)):
        for j in t1:
            Rep1 = Rep.loc[(Rep['filenum'] == f1[i])&(Rep[j]==1)]
            if Rep1.empty == 0:
                x1 = Rep1['frac_rep_graph']
                x3 = Rep1['frac_rep_nongraph']
                x2 = (Rep1.graph_nodes*Rep1.frac_rep_graph + Rep1.nongraph_nodes*Rep1.frac_rep_nongraph)/(Rep1.graph_nodes + Rep1.nongraph_nodes)
                X[count][0] = sum((x1-x2)>0)/len(x1)
                X[count][1] = sum((x3-x2)>0)/len(x3)
            count = count+1
    col = list(np.array(Set1_9.mpl_colors)[[6, 3, 2, 0, 1]]) + ['grey']
   #col = col*len(f1)
    plt.figure();sns.boxplot(x=T1, y=X[:,1],palette = col);plt.ylim((-0.1,1.1))#lower triangle fraction
    plt.figure();sns.boxplot(x=T1, y=X[:,0],palette = col);plt.ylim((-0.1,1.1))#upper triangle fraction


def get_acyclic_regeneration_swarms(Rep):
    unicel =(Rep.scc + Rep.cyc + Rep.chain + Rep.dag + Rep.tree) == 0
    Rep['unicell'] = unicel
    Rep['cyclic'] = Rep.cyc + Rep.scc
    f1 = np.array(list(set(Rep['filenum'])))
    t1 = ['unicell','cyclic','chain','tree','dag']
    X = np.array([[0.0,0.0]]*(5*len(f1)))
    T1 = t1*len(f1)
    count = 0
    for i in range(len(f1)):
        for j in t1:
            Rep1 = Rep.loc[(Rep['filenum'] == f1[i])&(Rep[j]==1)]
            if Rep1.empty == 0:
                x1 = Rep1['frac_rep_graph']
                x3 = Rep1['frac_rep_nongraph']
                x2 = (Rep1.graph_nodes*Rep1.frac_rep_graph + Rep1.nongraph_nodes*Rep1.frac_rep_nongraph)/(Rep1.graph_nodes + Rep1.nongraph_nodes)
                X[count][0] = sum((x1-x2)>0)/len(x1)
                X[count][1] = sum((x3-x2)>0)/len(x3)
            count = count+1
    col = ['xkcd:gunmetal', 'xkcd:medium grey','xkcd:pale red','k','xkcd:sky blue']
    edgecol = ['k', 'k', 'k', 'xkcd:medium grey', 'k']
   #col = col*len(f1)
    plt.figure();ax = sns.boxplot(x=T1, y=X[:,1],palette = col);plt.ylim((-0.1,1.1))#lower triangle fraction
    for i,box in enumerate(ax.artists):
        box.set_edgecolor(edgecol[i])
        for j in range(6*i,6*(i+1)):
            ax.lines[j].set_color(edgecol[i])
    plt.figure();ax = sns.boxplot(x=T1, y=X[:,0],palette = col);plt.ylim((-0.1,1.1))#upper triangle fraction
    for i,box in enumerate(ax.artists):
        box.set_edgecolor(edgecol[i])
        for j in range(6*i,6*(i+1)):
            ax.lines[j].set_color(edgecol[i])


def get_acyclic_acyclized_regeneration_swarms(Rep, edge_frac2, num_links2, num_chains2):
    Rep2 = pd.DataFrame(columns = ['unicell', 'cyclic', 'chain', 'tree', 'dag'])
    print([sum(Rep.tree), sum(Rep.chain), sum(Rep.dag)])
    acyc_thresh = 0.5
    treeness_thresh = 0.5

    Rep2['unicell'] = Rep.unicell
    Rep2['cyclic'] = (Rep.cyc + Rep.scc)*(edge_frac2<acyc_thresh)
    Rep2['chain'] = (edge_frac2>=acyc_thresh)&(num_links2==0)&(num_chains2==0)
    
    Rep2['tree'] = (num_chains2>0)&(edge_frac2 >= acyc_thresh)
    Rep2['tree'][num_links2 > num_chains2 * treeness_thresh] = 0

    Rep2['dag']= (num_chains2>0)&(edge_frac2 >= acyc_thresh)
    Rep2['dag'][num_links2 <= num_chains2 * treeness_thresh] = 0    
    print(sum(Rep2['tree']), sum(Rep2.chain), sum(Rep2['dag']))    

    f1 = np.array(list(set(Rep['filenum'])))
    t1 = ['unicell','cyclic','chain','tree','dag']
    X = np.array([[0.0,0.0]]*(5*len(f1)))
    T1 = t1*len(f1)
    count = 0
    for i in range(len(f1)):
        for j in t1:
            Rep1 = Rep.loc[(Rep['filenum'] == f1[i])&(Rep2[j]==1)]
            if Rep1.empty == 0:
                x1 = Rep1['frac_rep_graph']
                x3 = Rep1['frac_rep_nongraph']
                x2 = (Rep1.graph_nodes*Rep1.frac_rep_graph + Rep1.nongraph_nodes*Rep1.frac_rep_nongraph)/(Rep1.graph_nodes + Rep1.nongraph_nodes)
                X[count][0] = sum((x1-x2)>0)/len(x1)
                X[count][1] = sum((x3-x2)>0)/len(x3)
            count = count+1
    col = ['xkcd:gunmetal', 'xkcd:medium grey','xkcd:pale red','k','xkcd:sky blue']
    edgecol = ['k', 'k', 'k', 'xkcd:medium grey', 'k']
   #col = col*len(f1)
    plt.figure();ax = sns.boxplot(x=T1, y=X[:,1],palette = col);plt.ylim((-0.1,1.1))#lower triangle fraction
    for i,box in enumerate(ax.artists):
        box.set_edgecolor(edgecol[i])
        for j in range(6*i,6*(i+1)):
            ax.lines[j].set_color(edgecol[i])
    plt.figure();ax = sns.boxplot(x=T1, y=X[:,0],palette = col);plt.ylim((-0.1,1.1))#upper triangle fraction
    for i,box in enumerate(ax.artists):
        box.set_edgecolor(edgecol[i])
        for j in range(6*i,6*(i+1)):
            ax.lines[j].set_color(edgecol[i])
    return X

def get_topology_regenerative_capacity_box(Rep):
    unicel =(Rep.scc + Rep.cyc + Rep.chain + Rep.dag + Rep.tree) == 0
    Rep['unicell'] = unicel
    #f1 = np.array(list(set(Rep['filenum'])))
    T1 = ['unicell','scc','cyc','chain','dag','tree']
    order = ['u', 's', 'C', 'c', 'd', 't']
    topology = np.array(['']*len(Rep))
    x2 = (Rep.graph_nodes*Rep.frac_rep_graph + Rep.nongraph_nodes*Rep.frac_rep_nongraph)/(Rep.graph_nodes + Rep.nongraph_nodes)
    reg_cap = Rep['frac_rep_graph']/x2
    for i,j in enumerate(T1):
        topology[np.where(Rep[j]==1)[0]]=order[i]
    col = list(np.array(Set1_9.mpl_colors)[[6, 3, 2, 0, 1]]) + ['grey']
    plt.figure();sns.boxplot(x=topology, y=np.log10(reg_cap),palette = col,order=order);plt.ylim((-1.2,2.1))#regenerative capacity box plot
    plt.yticks([-1,0,1,2],['0.1','1','10','100'])    

def get_acyclic_regenerative_capacity_box(Rep):
    unicel =(Rep.scc + Rep.cyc + Rep.chain + Rep.dag + Rep.tree) == 0
    Rep['unicel'] = unicel
    Rep['cyclic'] = Rep.cyc + Rep.scc
    #f1 = np.array(list(set(Rep['filenum'])))
    T1 = ['unicel', 'cyclic','chain','tree','dag']
    order = ['u', 'C', 'c', 'd', 't']
    topology = np.array(['']*len(Rep))
    x2 = (Rep.graph_nodes*Rep.frac_rep_graph + Rep.nongraph_nodes*Rep.frac_rep_nongraph)/(Rep.graph_nodes + Rep.nongraph_nodes)
    reg_cap = Rep['frac_rep_graph']/x2
    for i,j in enumerate(T1):
        topology[np.where(Rep[j]==1)[0]]=order[i]
    col = ['xkcd:gunmetal grey', 'lightgrey', 'xkcd: pale red', 'k', 'xkcd:sky blue']
    plt.figure();sns.boxplot(x=topology, y=np.log10(reg_cap), palette = col, order=order); plt.ylim((-1.2,2.1))#regenerative capacity box plot
    plt.yticks([-1,0,1,2],['0.1','1','10','100'])    

def stacked_burger_plot(Rep):
# get stacked bar plot for graph topology abundance
    unicel =(Rep.scc + Rep.cyc + Rep.chain + Rep.dag + Rep.tree) == 0
    Cyc = (Rep.scc|Rep.cyc)>0
    Rep['unicell'] = unicel
    Rep['Cyc'] = Cyc
    #order = ['u', 's', 'C', 'c', 'd', 't']
    order = ['u', 'C', 'c', 't', 'd']
    #top = ['unicell','scc','cyc','chain','dag','tree']
    top = ['unicell', 'Cyc', 'chain', 'tree', 'dag']
    topology = np.array(['']*len(Rep))
    for i,o in enumerate(top):
        topology[np.where(Rep[o]==1)[0]]=order[i]# only takes in a single alphabet
    Rep['topology'] = topology

    n_counts = [Counter(Rep.loc[Rep.N==n, 'topology']) for n in range(3,8)]

    X = range(3,8)
    Y = np.array([[n_counts[i][o] for i in range(5)] for o in order], dtype=float)
    for i in range(5):
        Y[:,i] = Y[:,i] / Y[:,i].sum()
    width = 0.75
    #col = list(np.array(Set1_9.mpl_colors)[[6, 3, 2, 0, 1]])+['k']
    col = ['xkcd:gunmetal', 'xkcd:medium grey', 'xkcd:pale red', 'k', 'xkcd:sky blue']
# col = ['xkcd:gunmetal grey', 'lightgrey', 'xkcd:pale red', 'k', 'xkcd:sky blue']
    fig, ax = plt.subplots()
    bottom = np.zeros(5)
    for i in range(len(order)):
        print(col[i])
        ax.bar(X, Y[i], width, bottom=bottom, color=col[i], edgecolor='white', linewidth=0.5)
        bottom += Y[i]
    plt.xticks([3,4,5,6,7],('3','4','5','6','7'))
    plt.xlim((2.5,7.5))

def stacked_burger_plot_genomewise(Rep):
# get stacked bar plot for graph topology abundance
    Rep = Rep[Rep.N==7].reset_index(drop=True)
    unicel =(Rep.scc + Rep.cyc + Rep.chain + Rep.dag + Rep.tree) == 0
    Rep['unicell'] = unicel
    order = ['u', 's', 'C', 'c', 'd', 't']
    top = ['unicell','scc','cyc','chain','dag','tree']
    topology = np.array(['']*len(Rep))
    for i,o in enumerate(top):
        topology[np.where(Rep[o]==1)[0]]=order[i]# only takes in a single alphabet
    Rep['topology'] = topology
    fnum = list(set(Rep.filenum))# list of distinct genomes for N=7
    n_counts = [Counter(Rep.loc[Rep.filenum==n, 'topology']) for n in fnum]

    X = range(len(fnum))
    Y = np.array([[n_counts[i][o] for i in range(len(fnum))] for o in order], dtype=float)
    for i in range(len(fnum)):
        Y[:,i] = Y[:,i] / Y[:,i].sum()
    width = 0.75

    col = list(np.array(Set1_9.mpl_colors)[[6, 3, 2, 0, 1]]) + ['k']
    fig, ax = plt.subplots()
    bottom = np.zeros(len(fnum))
    for i in range(len(order)):
        ax.bar(X, Y[i], width, bottom=bottom, color=col[i], edgecolor='white', linewidth=0.5)
        bottom += Y[i]
   # plt.xticks([3,4,5,6,7],('3','4','5','6','7'))
   # plt.xlim((2.5,7.5))

def stacked_acyclic_plot(Rep):
    # stacked bar for acyclic graph types
    unicel =(Rep.scc + Rep.cyc + Rep.chain + Rep.dag + Rep.tree) == 0
    Rep['unicel']=unicel
    col = ['xkcd:gunmetal', 'lightgrey', 'xkcd: pale red', 'k', 'xkcd:sky blue']
    n_counts1 = np.array([sum(Rep.unicel>0), sum((Rep.cyclic>0)|(Rep.scc>0)), sum(Rep.chain>0), sum(Rep.tree>0), sum(Rep.dag>0)])
    fig, ax = plt.subplots()
    for i in range(len(n_counts1)):
        ax.bar(0, np.sum(n_counts1[i:])/np.sum(n_counts1), color=col[i], edgecolor='white', linewidth=0.5)
    plt.xlim((-0.1,0.1))
    
    col2 = [[0.8,0.8,0.8],[0.5,0.5,0.5],[0.3,0.3,0.3]]#[acyclic, cyclic, unicellular]
    n_counts2 = np.array([np.sum(n_counts1), sum((Rep.cyc>0)|(Rep.scc>0)), sum((Rep.cyc==0)&(Rep.scc==0)&(Rep.tree==0)&(Rep.dag==0)&(Rep.chain==0))])
    fig, ax = plt.subplots()
    for i in range(len(n_counts2)):
        ax.bar(0, np.sum(n_counts2[i:])/np.sum(n_counts2), color=col2[i], edgecolor='white', linewidth=0.5)
    plt.xlim((-0.1,0.1))
    return n_counts1/np.sum(n_counts1), n_counts2/np.sum(n_counts2)



def topology_parm_heatmaps(Rep):

    tree_param = Rep.loc[Rep['tree'] == 1,['sig','asym','adj']]
    dag_param = Rep.loc[Rep['dag'] == 1,['sig','asym','adj']]
    chain_param = Rep.loc[Rep['chain'] == 1,['sig','asym','adj']]
    scc_param = Rep.loc[Rep['scc'] == 1,['sig','asym','adj']]
    cyc_param = Rep.loc[Rep['cyc'] == 1,['sig','asym','adj']]
    uni_param = Rep.loc[(Rep['tree']+Rep['dag']+Rep['chain']+Rep['scc']+Rep['cyc']) == 0,['sig','asym','adj']]

    # sig_adj
    tree_sig_adj = np.zeros((11,11))
    dag_sig_adj = np.zeros((11,11))
    chain_sig_adj = np.zeros((11,11))
    scc_sig_adj = np.zeros((11,11))
    cyc_sig_adj = np.zeros((11,11))
    uni_sig_adj = np.zeros((11,11))
    # sig_asym
    tree_sig_asym = np.zeros((11,11))
    dag_sig_asym = np.zeros((11,11))
    chain_sig_asym = np.zeros((11,11))
    scc_sig_asym = np.zeros((11,11))
    cyc_sig_asym = np.zeros((11,11))
    uni_sig_asym = np.zeros((11,11))
    #adj_asym
    tree_adj_asym = np.zeros((11,11))
    dag_adj_asym = np.zeros((11,11))
    chain_adj_asym = np.zeros((11,11))
    scc_adj_asym = np.zeros((11,11))
    cyc_adj_asym = np.zeros((11,11))
    uni_adj_asym = np.zeros((11,11))

    for i in range(11):
        for j in range(11):
            #sig_adj
            tree_sig_adj[i,j] = len(tree_param.loc[(tree_param['sig'] == i+1)&(tree_param['adj'] == j+1)])/len(tree_param)
            dag_sig_adj[i,j] = len(dag_param.loc[(dag_param['sig'] == i+1)&(dag_param['adj'] == j+1)])/len(dag_param)
            chain_sig_adj[i,j] = len(chain_param.loc[(chain_param['sig'] == i+1)&(chain_param['adj'] == j+1)])/len(chain_param)
            scc_sig_adj[i,j] = len(scc_param.loc[(scc_param['sig'] == i+1)&(scc_param['adj'] == j+1)])/len(scc_param)
            cyc_sig_adj[i,j] = len(cyc_param.loc[(cyc_param['sig'] == i+1)&(cyc_param['adj'] == j+1)])/len(cyc_param)
            uni_sig_adj[i,j] = len(uni_param.loc[(uni_param['sig'] == i+1)&(uni_param['adj'] == j+1)])/len(uni_param)
            #sig_asym
            tree_sig_asym[i,j] = len(tree_param.loc[(tree_param['sig'] == i+1)&(tree_param['asym'] == j+1)])/len(tree_param)
            dag_sig_asym[i,j] = len(dag_param.loc[(dag_param['sig'] == i+1)&(dag_param['asym'] == j+1)])/len(dag_param)
            chain_sig_asym[i,j] = len(chain_param.loc[(chain_param['sig'] == i+1)&(chain_param['asym'] == j+1)])/len(chain_param)
            scc_sig_asym[i,j] = len(scc_param.loc[(scc_param['sig'] == i+1)&(scc_param['asym'] == j+1)])/len(scc_param)
            cyc_sig_asym[i,j] = len(cyc_param.loc[(cyc_param['sig'] == i+1)&(cyc_param['asym'] == j+1)])/len(cyc_param)
            uni_sig_asym[i,j] = len(uni_param.loc[(uni_param['sig'] == i+1)&(uni_param['asym'] == j+1)])/len(uni_param)
            #adj_asym
            tree_adj_asym[i,j] = len(tree_param.loc[(tree_param['adj'] == i+1)&(tree_param['asym'] == j+1)])/len(tree_param)
            dag_adj_asym[i,j] = len(dag_param.loc[(dag_param['adj'] == i+1)&(dag_param['asym'] == j+1)])/len(dag_param)
            chain_adj_asym[i,j] = len(chain_param.loc[(chain_param['adj'] == i+1)&(chain_param['asym'] == j+1)])/len(chain_param)
            scc_adj_asym[i,j] = len(scc_param.loc[(scc_param['adj'] == i+1)&(scc_param['asym'] == j+1)])/len(scc_param)
            cyc_adj_asym[i,j] = len(cyc_param.loc[(cyc_param['adj'] == i+1)&(cyc_param['asym'] == j+1)])/len(cyc_param)
            uni_adj_asym[i,j] = len(uni_param.loc[(uni_param['adj'] == i+1)&(uni_param['asym'] == j+1)])/len(uni_param)

    plt.figure();plt.imshow(tree_sig_asym,cmap = 'hot');plt.colorbar();plt.clim(0,0.04);plt.title('x=sig,y=asym')
    plt.figure();plt.imshow(tree_sig_adj,cmap = 'hot');plt.colorbar();plt.clim(0,0.04);plt.title('x=sig,y=adj')
    plt.figure();plt.imshow(tree_adj_asym,cmap = 'hot');plt.colorbar();plt.clim(0,0.04);plt.title('x=adj,y=asym')    

    plt.figure();plt.imshow(dag_sig_asym,cmap = 'hot');plt.colorbar();plt.clim(0,0.02);plt.title('x=sig,y=asym')
    plt.figure();plt.imshow(dag_sig_adj,cmap = 'hot');plt.colorbar();plt.clim(0,0.02);plt.title('x=sig,y=adj')
    plt.figure();plt.imshow(dag_adj_asym,cmap = 'hot');plt.colorbar();plt.clim(0,0.02);plt.title('x=adj,y=asym')    

    plt.figure();plt.imshow(chain_sig_asym,cmap = 'hot');plt.colorbar();plt.clim(0,0.02);plt.title('x=sig,y=asym')
    plt.figure();plt.imshow(chain_sig_adj,cmap = 'hot');plt.colorbar();plt.clim(0,0.02);plt.title('x=sig,y=adj')
    plt.figure();plt.imshow(chain_adj_asym,cmap = 'hot');plt.colorbar();plt.clim(0,0.02);plt.title('x=adj,y=asym')

    plt.figure();plt.imshow(scc_sig_asym,cmap = 'hot');plt.colorbar();plt.clim(0,0.02);plt.title('x=sig,y=asym')
    plt.figure();plt.imshow(scc_sig_adj,cmap = 'hot');plt.colorbar();plt.clim(0,0.02);plt.title('x=sig,y=adj')
    plt.figure();plt.imshow(scc_adj_asym,cmap = 'hot');plt.colorbar();plt.clim(0,0.02);plt.title('x=adj,y=asym')

    plt.figure();plt.imshow(cyc_sig_asym,cmap = 'hot');plt.colorbar();plt.clim(0,0.02);plt.title('x=sig,y=asym')
    plt.figure();plt.imshow(cyc_sig_adj,cmap = 'hot');plt.colorbar();plt.clim(0,0.02);plt.title('x=sig,y=adj')
    plt.figure();plt.imshow(cyc_adj_asym,cmap = 'hot');plt.colorbar();plt.clim(0,0.02);plt.title('x=adj,y=asym')

    plt.figure();plt.imshow(uni_sig_asym,cmap = 'hot');plt.colorbar();plt.clim(0,0.02);plt.title('x=sig,y=asym')
    plt.figure();plt.imshow(uni_sig_adj,cmap = 'hot');plt.colorbar();plt.clim(0,0.02);plt.title('x=sig,y=adj')
    plt.figure();plt.imshow(uni_adj_asym,cmap = 'hot');plt.colorbar();plt.clim(0,0.02);plt.title('x=adj,y=asym')


def get_indeg_outdeg_plots(indeg2,outdeg2,cyc, scc, unicel, dag, tree, chain, numcomp):
    t1=np.array(tree)[np.array(numcomp)==1]
    d1=np.array(dag)[np.array(numcomp)==1]
    c1=np.array(chain)[np.array(numcomp)==1]
    c2=np.array(cyc)[np.array(numcomp)==1]
    s1=np.array(scc)[np.array(numcomp)==1]
    u1=np.array(unicel)[np.array(numcomp)==1]
    indeg3 = [indeg2[i] for i in list(np.where(np.array(numcomp)==1)[0])]
    outdeg3 = [outdeg2[i] for i in list(np.where(np.array(numcomp)==1)[0])]


    tree_where = np.where(t1==1)[0][0:25000]
    dag_where = np.where(d1==1)[0][0:25000]
    chain_where = np.where(c1==1)[0][0:25000]
    cyc_where = np.where(c2==1)[0][0:25000]
    scc_where = np.where(s1==1)[0][0:25000]
    unicel_where = np.where(u1==1)[0][0:25000]

    indeg_tree = list(itertools.chain(*[indeg3[i] for i in list(tree_where)]))
    outdeg_tree =list(itertools.chain(*[outdeg3[i] for i in list(tree_where)]))
    indeg_dag = list(itertools.chain(*[indeg3[i] for i in list(dag_where)]))
    outdeg_dag =list(itertools.chain(*[outdeg3[i] for i in list(dag_where)]))
    indeg_chain = list(itertools.chain(*[indeg3[i] for i in list(chain_where)]))
    outdeg_chain =list(itertools.chain(*[outdeg3[i] for i in list(chain_where)]))
    indeg_scc = list(itertools.chain(*[indeg3[i] for i in list(scc_where)]))
    outdeg_scc =list(itertools.chain(*[outdeg3[i] for i in list(scc_where)]))
    indeg_cyc = list(itertools.chain(*[indeg3[i] for i in list(cyc_where)]))
    outdeg_cyc =list(itertools.chain(*[outdeg3[i] for i in list(cyc_where)]))
    indeg_uni = [1]*25000
    outdeg_uni = [1]*25000
             
    col = list(np.array(Set1_9.mpl_colors)[[6, 3, 2, 0, 1]]) + ['k']
    plt.figure();plt.scatter(indeg_tree + np.random.randn(len(indeg_tree))*0.15,outdeg_tree+np.random.randn(len(indeg_tree))*0.15,s=0.03, c =col[5]);plt.xlim((-0.5,5.5));plt.ylim((-0.5,5.5))
    plt.figure();plt.scatter(indeg_dag+np.random.randn(len(indeg_dag))*0.15,outdeg_dag+np.random.randn(len(indeg_dag))*0.15,s=0.03, c =col[4]);plt.xlim((-0.5,5.5));plt.ylim((-0.5,5.5))
    plt.figure();plt.scatter(indeg_chain+np.random.randn(len(indeg_chain))*0.15,outdeg_chain+np.random.randn(len(indeg_chain))*0.15,s=0.03, c =col[3]);plt.xlim((-0.5,5.5));plt.ylim((-0.5,5.5))
    plt.figure();plt.scatter(indeg_cyc+np.random.randn(len(indeg_cyc))*0.15,outdeg_cyc+np.random.randn(len(indeg_cyc))*0.15,s=0.03, c =col[2]);plt.xlim((-0.5,5.5));plt.ylim((-0.5,5.5))
    plt.figure();plt.scatter(indeg_scc+np.random.randn(len(indeg_scc))*0.15,outdeg_scc+np.random.randn(len(indeg_scc))*0.15,s=0.03, c =col[1]);plt.xlim((-0.5,5.5));plt.ylim((-0.5,5.5))
    plt.figure();plt.scatter(indeg_uni+np.random.randn(len(indeg_uni))*0.15,outdeg_uni+np.random.randn(len(indeg_uni))*0.15,s=0.03, c =col[0]);plt.xlim((-0.5,5.5));plt.ylim((-0.5,5.5))
    return 0

def get_celltype_stack(Ind,P):
    graph_nodes = Ind.Ab + Ind.ABCDE + Ind.ABcDE + Ind.ABcdE + Ind.ABCDe + Ind.ABcDe + Ind.ABcde + Ind.ab + Ind.aBCDE + Ind.aBcDE + Ind.aBcdE + Ind.aBCDe + Ind.aBcDe + Ind.aBcde

# get stacked bar plot for pluripotent cell types
    # according to model param
    Ind['AB'] = Ind.ABCDe + Ind.ABCDE + Ind.ABcDe + Ind.ABcde + Ind.ABcDE + Ind.ABcdE
    Ind['aB'] = Ind.aBCDe + Ind.aBCDE + Ind.aBcDE + Ind.aBcDe + Ind.aBcdE + Ind.aBcde
    top = ['ab','aB','Ab','AB']
    X = range(1,12)
    Y = np.zeros((4,11))
    for x in [0,1,2,3]:
        for y in list(range(11)):
            Y[x,y] = np.sum(Ind[top[x]][Ind[P] == y+1])# change in this line for asym or adj or sig
    for i in range(11):
        Y[:,i] = Y[:,i] / Y[:,i].sum()
    width = 0.75

    #col = ['xkcd:black','xkcd:dark grey','xkcd:medium grey','xkcd:brown','xkcd:indigo','xkcd:teal']
    col = ['xkcd:black','xkcd:medium grey','xkcd:brown','xkcd:teal']
    fig, ax = plt.subplots()
    bottom = np.zeros(11)
    for i in range(len(top)):
        ax.bar(X, Y[i], width, bottom=bottom, color=col[i], edgecolor='white', linewidth=0.5)
        bottom += Y[i]
    plt.xticks([1,2,3,4,5,6,7,8,9,10,11],('0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0'))
    plt.xlim((0.5,11.5))

def stacked_independance_plot_graphtop(Ind,Rep):
    unicel =(Rep.scc + Rep.cyc + Rep.chain + Rep.dag + Rep.tree) == 0
    Rep['unicell'] = unicel
    Ind['AB'] = Ind.ABCDe + Ind.ABCDE + Ind.ABcDe + Ind.ABcde + Ind.ABcDE + Ind.ABcdE 
    Ind['aB'] = Ind.aBCDe + Ind.aBCDE + Ind.aBcDE + Ind.aBcDe + Ind.aBcdE + Ind.aBcde 
    top = ['ab','aB','Ab','AB'] # [not reg- not indep, indep | reg- not indep, indep]

    graph_top = ['unicell','scc','cyc','chain','dag','tree'] 
    X = range(1,7) 
    Y = np.zeros((4,6)) 
    for x in [0,1,2,3]: 
        for y in list(range(6)): 
            Y[x,y] = np.sum(Ind[top[x]][Rep[graph_top[y]] == 1])
   
     
    for i in range(6): 
       #print(Y[:,i].sum())
        Y[:,i] = Y[:,i] / Y[:,i].sum() 
    width = 0.75 
   
    col = ['xkcd:black','xkcd:medium grey','xkcd:brown','xkcd:teal']
    fig, ax = plt.subplots() 
    bottom = np.zeros(6) 
    for i in range(len(top)): 
        ax.bar(X, Y[i], width, bottom=bottom, color=col[i], edgecolor='white', linewidth=0.5) 
        bottom += Y[i] 
    plt.xticks([1,2,3,4,5,6,7,8,9,10,11],('0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0')) 
    plt.xlim((0.5,6.5)) 
    return 0

def stacked_independance_plot_acyclic(Ind,Rep):
    cyclic =(Rep.chain + Rep.dag + Rep.tree) == 0
    Rep['cyclic'] = cyclic
    Ind['AB'] = Ind.ABCDe + Ind.ABCDE + Ind.ABcDe + Ind.ABcde + Ind.ABcDE + Ind.ABcdE 
    Ind['aB'] = Ind.aBCDe + Ind.aBCDE + Ind.aBcDE + Ind.aBcDe + Ind.aBcdE + Ind.aBcde 
    top = ['ab','aB','Ab','AB'] # [not reg- not indep, indep | reg- not indep, indep]

    graph_top = ['cyclic','chain','tree','dag'] 
    X = range(1,5) 
    Y = np.zeros((4,4)) 
    for x in [0,1,2,3]: 
        for y in list(range(4)): 
            Y[x,y] = np.sum(Ind[top[x]][Rep[graph_top[y]] == 1])
   
     
    for i in range(4): 
       #print(Y[:,i].sum())
        Y[:,i] = Y[:,i] / Y[:,i].sum() 
    width = 0.75 
   
    col = ['xkcd:black','xkcd:medium grey','xkcd:brown','xkcd:puce']
    fig, ax = plt.subplots() 
    bottom = np.zeros(4) 
    for i in range(len(top)): 
        ax.bar(X, Y[i], width, bottom=bottom, color=col[i], edgecolor='white', linewidth=0.5) 
        bottom += Y[i] 
    plt.xticks([1,2,3,4],('0','0.1','0.2','0.3')) 
    plt.xlim((0.5,4.5)) 
    return 0

def intrinsic_independance_heatmap(Int,Ind):
    int_indep_nongraph = Int.CDE+Int.CDe
    notint_indep_nongraph = Int.cde+Int.cdE+Int.cDe+Int.cDE
    x1 = int_indep_nongraph/(int_indep_nongraph+notint_indep_nongraph)

    graph_nodes = Ind.Ab + Ind.ABCDE + Ind.ABcDE + Ind.ABcdE + Ind.ABCDe + Ind.ABcDe + Ind.ABcde + Ind.ab + Ind.aBCDE + Ind.aBcDE + Ind.aBcdE + Ind.aBCDe + Ind.aBcDe + Ind.aBcde
    int_indep_graph = Ind.ABCDE + Ind.ABCDe + Ind.aBCDE +  Ind.aBCDe
    y1 = int_indep_graph/graph_nodes

    x2 = (int_indep_nongraph + int_indep_graph)/(graph_nodes + int_indep_nongraph + notint_indep_nongraph)

    HM = np.zeros((10,10));
    frac = np.linspace(0,1,11)
    for i in range(10):
        for j in range(10):
            HM[i,j] = sum((x1>frac[int(i)])&(x1<=frac[int(i+1)])&(y1>frac[int(j)])&(y1<=frac[int(j+1)]))/len(x1)
            print([i,j])
    HM[0,0] = HM[0,0] + sum((x1==0)&(y1==0))/len(x1)
    plt.figure();plt.imshow(HM,cmap = 'hot');plt.colorbar();plt.show()
    
def get_intrinsic_indep_stack(Ind,Int,P):
    graph_nodes = Ind.Ab + Ind.ABCDE + Ind.ABcDE + Ind.ABcdE + Ind.ABCDe + Ind.ABcDe + Ind.ABcde + Ind.ab + Ind.aBCDE + Ind.aBcDE + Ind.aBcdE + Ind.aBCDe + Ind.aBcDe + Ind.aBcde

# get stacked bar plot for pluripotent cell types
    # according to model param
    Ind['int_indep_nongraph'] = np.array(Int.CDE)+np.array(Int.CDe)
    Ind['notint_indep_nongraph'] = np.array(Int.cde) + np.array(Int.cdE) + np.array(Int.cDe) + np.array(Int.cDE) 
    Ind['int_indep_graph'] = np.array(Ind.ABCDE) + np.array(Ind.ABCDe) + np.array(Ind.aBCDE) +  np.array(Ind.aBCDe)
    Ind['notint_indep_graph'] = np.array(Ind.Ab) + np.array(Ind.ABcDE) + np.array(Ind.ABcdE) + np.array(Ind.ABcDe) + np.array(Ind.ABcde) + np.array(Ind.ab) + np.array(Ind.aBcDE) + np.array(Ind.aBcdE) + np.array(Ind.aBcDe) + np.array(Ind.aBcde) 
    
    top = ['notint_indep_nongraph','int_indep_nongraph','notint_indep_graph','int_indep_graph']
    X = range(1,12)
    Y = np.zeros((4,11))
    for x in [0,1,2,3]:
        for y in list(range(11)):
            Y[x,y] = np.sum(Ind[top[x]][np.array(Ind[P]) == y+1])# change in this line for asym or adj or sig
    for i in range(11):
        Y[:,i] = Y[:,i] / Y[:,i].sum()
    width = 0.75

    col = ['xkcd:maroon','xkcd:scarlet','xkcd:dark olive','xkcd:sage green']
    fig, ax = plt.subplots()
    bottom = np.zeros(11)
    for i in range(len(top)):
        ax.bar(X, Y[i], width, bottom=bottom, color=col[i], edgecolor='white', linewidth=0.5)
        bottom += Y[i]
    plt.xticks([1,2,3,4,5,6,7,8,9,10,11],('0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0'))
    plt.xlim((0.5,11.5))

def get_intrinsic_indep_stack_N(Ind,Int):
    graph_nodes = Ind.Ab + Ind.ABCDE + Ind.ABcDE + Ind.ABcdE + Ind.ABCDe + Ind.ABcDe + Ind.ABcde + Ind.ab + Ind.aBCDE + Ind.aBcDE + Ind.aBcdE + Ind.aBCDe + Ind.aBcDe + Ind.aBcde

# get stacked bar plot for pluripotent cell types
    # according to model param
    Ind['int_indep_nongraph'] = np.array(Int.CDE)+np.array(Int.CDe)
    Ind['notint_indep_nongraph'] = np.array(Int.cde) + np.array(Int.cdE) + np.array(Int.cDe) + np.array(Int.cDE)
    Ind['int_indep_graph'] = np.array(Ind.ABCDE) + np.array(Ind.ABCDe) + np.array(Ind.aBCDE) +  np.array(Ind.aBCDe)
    Ind['notint_indep_graph'] = np.array(Ind.Ab) + np.array(Ind.ABcDE) + np.array(Ind.ABcdE) + np.array(Ind.ABcDe) + np.array(Ind.ABcde) + np.array(Ind.ab) + np.array(Ind.aBcDE) + np.array(Ind.aBcdE) + np.array(Ind.aBcDe) + np.array(Ind.aBcde)

    top = ['notint_indep_nongraph','int_indep_nongraph','notint_indep_graph','int_indep_graph']
    X = range(1,6)
    Y = np.zeros((4,5))
    for x in [0,1,2,3]:
        for y in list(range(5)):
            Y[x,y] = np.sum(Ind[top[x]][np.array(Ind['N']) == y+3])
    for i in range(5):
        Y[:,i] = Y[:,i] / Y[:,i].sum()
    width = 0.75

    col = ['xkcd:maroon','xkcd:scarlet','xkcd:dark olive','xkcd:sage green']
    fig, ax = plt.subplots()
    bottom = np.zeros(5)
    for i in range(len(top)):
        ax.bar(X, Y[i], width, bottom=bottom, color=col[i], edgecolor='white', linewidth=0.5)
        bottom += Y[i]
    plt.xticks([1,2,3,4,5],('3','4','5','6','7'))
    plt.xlim((0.5,5.5))


# SUPPLEMENT FIGURES

def get_dag_loopnumber(Rep,indeg2):
    indeg_dag = [indeg2[i] for i in list(np.where(Rep.dag==1)[0])]
    df = Rep[Rep.dag==1]
    S = np.array([sum(indeg_dag[i]) for i in np.arange(len(indeg_dag))])
    L = np.array([len(indeg_dag[i]) for i in np.arange(len(indeg_dag))])
    loopfrac = (S-L+1)/S

    x1 = loopfrac[df.N == 3]
    x2 = loopfrac[df.N == 4]
    x3 = loopfrac[df.N == 5]
    x4 = loopfrac[df.N == 6]
    x5 = loopfrac[df.N == 7]

    bins = np.linspace(0,1,11)
    f1 = plt.hist([x1,x2,x3,x4,x5],bins,density=True,color=['xkcd:medium green','xkcd:denim blue','xkcd:gold','xkcd:pale red','black'],label=['N=3','4','5','6','7'])
    plt.legend(loc='best')

    plt.figure();sns.lineplot(df.asym,loopfrac)
    plt.figure();sns.lineplot(df.adj,loopfrac)
    plt.figure();sns.lineplot(df.sig,loopfrac) 
    plt.figure();sns.lineplot(df.graph_nodes,loopfrac)
    return loopfrac

def get_tree_convergence(indeg2,outdeg2,Rep):
    indeg_tree = [indeg2[i] for i in list(np.where(Rep.tree==1)[0])]
    outdeg_tree = [outdeg2[i] for i in list(np.where(Rep.tree==1)[0])]
    ne = [sum(indeg_tree[i]) for i in range(len(indeg_tree))]
    df = Rep[Rep.tree==1]    
 
    conv_in = np.array([0 if (sum(indeg_tree[i]>1)==0) else sum(indeg_tree[i][indeg_tree[i]>1]) for i in np.arange(len(indeg_tree))])
    div_out = np.array([0 if (sum(outdeg_tree[i]>1)==0) else sum(outdeg_tree[i][outdeg_tree[i]>1]) for i in np.arange(len(outdeg_tree))])
    div_deg = (div_out - conv_in)/ne
    # div_deg = (total outdegree of divergent nodes - total indegree of convergent nodes)/total number of edges 
    # for completely divergent trees, div_deg = 1. completely convergent trees, div_deg = -1. 
    x1 = div_deg[df.N == 3]
    x2 = div_deg[df.N == 4]
    x3 = div_deg[df.N == 5]
    x4 = div_deg[df.N == 6]
    x5 = div_deg[df.N == 7]
    
    f1 = plt.hist([x1,x2,x3,x4,x5],density=True,color=['xkcd:medium green','xkcd:denim blue','xkcd:gold','xkcd:pale red','black'])
    #plt.legend(loc='best')
    ylim = np.sum(f1[0][0])
    plt.yticks(np.linspace(0,ylim,3),('0','0.5','1'))
    
    plt.figure();sns.lineplot(df.asym,div_deg)
    plt.figure();sns.lineplot(df.adj,div_deg)
    plt.figure();sns.lineplot(df.sig,div_deg)
    plt.figure();sns.lineplot(df.graph_nodes,div_deg)

    return div_deg

#ef get_div_tree_levels(indeg2,outdeg2,Rep,div_deg):
#   indeg_tree = [indeg2[i] for i in list(np.where(Rep.tree==1)[0]) 
#   outdeg_tree = [outdeg2[i] for i in list(np.where(Rep.tree==1)[0])]    
#   
#   indeg_tree_div = [indeg_tree[i] for i in list(np.where(np.array(div_deg)==1)[0])]    
#   outdeg_tree_div = [outdeg_tree[i] for i in list(np.where(np.array(div_deg)==1)[0])]

#   levels = np.zeros((1,len(indeg_tree_div)))

#   for i in range(len(levels)):
#       in1 = indeg_tree_div[i]
#       out1 = outdeg_tree_div[i]
        


def non_random_edge_dist(Rep):
    x = [1,2,3,4,5,6,7,8,9,10]
    e1 = [Rep.num_edges[Rep.graph_nodes==i] for i in x]
    e1 = list(itertools.chain.from_iterable(e1))
    n1 = [Rep.graph_nodes[Rep.graph_nodes==i] for i in x]
    n1 = list(itertools.chain.from_iterable(n1))
    e2 = [Rep.num_edges[Rep.graph_nodes==2*i] for i in x]
    e2 = list(itertools.chain.from_iterable(e2))
    n2 = [Rep.graph_nodes[Rep.graph_nodes==2*i] for i in x]
    n2 = list(itertools.chain.from_iterable(n2))
    e2_exp = [4*e1[i] for i in range(len(e1))]#expected for erdos-renyi graphs
    
    sns.lineplot(n2,e2,color = 'k')                                                                                                               
    sns.lineplot(2*np.array(n1),e2_exp,color = 'grey')                                                                                            
    plt.xticks(np.arange(2,22,2),('2','4','6','8','10','12','14','16','18','20'))

