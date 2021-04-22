#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 2021

@author: Palash Sashittal
"""

import pandas as pd
import sys
import argparse
import itertools
import math
import numpy as np

import networkx as nx
import itertools

def main(args):

    np.random.seed(args.s)

    S = nx.DiGraph() # snv tree
    T = nx.DiGraph() # cna tree
    G = nx.DiGraph() # clone tree

    # add root nodes
    S.add_node(0)
    T.add_node(0)
    G.add_node((0,0))

    # build tree
    nmuts = args.m + args.d - 2
    mut_order = np.random.permutation(nmuts)
    snv_counter = 0
    cna_counter = 0
    for mut in mut_order:
        #print(f'mut is {mut}')
        # choose clone node at random
        parent_node = list(G.nodes)[np.random.randint(len(G.nodes))]
        if mut < args.m - 1:
            # SNV mutation
            snv_counter += 1
            G.add_edge(parent_node, (snv_counter, parent_node[1]))
            S.add_edge(parent_node[0], snv_counter)
            #print('added SNV mut')
            #print(f'{parent_node}, ({snv_counter}, {parent_node[1]})')
        else:
            # CNA mutation
            cna_counter += 1
            G.add_edge(parent_node, (parent_node[0], cna_counter))
            T.add_edge(parent_node[1], cna_counter)
            #print('added CNA mut')
            #print(f'{parent_node}, ({parent_node[0]}, {cna_counter})')

#    print(G.nodes)
#    print(G.edges)
#    print('-'*50)
#    print(S.nodes)
#    print(S.edges)
#    print('-'*50)
#    print(T.nodes)
#    print(T.edges)

    gamma = len(G.nodes)
    # get proportions of the clones
    clone_props = np.random.dirichlet([1]*gamma, args.n).transpose()

#    # construct noise
#    plus_noise = args.t * np.random.dirichlet([1]*gamma, args.n).transpose()
#    minus_noise = args.t * np.random.dirichlet([1]*gamma, args.n).transpose()
#    total_noise = plus_noise - minus_noise
#
#    noisy_clone_props = clone_props + total_noise

    # build clone dataframe
    data_clone = []
    for idx, clone in enumerate(list(G.nodes)):
        data_clone.append([clone, clone[0], clone[1]] + list(clone_props[idx, :]))

    df_clone = pd.DataFrame(data_clone, columns=['clone', 'snv', 'cna'] + [f'sample_{idx}' for idx in range(args.n)])

#    # build noisy clone dataframe
#    data_noisy_clone = []
#    for idx, clone in enumerate(list(G.nodes)):
#        data_noisy_clone.append([clone, clone[0], clone[1]] + list(noisy_clone_props[idx, :]))
#
#    df_noisy_clone = pd.DataFrame(data_noisy_clone, columns=['clone', 'snv', 'cna'] + [f'sample_{idx}' for idx in range(args.n)])
#

    # write clone tree and dataframe
    df_clone.to_csv(f'{args.o}_clone.out', sep='\t', index=False)
    nx.write_edgelist(G, f'{args.o}_clone_tree.tsv', data=False, delimiter='\t')

#    # write noisy clone dataframe
#    df_noisy_clone.to_csv(f'{args.o}_clone_noisy.out', sep='\t', index=False)

    df_snv = df_clone.groupby('snv').sum().drop('cna', axis=1)
    df_cna = df_clone.groupby('cna').sum().drop('snv', axis=1)
    df_snv.index.names = ['genotypes']
    df_cna.index.names = ['genotypes']

    # write SNV and CNA tree and dataframe
    df_snv.to_csv(f'{args.o}_snv.csv')
    df_cna.to_csv(f'{args.o}_cna.csv')
    nx.write_edgelist(S, f'{args.o}_snv_tree.csv', data=False, delimiter=',')
    nx.write_edgelist(T, f'{args.o}_cna_tree.csv', data=False, delimiter=',')

    # get the transitive closures of S and T and write those edges
    S_dag = nx.transitive_closure_dag(S)
    T_dag = nx.transitive_closure_dag(T)
    nx.write_edgelist(S_dag, f'{args.o}_snv_dag.csv', data=False, delimiter=',')
    nx.write_edgelist(T_dag, f'{args.o}_cna_dag.csv', data=False, delimiter=',')

    # construct noise for SNV
    noisy_values = (1 - args.t) * df_snv.values + args.t * np.random.dirichlet([1]*args.m, args.n).transpose()
    df_snv_noisy = pd.DataFrame(noisy_values, index=df_snv.index, columns=df_snv.columns)

    #base_noise = args.t * np.random.dirichlet([1]*args.d, args.n).transpose()
    noisy_values = (1 - args.t) * df_cna.values + args.t * np.random.dirichlet([1]*args.d, args.n).transpose()
    df_cna_noisy = pd.DataFrame(noisy_values, index=df_cna.index, columns=df_cna.columns)

    df_snv_noisy.to_csv(f'{args.o}_snv_noisy.csv')
    df_cna_noisy.to_csv(f'{args.o}_cna_noisy.csv')

def str2bool(v):
    if isinstance(v, bool):
       return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', type=int, help='number of samples [1]', default = 1)
    parser.add_argument('-m', type=int, help='number of SNV genotypes [5]', default = 5)
    parser.add_argument('-d', type=int, help='number of CNA genotypes [4]', default = 4)
    parser.add_argument('-o', type=str, help='output prefix', default='sample')
    parser.add_argument('-s', type=int, help='seed [0]', default = 0)
    parser.add_argument('-t', type=float, help='noise threshold [0.05]', default = 0.05)
    args = parser.parse_args(None if sys.argv[1:] else ['-h'])

    main(args)
