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

from reconciliation import solvePCR
from reconciliation import solveMCTPCR

def getAdjacencyMatrix(edgefile, nodes):
    adjMat = np.zeros((len(nodes), len(nodes)))
    with open(edgefile, 'r') as inp:
        for line in inp:
            data = line.rstrip('\n').split(',')
            node_out = data[0]
            node_in = data[1]

            if node_out == node_in:
                continue

            idx_out = nodes.index(node_out)
            idx_in = nodes.index(node_in)

            adjMat[idx_out, idx_in] = 1
            adjMat[idx_in, idx_out] = -1

    return adjMat

def getTreeEdges(edgefile, nodes):
    edgeList = []
    nodes = [str(node) for node in nodes]
    with open(edgefile, 'r') as inp:
        for line in inp:
            data = line.rstrip('\n').split(',')
            node_out = data[0]
            node_in = data[1]

            if node_out == node_in:
                continue

            idx_out = nodes.index(node_out)
            idx_in = nodes.index(node_in)

            edgeList.append((idx_out, idx_in))

    return edgeList


def main(args):

    if args.fsnv:
        df_fsnv = pd.read_csv(args.fsnv, sep=',', index_col = 'genotypes')

    if args.fcna:
        df_fcna = pd.read_csv(args.fcna, sep=',', index_col = 'genotypes')

    if not args.snv_tree and not args.cna_tree:

        solver = solvePCR(df_fsnv.values, df_fcna.values)
        solver.solve()
        solver.writeCloneFile(f'{args.o}_clone_prediction.out', snv_clones = list(df_fsnv.index), cna_clones = list(df_fcna.index))

    else:

        #solver.solveTree()
        snv_tree_edges = getTreeEdges(args.snv_tree, list(df_fsnv.index))
        cna_tree_edges = getTreeEdges(args.cna_tree, list(df_fcna.index))

        print('tree edges are ')
        print(getTreeEdges(args.snv_tree, list(df_fsnv.index)), getTreeEdges(args.cna_tree, list(df_fcna.index)), sep='\n')

        solver = solveMCTPCR(df_fsnv.values, df_fcna.values, snv_tree_edges, cna_tree_edges)
        solver.solve()

        #solver.solve()
        solver.writeCloneFile(f'{args.o}_clone_prediction.out', snv_clones = list(df_fsnv.index), cna_clones = list(df_fcna.index))
        solver.writeCloneTree(f'{args.o}_clone_tree_prediction.tsv', snv_clones = list(df_fsnv.index), cna_clones = list(df_fcna.index))
        solver.writeDOT(f'{args.o}_clone_tree.dot', snv_clones = list(df_fsnv.index), cna_clones = list(df_fcna.index))

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
    parser.add_argument('--fsnv', type=str, help='csv file with abundance of each SNV genotype in each sample', required=True)
    parser.add_argument('--fcna', type=str, help='csv file with abundance of each CNA genotype in each sample', required=True)
    parser.add_argument('--snv_tree', type=str, help='csv file containing the edges of SNV tree')
    parser.add_argument('--cna_tree', type=str, help='csv file containing the edges of CNA tree')
    parser.add_argument('-o', type=str, help='output prefix', required=True)

    args = parser.parse_args(None if sys.argv[1:] else ['-h'])

    main(args)
