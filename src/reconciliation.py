#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 2021

@author: Palash Sashittal
"""

import gurobipy as gp
import numpy as np
import pandas as pd
import networkx as nx
import itertools

# parsimonious clone reconciliation problem
class solvePCR:

    def __init__(self, snv_mat, cna_mat, threads = 1, timelimit = None, verbose = True):
        self.snv_mat = snv_mat
        self.cna_mat = cna_mat
        self.threads = threads
        self.timelimit = timelimit
        self.verbose = verbose

        self.nsamples = self.snv_mat.shape[1]
        self.sol_clones = None
        self.sol_props = None

    def solve(self):

        nsamples = self.snv_mat.shape[1]
        assert nsamples == self.cna_mat.shape[1], 'SNV and CNA matrix sizes do not match up.'

        nsnv = self.snv_mat.shape[0]
        ncna = self.cna_mat.shape[0]

        model = gp.Model('solvePCR')
        x = model.addVars(nsnv, ncna, vtype=gp.GRB.BINARY, name='x')
        w = model.addVars(nsamples, nsnv, ncna, vtype = gp.GRB.CONTINUOUS, lb = 0, ub = 1, name = 'w')
        y = model.addVars(nsamples, nsnv, ncna, vtype = gp.GRB.CONTINUOUS, lb = 0, ub = 1, name = 'y')

        for i in range(nsamples):
            for j in range(nsnv):
                for k in range(ncna):
                    model.addConstr(w[i,j,k] <= y[i,j,k])
                    model.addConstr(w[i,j,k] <= x[j,k])
                    model.addConstr(w[i,j,k] >= x[j,k] + y[i,j,k] - 1)

        for i in range(nsamples):
            for j in range(nsnv):
                sum = gp.LinExpr()
                for k in range(ncna):
                    sum += w[i,j,k]
                model.addConstr(sum == self.snv_mat[j,i])

        for i in range(nsamples):
            for k in range(ncna):
                sum = gp.LinExpr()
                for j in range(nsnv):
                    sum += w[i,j,k]
                model.addConstr(sum == self.cna_mat[k,i])

        for i in range(nsamples):
            sum = gp.LinExpr()
            for j in range(nsnv):
                for k in range(ncna):
                    sum += w[i,j,k]
            model.addConstr(sum == 1)

        obj_sum = gp.LinExpr()
        for j in range(nsnv):
            for k in range(ncna):
                obj_sum += x[j,k]
        model.setObjective(obj_sum, gp.GRB.MINIMIZE)

        model.setParam(gp.GRB.Param.Threads, self.threads)
        model.optimize()

        if model.status == gp.GRB.OPTIMAL:
            solx = model.getAttr('x', x)
            self.sol_clones = [key for key, val in solx.items() if val >= 0.5]
            self.sol_props = model.getAttr('x', w)

    def writeCloneFile(self, clone_file, snv_clones = None, cna_clones = None):
        clone_data = []
        for clone in self.sol_clones:
            #for sample in range(self.nsamples):
            if snv_clones:
                snv_clone = snv_clones[clone[0]]
            else:
                snv_clone = clone[0]
            if cna_clones:
                cna_clone = cna_clones[clone[1]]
            else:
                cna_clone = clone[1]

            clone_data.append([clone, snv_clone, cna_clone] + [self.sol_props[sample, clone[0], clone[1]] for sample in range(self.nsamples)])
        df_clone = pd.DataFrame(clone_data, columns=['clone', 'snv_clone', 'cna_clone'] + [f'sample_{idx}' for idx in range(self.nsamples)])
        df_clone.to_csv(clone_file, sep='\t', index=False)

# minimum correction tree parsimonious clone reconciliation problem
class solveMCTPCR:

    def __init__(self, snv_mat, cna_mat, snv_edges, cna_edges, threads = 1, timelimit = None, verbose = True):
        self.snv_mat = snv_mat
        self.cna_mat = cna_mat

        G = nx.DiGraph()
        G.add_edges_from(snv_edges)
        self.snv_dag = nx.algorithms.transitive_closure_dag(G)
        G.clear()
        G.add_edges_from(cna_edges)
        self.cna_dag = nx.algorithms.transitive_closure_dag(G)

        self.threads = threads
        self.timelimit = timelimit
        self.verbose = verbose

        self.nsamples = self.snv_mat.shape[1]
        self.nsnv = self.snv_mat.shape[0]
        self.ncna = self.cna_mat.shape[0]

        self.sol_clones = None
        self.sol_props = None

        # snv/cna parent dictionary
        self.snv_parent_dict = {}
        for edge in snv_edges:
            child = edge[1]
            parent = edge[0]
            if child not in self.snv_parent_dict.keys():
                self.snv_parent_dict[child] = [parent]
            else:
                self.snv_parent_dict[child].append(parent)
        self.cna_parent_dict = {}
        for edge in cna_edges:
            child = edge[1]
            parent = edge[0]
            if child not in self.cna_parent_dict.keys():
                self.cna_parent_dict[child] = [parent]
            else:
                self.cna_parent_dict[child].append(parent)

        self.snv_edges = snv_edges
        self.cna_edges = cna_edges

        for j in range(self.nsnv):
            if j not in self.snv_parent_dict.keys():
                self.snv_root = j
                break
        for k in range(self.ncna):
            if k not in self.cna_parent_dict.keys():
                self.cna_root = k
                break

        print(f'snv root is {self.snv_root} and cna root is {self.cna_root}')

    def solve(self):

        nsamples = self.snv_mat.shape[1]
        assert nsamples == self.cna_mat.shape[1], 'SNV and CNA matrix sizes do not match up.'

        nsnv = self.snv_mat.shape[0]
        ncna = self.cna_mat.shape[0]

        model = gp.Model('solveMCTPCR')
        x = model.addVars(nsnv, ncna, vtype=gp.GRB.BINARY, name='x')
        w = model.addVars(nsamples, nsnv, ncna, vtype = gp.GRB.CONTINUOUS, lb = 0, ub = 1, name = 'w')
        y = model.addVars(nsamples, nsnv, ncna, vtype = gp.GRB.CONTINUOUS, lb = 0, ub = 1, name = 'y')
        z_snv = model.addVars(nsnv-1, ncna, vtype = gp.GRB.CONTINUOUS, lb = 0, ub = 1, name = 'z_snv')
        z_cna = model.addVars(nsnv, ncna-1, vtype = gp.GRB.CONTINUOUS, lb = 0, ub = 1, name = 'z_cna')
        d_snv = model.addVars(nsamples, nsnv, vtype = gp.GRB.CONTINUOUS, lb = 0, ub = 1, name = 'delta_snv')
        d_cna = model.addVars(nsamples, ncna, vtype = gp.GRB.CONTINUOUS, lb = 0, ub = 1, name = 'delta_cna')

        # encode product w[i,j,k] = y[i,j,k] * x[j,k]
        for i in range(nsamples):
            for j in range(nsnv):
                for k in range(ncna):
                    model.addConstr(w[i,j,k] <= y[i,j,k])
                    model.addConstr(w[i,j,k] <= x[j,k])
                    model.addConstr(w[i,j,k] >= x[j,k] + y[i,j,k] - 1)

        # encode abundance constraint for snv with correction
        for i in range(nsamples):
            for j in range(nsnv):
                sum = gp.LinExpr()
                for k in range(ncna):
                    sum += w[i,j,k]
                #model.addConstr(sum == self.snv_mat[j,i])
                model.addConstr(self.snv_mat[j,i] - sum <= d_snv[i,j])
                model.addConstr(sum - self.snv_mat[j,i] <= d_snv[i,j])

        # encode abundance constraint for cna
        for i in range(nsamples):
            for k in range(ncna):
                sum = gp.LinExpr()
                for j in range(nsnv):
                    sum += w[i,j,k]
                #model.addConstr(sum == self.cna_mat[k,i])
                model.addConstr(self.cna_mat[k,i] - sum <= d_cna[i,k])
                model.addConstr(sum - self.cna_mat[k,i] <= d_cna[i,k])

        # encode total abundance constraint
        for i in range(nsamples):
            sum = gp.LinExpr()
            for j in range(nsnv):
                for k in range(ncna):
                    sum += w[i,j,k]
            model.addConstr(sum == 1)

        # encode z_snv[j,k] = x[parent(j), k] * x[j,k]
        # encode z_cna[j,k] = x[j, parent(k)] * x[j,k]
        for edge_idx, edge in enumerate(self.snv_edges):
            parent = edge[0]
            child = edge[1]
            for k in range(ncna):
                model.addConstr(z_snv[edge_idx, k] <= x[parent, k])
                model.addConstr(z_snv[edge_idx, k] <= x[child, k])
                model.addConstr(z_snv[edge_idx, k] >= x[parent, k] + x[child, k] - 1)

        for edge_idx, edge in enumerate(self.cna_edges):
            parent = edge[0]
            child = edge[1]
            for j in range(nsnv):
                model.addConstr(z_cna[j, edge_idx] <= x[j, parent])
                model.addConstr(z_cna[j, edge_idx] <= x[j, child])
                model.addConstr(z_cna[j, edge_idx] >= x[j, parent] + x[j, child] - 1)

        # encode sum_{k} z_snv[j, k] == 1
        # encode sum_{j} z_cna[j, k] == 1
        for edge_idx in range(nsnv-1):
            sum = gp.LinExpr()
            for k in range(ncna):
                sum += z_snv[edge_idx, k]
            #model.addConstr(sum <= 1)
            model.addConstr(sum == 1)
        for edge_idx in range(ncna - 1):
            sum = gp.LinExpr()
            for j in range(nsnv):
                sum += z_cna[j, edge_idx]
            #model.addConstr(sum <= 1)
            model.addConstr(sum == 1)

        # encode x[j,k] <= x[parent(j), k] + x[j, parent(k)]
        for j in range(nsnv):
            for k in range(ncna):
                if j in self.snv_parent_dict.keys() or k in self.cna_parent_dict.keys():
                    sum = gp.LinExpr()
                    if j in self.snv_parent_dict.keys():
                        sum += x[self.snv_parent_dict[j][0], k]
                    if k in self.cna_parent_dict.keys():
                        sum += x[j, self.cna_parent_dict[k][0]]
                    model.addConstr(x[j,k] <= sum)

        # set objective function
        obj_sum = gp.LinExpr()
#        for j in range(nsnv):
#            for k in range(ncna):
#                obj_sum += x[j,k]
#        model.setObjective(obj_sum, gp.GRB.MINIMIZE)
        for i in range(nsamples):
            for j in range(nsnv):
                obj_sum += d_snv[i,j]
            for k in range(ncna):
                obj_sum += d_cna[i,k]
        model.setObjective(obj_sum, gp.GRB.MINIMIZE)

#        model.write('mctpcr.lp')

        model.setParam(gp.GRB.Param.Threads, self.threads)
        model.optimize()
        if model.status == gp.GRB.OPTIMAL:
            solx = model.getAttr('x', x)
            self.sol_clones = [key for key, val in solx.items() if val >= 0.5]
            self.sol_props = model.getAttr('x', w)


    def writeCloneFile(self, clone_file, snv_clones = None, cna_clones = None):
        clone_data = []
        for clone in self.sol_clones:
            #for sample in range(self.nsamples):
            if snv_clones:
                snv_clone = snv_clones[clone[0]]
            else:
                snv_clone = clone[0]
            if cna_clones:
                cna_clone = cna_clones[clone[1]]
            else:
                cna_clone = clone[1]

            clone_data.append([clone, snv_clone, cna_clone] + [self.sol_props[sample, clone[0], clone[1]] for sample in range(self.nsamples)])
        df_clone = pd.DataFrame(clone_data, columns=['clone', 'snv_clone', 'cna_clone'] + [f'sample_{idx}' for idx in range(self.nsamples)])
        df_clone.to_csv(clone_file, sep='\t', index=False)

    def writeCloneTree(self, clone_tree_file, snv_clones = None, cna_clones = None):
        clone_edges = []
        for clone1, clone2 in itertools.permutations(self.sol_clones, 2):
            if snv_clones:
                snv_clone1 = snv_clones[clone1[0]]
                snv_clone2 = snv_clones[clone2[0]]
            else:
                snv_clone1 = clone1[0]
                snv_clone2 = clone2[0]

            if cna_clones:
                cna_clone1 = cna_clones[clone1[1]]
                cna_clone2 = cna_clones[clone2[1]]
            else:
                cna_clone1 = clone1[1]
                cna_clone2 = clone2[1]

            if clone1[0] == clone2[0]:
                if clone1[1] in self.cna_parent_dict.keys():
                    if clone2[1] in self.cna_parent_dict[clone1[1]]:
                        clone_edges.append(((snv_clone2, cna_clone2), (snv_clone1, cna_clone1)))

            if clone1[1] == clone2[1]:
                if clone1[0] in self.snv_parent_dict.keys():
                    if clone2[0] in self.snv_parent_dict[clone1[0]]:
                        clone_edges.append(((snv_clone2, cna_clone2), (snv_clone1, cna_clone1)))

        with open(clone_tree_file, 'w') as output:
            for clone_edge in clone_edges:
                output.write(f'{clone_edge[0]}\t{clone_edge[1]}\n')

    def writeDOT(self, dot_file, snv_clones = None, cna_clones = None):

        with open(dot_file, 'w') as output:

            output.write(f'digraph N {{\n')
            output.write(f"\toverlap=\"false\"\n")
            output.write(f"\trankdir=\"TB\"\n")

            idx_dict = {}
            idx = 0
            for clone in self.sol_clones:
                if snv_clones:
                    snv_clone = snv_clones[clone[0]]
                else:
                    snv_clone = clone[0]
                if cna_clones:
                    cna_clone = cna_clones[clone[1]]
                else:
                    cna_clone = clone[1]

                idx_dict[clone] = idx
                output.write(f'\t{idx} [label=\"{snv_clone}, {cna_clone}\", style=\"bold\"];\n')

                idx += 1

            for clone1, clone2 in itertools.permutations(self.sol_clones, 2):
                if snv_clones:
                    snv_clone1 = snv_clones[clone1[0]]
                    snv_clone2 = snv_clones[clone2[0]]
                else:
                    snv_clone1 = clone1[0]
                    snv_clone2 = clone2[0]

                if cna_clones:
                    cna_clone1 = cna_clones[clone1[1]]
                    cna_clone2 = cna_clones[clone2[1]]
                else:
                    cna_clone1 = clone1[1]
                    cna_clone2 = clone2[1]

                if clone1[0] == clone2[0]:
                    if clone1[1] in self.cna_parent_dict.keys():
                        if clone2[1] in self.cna_parent_dict[clone1[1]]:
                            output.write(f"\t{idx_dict[clone2]} -> {idx_dict[clone1]} [style=\"bold\"];\n")

                if clone1[1] == clone2[1]:
                    if clone1[0] in self.snv_parent_dict.keys():
                        if clone2[0] in self.snv_parent_dict[clone1[0]]:
                            output.write(f"\t{idx_dict[clone2]} -> {idx_dict[clone1]} [style=\"bold\"];\n")

            output.write(f'}}')
