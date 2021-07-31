#!/usr/bin/env python3

import sys
import numpy as np
import RNA
import matplotlib.pyplot as plt

from pyfaidx import Fasta
from scipy.spatial.distance import pdist, squareform
from sklearn.manifold import MDS
from sklearn.cluster import KMeans

def get_bppm(seq, p_res=1, icSHAPE=None, temp=37.0,
             p_method="api", p_display="full"):

    if p_method == "api":
        md = RNA.md()
        md.uniq_ML = 1
        md.temperature = temp
        fc = RNA.fold_compound(seq, md)
        if icSHAPE:
            m, b = 1.8, -0.6
            fc.sc_add_SHAPE_deigan(icSHAPE, m, b)
        fc.pf()
        structs = np.array(fc.bpp()) # upper triangle

    n_bins = len(seq)//p_res
    data = np.zeros((n_bins+1, n_bins+1))
    for i in range(len(structs)):
        for j in range(len(structs)):

            if p_display in ("full", "lower"):
                data[i//p_res,j//p_res] += structs[i,j]
            if p_display in ("full", "upper"):
                data[j//p_res,i//p_res] += structs[i,j]

    return data

def get_matrix(seq, p_res=1, constraints=None, icSHAPE=None,
               no_of_folds=250, subopt=True,
               p_method="api", p_display="full"):

    if p_method == "api":
        md = RNA.md()
        md.uniq_ML = 1
        md.temperature = 37.0
        fc = RNA.fold_compound(seq, md)

        if icSHAPE:
            m, b = 1.8, -0.6
            fc.sc_add_SHAPE_deigan(icSHAPE, m, b)

        if subopt:
            fc.pf()
            structs = [ fc.pbacktrack() for i in range(no_of_folds) ]
        else:
            # MFE
            struct, fe = fc.mfe()
            structs = [struct]

    n_bins = len(seq)//p_res
    data = np.zeros((n_bins+1, n_bins+1))
    for struct in structs:
        data += parse_struct_into_contactmap(struct, p_res, p_display)
    return data

def get_structs(seq, constraints=None, SHAPE=None,
                p_shape_type="deigan",
                no_of_folds=250, temp=37.0, subopt=True):
    md = RNA.md()
    md.uniq_ML = 1
    md.temperature = temp
    fc = RNA.fold_compound(seq, md)
    if SHAPE is not None:
        if p_shape_type == "deigan":
            #m, b = 1.8, -0.6
            m, b = 2.6, -0.8
            fc.sc_add_SHAPE_deigan(SHAPE, m, b)
        elif p_shape_type == "zarringhalam":
            beta = 0.89  # beta: 0.5 to 1.5 check paper fig 5
            fc.sc_add_SHAPE_zarringhalam(SHAPE, beta, -1, "Z")
    if subopt:
        fc.pf()
        structs = [ fc.pbacktrack() for i in range(no_of_folds) ]
    else:  # MFE
        struct, fe = fc.mfe()
        structs = [struct]
    return structs

def parse_struct_into_vector(struct, p_encode="double"):
    if p_encode == "double":
        dict_symbol = {".": 0, ",": 0, "(": 1, ")": 1, "{": 1, "}": 1, "[": 1, "]": 1}
    elif p_encode == "single":
        dict_symbol = {".": 1, ",": 1, "(": 0, ")": 0, "{": 0, "}": 0, "[": 0, "]": 0}
    new_struct = [ dict_symbol[s] for s in list(struct) ]
    return new_struct

def parse_struct_into_contactmap(struct, p_res, p_display):
    n_bins = len(struct)//p_res
    data = np.zeros((n_bins+1, n_bins+1))
    ds = []
    for j, fold in enumerate(struct):
        if fold == "(":
            ds.append(j)
        elif fold == ")":
            i = ds.pop()
            bin_i = i//p_res
            bin_j = j//p_res
            if p_display in ("full", "lower"):
                data[bin_i][bin_j] += 1
            if p_display in ("full", "upper"):
                data[bin_j][bin_i] += 1
    return data

if __name__ == "__main__":

    pass
