#!/usr/bin/env python3

import matplotlib

import sys
import gzip
import argparse
import pypairix
import numpy as np
import scipy.stats as spstats

import pandas as pd
import pyBigWig
from pyfaidx import Fasta
import matplotlib.pyplot as plt
from .utils_suboptimal_structure import get_matrix, get_bppm

proj_dir = ""

def parse_structure(f_structure, p_res, p_display):
    struct = ""
    with open(f_structure, "r") as f:
        for line in f:
            struct += line.strip("\r\n")
    n_bins = len(struct)//p_res
    data = np.zeros((n_bins+1, n_bins+1))
    #print(struct.count("("), struct.count(")"), struct.count("."))
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

def get_size(f_size):
    sizes = {}
    with open(f_size, "r") as f:
        for line in f:
            row = line.strip("\r\n").split("\t")
            sizes[row[0]] = int(row[1])
    return sizes

def get_annotation(f_annotation):
    annotations = {}
    with open(f_annotation, "r") as f:
        for line in f:
            row = line.strip("\r\n").split("\t")
            if int(row[2]) == int(row[3]):  # entire CDS is UTR
                annotations[row[0]] = (int(row[1]), int(row[4]))
            else:
                annotations[row[0]] = (int(row[2]), int(row[3]))
    return annotations

def get_splice(f_bed12):
    splices = {}
    with open(f_bed12, "r") as f:
        for line in f:
            row = line.strip("\r\n").split("\t")
            iid = row[3]
            strand = row[5]
            exon_blocks = np.array(list(map(int, row[10].split(",")[:-1])))
            if strand == "+":
                splice_sites = exon_blocks.cumsum()
            elif strand == "-":
                splice_sites = exon_blocks[::-1].cumsum()

            splices[iid] = splice_sites[:-1]
    return splices

def get_longest_by_gene(f_index, gene):
    index = {}
    with open(f_index, "r") as f:
        for line in f:
            row = line.strip("\r\n").split("\t")
            try:
                index[row[1]].append(row)
            except KeyError:
                index[row[1]] = [row]

    try:
        query = index[gene][0]  # longest
        rid, start, end = (query[0], 0, int(query[2])-1)
    except KeyError:
        rid, start, end = (None, None, None)

    return rid, start, end

def get_by_rid(f_index, rid):
    index = {}
    with open(f_index, "r") as f:
        for line in f:
            row = line.strip("\r\n").split("\t")
            index[row[0]] = row

    try:
        query = index[rid]
        start, end = (0, int(query[2])-1)
    except KeyError:
        start, end = (None, None)

    return start, end

def load_loops(chrom, p_start, p_end, f_loop_path,
               p_resolution, p_point=False):

    loop = []
    with open(f_loop_path, "r") as f:
        for line in f:
            if line.startswith("track"):
                continue
            row = line.strip("\r\n").split("\t")
            #lchrom, lstart1, lend1, lstart2, lend2 = (row[0],
            #                                          int(row[1]), int(row[1])+10,
            #                                          int(row[2]), int(row[2])+10)

            if len(row) > 6:
                (lchrom1, lstart1, lend1,
                 lchrom2, lstart2, lend2, category) = (row[0], int(row[1]), int(row[2]),
                                                       row[3], int(row[4]), int(row[5]),
                                                       row[6])
            else:
                (lchrom1, lstart1, lend1,
                 lchrom2, lstart2, lend2) = (row[0], int(row[1]), int(row[2]),
                                             row[3], int(row[4]), int(row[5]))

            # TODO intra: check whether they are in the same chromosome
            # TODO inter: no checks
            #if lchrom != chrom:
            #    continue

            #if lstart1 >= p_start and lend2 <= p_end:
            if len(row) > 6:
                if p_point:
                    p_pres = p_resolution * 2
                    loop.append((lstart1-p_pres, lstart1+p_pres,
                                 lstart2-p_pres, lstart2+p_pres, category))
                else:
                    loop.append((lstart1, lend1, lstart2, lend2, category))
            else:
                if p_point:
                    p_pres = p_resolution * 2
                    loop.append((lstart1-p_pres, lstart1+p_pres,
                                 lstart2-p_pres, lstart2+p_pres))
                else:
                    loop.append((lstart1, lend1, lstart2, lend2))
    return loop

def tabulate(data, f_pair, loci, p_start1, p_start2,
             p_end1, p_end2,
             p_chrom1, p_chrom2, p_resolution,
             p_downsample=1.0):

    np.random.seed(386)

    for no, locus in enumerate(loci):

        pairs = pypairix.open(f_pair)
        iterator = pairs.querys2D(locus)

        for row in iterator:

            if np.random.random() > p_downsample:
                continue

            chrom1, start1, chrom2, start2 = (row[1], int(row[2]),
                                              row[3], int(row[4]))

            if p_start1 and p_end1 and p_start2 and p_end2:
                binth1 = (start1-p_start1)//p_resolution
                binth2 = (start2-p_start2)//p_resolution

                data[binth1][binth2] += 1

            elif p_start1 and p_end1:
                binth1 = (start1-p_start1)//p_resolution
                binth2 = (start2-p_start1)//p_resolution

                data[binth1][binth2] += 1
                data[binth2][binth1] += 1

            elif p_chrom1 and p_chrom2:  # flip

                binth1 = start1//p_resolution
                binth2 = start2//p_resolution

                if no == 0:
                    data[binth1][binth2] += 1
                elif no == 1:
                    data[binth2][binth1] += 1

            elif p_chrom1:
                binth1 = start1//p_resolution
                binth2 = start2//p_resolution

                data[binth1][binth2] += 1
                data[binth2][binth1] += 1

    return data

def normalize_by_diag(X):
    pseudocount = 1e-10
    N = np.zeros(X.shape)
    for i in range(X.shape[0]):
        pos_indices = np.where(np.eye(X.shape[0], k=i)==1)
        neg_indices = np.where(np.eye(X.shape[0], k=-i)==1)
        N[pos_indices] = np.log2((X[pos_indices] + pseudocount) / (np.diagonal(X, i).mean() + pseudocount))
        N[neg_indices] = np.log2((X[neg_indices] + pseudocount) / (np.diagonal(X, -i).mean() + pseudocount))
    return N

def normalize_diff(X0, X1):
    pseudocount = 1e-10
    Xm = (X0 + X1)/2
    Xd = X0 - X1
    print(Xm.shape)
    print(Xd.shape)
    N = np.zeros(Xm.shape)
    for i in range(Xm.shape[0]):
        pos_indices = np.where(np.eye(Xm.shape[0], k=i)==1)
        neg_indices = np.where(np.eye(Xm.shape[0], k=-i)==1)
        N[pos_indices] = (Xd[pos_indices] + pseudocount) / (np.diagonal(Xm, i).mean() + pseudocount)
        N[neg_indices] = (Xd[neg_indices] + pseudocount) / (np.diagonal(Xm, -i).mean() + pseudocount)
    return N

def get_ref(p_genome):
    dict_sizes = {"mm10": proj_dir + "data/reference/mm10/mm10_refseq.sizes",
                  "hg19": proj_dir + "data/reference/hg19/hg19_refseq_wo_version.sizes"}
    dict_fasta = {"mm10": proj_dir + "data/reference/mm10/mm10.fa",
                  "hg19": proj_dir + "data/reference/hg19/hg19_refseq_wo_version.fasta"}
    try:
        f_size = dict_sizes[p_genome]
        f_fasta = dict_fasta[p_genome]
    except KeyError:
        f_size = None
        f_fasta = None
    return f_size, f_fasta

def get_loci(p_chrom1=None, p_chrom2=None,
             p_start1=None, p_start2=None,
             p_end1=None, p_end2=None, p_gene1=None, p_gene2=None,
             p_resolution=100, p_genome=None, f_size=None, f_index=None, **kwargs):

    if p_genome:
        f_size, f_fasta = get_ref(p_genome)

    sizes = get_size(f_size)

    bins1 = None
    bins2 = None

    if p_chrom1 and (p_start1 is None) and (p_end1 is None):
        p_start1, p_end1 = 0, sizes[p_chrom1]-1
    if p_chrom2 and (p_start2 is None) and (p_end2 is None):
        p_start2, p_end2 = 0, sizes[p_chrom2]-1
    
    if (p_start1 is not None) and (p_end1 is not None) and (p_start2 is not None) and (p_end2 is not None):
        bins1 = (p_end1 - p_start1)//p_resolution
        bins2 = (p_end2 - p_start2)//p_resolution
        locus = "%s:%s-%s|%s:%s-%s" % (p_chrom1, p_start1, p_end1, p_chrom2, p_start2, p_end2)
        loci = [locus]
        m_bins = bins1+1
        n_bins = bins2+1
        #DEBUG print("A")
    elif (p_start1 is not None) and (p_end1 is not None):
        p_start1, p_end1 = (int(p_start1), int(p_end1))
        bins1 = (p_end1 - p_start1)//p_resolution
        bins2 = None
        locus = "%s:%s-%s|%s:%s-%s" % (p_chrom1, p_start1, p_end1, p_chrom1, p_start1, p_end1)
        loci = [locus]
        m_bins = bins1+1
        n_bins = bins1+1
        #DEBUG print("B")
    elif (p_chrom1 is not None) and (p_chrom2 is not None):
        bins1 = sizes[p_chrom1]//p_resolution
        bins2 = sizes[p_chrom2]//p_resolution
        locus1 = "%s|%s" % (p_chrom1, p_chrom2)
        locus2 = "%s|%s" % (p_chrom2, p_chrom1)
        loci = [locus1, locus2]
        m_bins = bins1+1
        n_bins = bins2+1
        #DEBUG print("C")

    return (loci, p_chrom1, p_start1, p_end1,
            p_chrom2, p_start2, p_end2,
            m_bins, n_bins, bins1, bins2)

def get_contactmap(f_pairs, p_chrom1=None, p_chrom2=None,
                   p_start1=None, p_start2=None,
                   p_end1=None, p_end2=None, p_gene1=None, p_gene2=None,
                   p_resolution=100, p_multi=False, p_diff=False,
                   p_genome=None, f_size=None, f_index=None, p_downsample=1.0,
                   p_point=False, **kwargs):

    (loci, p_chrom1, p_start1, p_end1,
     p_chrom2, p_start2, p_end2,
     m_bins, n_bins, bins1, bins2) = get_loci(p_chrom1=p_chrom1, p_chrom2=p_chrom2,
                                              p_start1=p_start1, p_start2=p_start2,
                                              p_end1=p_end1, p_end2=p_end2, p_gene1=p_gene1,
                                              p_gene2=p_gene2,
                                              p_resolution=p_resolution, p_genome=p_genome,
                                              f_size=f_size, f_index=f_index)

    raw = []
    if p_multi:

        for f_pair in f_pairs:
            for f_split in f_pair.split(","):
                datum = np.zeros((m_bins, n_bins))
                datum = tabulate(datum, f_split, loci,
                                 p_start1, p_start2,
                                 p_end1, p_end2,
                                 p_chrom1, p_chrom2, p_resolution,
                                 p_downsample)
            raw.append(datum)

    elif p_diff:

        if len(f_pairs) != 2:
            assert False

        for f_pair in f_pairs:
            for f_split in f_pair.split(","):
                datum = np.zeros((m_bins, n_bins))
                datum = tabulate(datum, f_split, loci,
                                 p_start1, p_start2,
                                 p_end1, p_end2,
                                 p_chrom1, p_chrom2, p_resolution,
                                 p_downsample)
            raw.append(datum)

    else:  # accumulate

        datum = np.zeros((m_bins, n_bins))
        for f_pair in f_pairs:
            for f_split in f_pair.split(","):
                datum = tabulate(datum, f_split, loci,
                                 p_start1, p_start2,
                                 p_end1, p_end2,
                                 p_chrom1, p_chrom2, p_resolution,
                                 p_downsample)

        raw.append(datum)

    return raw

def normalize_contactmap(raw, p_mode="none"):
    norm = []
    for datum in raw:
        if p_mode == "self-relative":
            depth = datum.sum()//2
            datum = datum/float(depth)
            print("relative using self depth")
        elif p_mode == "KR":
            import scipy.sparse as sps
            from .jupyter_KR import knightRuizAlg, removeZeroDiagonalCSR, addZeroes

            #test = np.random.randint(0,20+1, size=(5,5))
            #test = (test + test.T)/2
            #test = sps.csr_matrix(test)
            #print(test.toarray())
            #DEBUG print(datum.sum(axis=0))
            sdatum = sps.csr_matrix(datum) #+1e-5)

            # remove sparse or diagonals with zeros
            percentOfSparseToRemove = 0
            mtxAndRemoved = removeZeroDiagonalCSR(sdatum, percentOfSparseToRemove)
            initialSize = sdatum.shape[0]
            sdatum, removed = mtxAndRemoved[0], mtxAndRemoved[1]
            newSize = sdatum.shape[0]

            #DEBUG print(sdatum.toarray().sum(axis=0))
            #DEBUG print(removed)

            CC = ((sdatum.sum())/(sdatum.shape[0]*2))
            CCother = sdatum.sum()/sdatum.size
            result = knightRuizAlg(sdatum)
            col = result[0]
            x = sps.diags(col.flatten(), 0, format='csr')
            mtx = x.dot(sdatum.dot(x))

            #print(mtx.toarray())
            CCmtx = CC * mtx
            CCothermtx = CCother * mtx

            # restore zeros
            datum = addZeroes(CCmtx.toarray(), removed)
            #datum = CCmtx.toarray()
            #datum = CCothermtx.toarray()

            print("balanced")
        elif p_mode == "VC":
            rowsum = datum.sum(axis=0)
            datum = datum / np.sqrt(np.outer(rowsum, rowsum))
            datum = np.nan_to_num(datum)
            print("VC")
        elif p_mode == "normalize":
            #depth = datum.sum()/2
            #datum = datum/float(depth)
            datum = normalize_by_diag(datum)
        norm.append(datum)
    return norm

def compute_differential(norm, p_diff="norm_diff"):
    pseudocount = 1e-10
    if p_diff == "log2FC":
        results = np.log2((norm[0]+pseudocount)/(norm[1]+pseudocount))
    elif p_diff == "diff":
        results = (norm[0]+pseudocount)-(norm[1]+pseudocount)
    elif p_diff == "pct_diff":
        results = ((norm[0]+pseudocount)-(norm[1]+pseudocount))/(norm[1]+pseudocount)
    elif p_diff == "norm_diff":
        results = normalize_diff(norm[0], norm[1])
    return results

def set_n_max(norm, p_max, p_mode):
    if p_max is not None:
        n_max = p_max
    else:
        if p_mode == "relative":
            n_max = max([datum.mean()*3 for datumn in norm])
        elif p_mode == "VC":
            n_max = 0.3
        else:
            n_max = 1
    return n_max

def set_p_mode(p_relative, p_self_relative, p_balance, p_vc, p_normalize, **kwargs):
    p_mode = kwargs["p_mode"]
    if p_relative:
        p_mode = "relative"
    if p_self_relative:
        p_mode = "self-relative"
    if p_balance:
        p_mode = "KR"
    if p_vc:
        p_mode = "VC"
    if p_normalize:
        p_mode = "normalize"
    return p_mode

def set_dims(f_pairs, p_col, p_multi, **kwargs):
    if p_multi:
        n_rows = ((len(f_pairs)-1)//p_col)+1
        n_cols = ((len(f_pairs)-1)%p_col)+1
    else:
        n_rows = 1
        n_cols = 1
    return n_rows, n_cols

def add_structure_overlay(ax,
                          f_dot, p_resolution, p_display, p_structure, f_fasta,
                          p_chrom1, p_start1, p_end1, p_dot_size, p_col, **kwargs):
        
    datum_structs = []
    color_structs = ["#256baf", "red"]
    #color_structs = ["purple", "red"]

    if f_dot:
        datum_structs.append(parse_structure(f_dot, p_resolution, p_display))
    if p_structure:
        o_fasta = Fasta(f_fasta)
        seq = str(o_fasta[p_chrom1][p_start1:p_end1]).upper()
        datum_structs.append(get_bppm(seq, p_res=p_resolution, p_display=p_display))

    for datum_struct, color in zip(datum_structs, color_structs):
        X = []
        Y = []
        size = []
        for i in range(datum_struct.shape[0]-1):
            for j in range(datum_struct.shape[1]-1):
                x = datum_struct[i,j]
                X.append(i)
                Y.append(j)
                size.append(x)
        size = np.array(size)
        size = size * p_dot_size

        # filter out super small points in the overlay
        size[size<0.001] = 0.0
        
        for i in range(ax.shape[0]*ax.shape[1]):
            r, c = i//p_col, i%p_col
            ax[r][c].scatter(X, Y, marker="o", s=size, c=color, alpha=1.0)
        
    return ax

def add_loop_overlay(ax,
                     p_chrom1, p_start1, p_end1,
                     f_loop, p_resolution, p_point, **kwargs):
    dict_color = {"loops": "black", "lstripe": "blue", "rstripe": "green"}
    if f_loop is not None:
        p_boxwidth = 1.0

        loop = load_loops(p_chrom1, p_start1, p_end1, f_loop,
                          p_resolution, p_point=p_point)

        for l in loop:
            if len(l) == 5:
                lstart1, lend1, lstart2, lend2, category = l
                color = dict_color.get(category, "orange")
            else:
                lstart1, lend1, lstart2, lend2 = l
                color = "black" #"blue"

            loop_start1 = ((lstart1-p_start1)//p_resolution)-0.5
            loop_end1 = ((lend1-p_start1)//p_resolution)-0.5
            loop_start2 = ((lstart2-p_start1)//p_resolution)-0.5
            loop_end2 = ((lend2-p_start1)//p_resolution)-0.5

            def define_boxes(ax):
                ax.hlines(loop_start1-p_boxwidth, loop_start2-p_boxwidth,
                          loop_end2+p_boxwidth, color=color, alpha=1.0, linewidth=1.0)
                ax.hlines(loop_end1+p_boxwidth, loop_start2-p_boxwidth,
                          loop_end2+p_boxwidth, color=color, alpha=1.0, linewidth=1.0)

                ax.vlines(loop_start2-p_boxwidth, loop_start1-p_boxwidth,
                          loop_end1+p_boxwidth, color=color, alpha=1.0, linewidth=1.0)
                ax.vlines(loop_end2+p_boxwidth, loop_start1-p_boxwidth,
                          loop_end1+p_boxwidth, color=color, alpha=1.0, linewidth=1.0)
                return ax
            
            for i in range(ax.shape[0]*ax.shape[1]):
                r, c = i//p_col, i%p_col
                ax[r][c] = define_boxes(ax[r][c])
            
    return ax

def label_axes(ax, sizes, p_multi, p_label_res,
               p_chrom1, p_start1, p_end1,
               p_chrom2, p_start2, p_end2,
               bins1, bins2):
    
    xloc = ax[0][0].get_xticks()
    yloc = ax[0][0].get_yticks()
    xlabel = ax[0][0].get_xticklabels()
    ylabel = ax[0][0].get_yticklabels()

    if p_start1 and p_end1 and p_start2 and p_end2:
        xnewloc = np.linspace(0, bins2, len(xloc)-1)
        ynewloc = np.linspace(0, bins1, len(yloc)-1)
        xnewlabel = [ "%.3f" % (i/p_label_res) for i in np.linspace(p_start2, p_end2, len(xloc)-1) ]
        ynewlabel = [ "%.3f" % (i/p_label_res) for i in np.linspace(p_start1, p_end1, len(yloc)-1) ]
    elif p_start1 and p_end1:
        xnewloc = np.linspace(0, bins1, len(xloc)-1)
        ynewloc = np.linspace(0, bins1, len(yloc)-1)
        xnewlabel = [ "%.1f" % (i/p_label_res) for i in np.linspace(p_start1, p_end1, len(xloc)-1) ]
        ynewlabel = [ "%.1f" % (i/p_label_res) for i in np.linspace(p_start1, p_end1, len(yloc)-1) ]
    elif p_chrom1 and p_chrom2:
        xnewloc = np.linspace(0, bins2, len(xloc)-1)
        ynewloc = np.linspace(0, bins1, len(yloc)-1)
        xnewlabel = [ "%.1f" % (i/p_label_res) for i in np.linspace(0, sizes[p_chrom2], len(xloc)-1) ]
        ynewlabel = [ "%.1f" % (i/p_label_res) for i in np.linspace(0, sizes[p_chrom1], len(yloc)-1) ]
    elif p_chrom1:
        xnewloc = np.linspace(0, bins1, len(xloc)-1)
        ynewloc = np.linspace(0, bins1, len(yloc)-1)
        xnewlabel = [ "%.1f" % (i/p_label_res) for i in np.linspace(0, sizes[p_chrom1], len(xloc)-1) ]
        ynewlabel = [ "%.1f" % (i/p_label_res) for i in np.linspace(0, sizes[p_chrom1], len(yloc)-1) ]

    plt.setp(ax, xticks=xnewloc, yticks=ynewloc,
             xticklabels=xnewlabel, yticklabels=ynewlabel)

    return ax

def create_axes(n_rows, n_cols, p_ax, p_figsize, p_dpi, **kwargs):
    if p_ax is None:
        fig, ax = plt.subplots(n_rows, n_cols,
                               figsize=((p_figsize[0]*n_cols)+2, p_figsize[1]*n_rows),
                               dpi=p_dpi, squeeze=0)
    else:
        fig, ax = p_ax
        ax = np.array([[ax]])  # squeeze
    return fig, ax


def plot_contactmap(f_pairs, p_chrom1=None, p_chrom2=None,
                    p_start1=None, p_start2=None,
                    p_end1=None, p_end2=None, p_gene1=None, p_gene2=None,
                    p_resolution=100, p_max=None,
                    p_balance=False, p_relative=False, p_self_relative=False,
                    p_normalize=False, p_mode="none",
                    p_vc=False, p_multi=False, p_col=3,
                    p_genome=None, p_structure=False,
                    p_display="full",
                    f_size=None, f_out=None,
                    f_loop=None, f_dot=None, f_fasta=None,
                    f_index=None,
                    p_downsample=1.0, p_figsize=(4, 4), p_dot_size=5.0, p_point=False,
                    p_label_res=1000, p_dpi=100, p_ax=None, p_verbose=True):
    
    if p_genome:
        f_size, f_fasta = get_ref(p_genome)
    
    dict_kwargs = {"f_pairs": f_pairs, "p_chrom1": p_chrom1, "p_chrom2": p_chrom2,
                   "p_start1": p_start1, "p_start2": p_start2,
                   "p_end1": p_end1, "p_end2": p_end2,
                   "p_gene1": p_gene1, "p_gene2": p_gene2,
                   "p_resolution": p_resolution, "p_max": p_max,
                   "p_balance": p_balance, "p_relative": p_relative, "p_self_relative": p_self_relative,
                   "p_normalize": p_normalize, "p_mode": p_mode,
                   "p_vc": p_vc, "p_multi": p_multi, "p_col": p_col, "p_genome": p_genome, "p_structure": p_structure,
                   "p_display": p_display, "f_size": f_size, "f_out": f_out, "f_loop": f_loop,
                   "f_dot": f_dot, "f_fasta": f_fasta, "f_index": f_index,
                   "p_downsample": p_downsample, "p_figsize": p_figsize,
                   "p_dot_size": p_dot_size, "p_point": p_point, "p_label_res": p_label_res,
                   "p_dpi": p_dpi, "p_ax": p_ax, "p_verbose": p_verbose}

    (loci, p_chrom1, p_start1, p_end1,
     p_chrom2, p_start2, p_end2,
     m_bins, n_bins, bins1, bins2) = get_loci(**dict_kwargs)
    
    sizes = get_size(f_size)    
    raw = get_contactmap(**dict_kwargs)
    p_mode = set_p_mode(**dict_kwargs)
    norm = normalize_contactmap(raw, p_mode=p_mode)
    n_max = set_n_max(norm, p_max, p_mode)
    n_rows, n_cols = set_dims(**dict_kwargs)
    fig, ax = create_axes(n_rows, n_cols, **dict_kwargs)
    
    if p_verbose:
        print(ax.shape)
        print(loci)        
        for datum in raw:
            print("Dim: %s | Total Counts: %d" % (datum.shape, datum.sum()//2))
        print(n_rows, n_cols)

    # ================
    # Plot Contactmap
    # ================
    for i in range(ax.shape[0]*ax.shape[1]):
        r, c = i//p_col, i%p_col
        indices = (r*p_col)+c
        if p_normalize:
            im = ax[r][c].imshow(norm[indices], cmap="coolwarm",
                                 vmin=-abs(float(n_max)), vmax=abs(float(n_max)))
        else:
            im = ax[r][c].imshow(norm[indices], cmap = "Reds", vmin=0, vmax=float(n_max))

        for tick in ax[r][c].get_xticklabels():
            tick.set_rotation(45)

    # ============================
    # Adjust fig and add colorbar
    # ============================
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.02, 0.7])
    fig.colorbar(im, cax=cbar_ax)
        
    # =======================
    # Label axes and ticks
    # =======================
    ax = label_axes(ax, sizes, p_multi, p_label_res,
                    p_chrom1, p_start1, p_end1,
                    p_chrom2, p_start2, p_end2,
                    bins1, bins2)
    
    # ======================================
    # Set up structure overlay onto heatmap
    # ======================================
    ax = add_structure_overlay(ax, **dict_kwargs)

    # ============
    # Plot Loops
    # ============
    ax = add_loop_overlay(ax, **dict_kwargs)

    if f_out:
        fig.savefig(f_out)

    return fig, ax


def plot_diffmap(f_pairs, p_chrom1=None, p_chrom2=None,
                 p_start1=None, p_start2=None,
                 p_end1=None, p_end2=None, p_gene1=None, p_gene2=None,
                 p_resolution=100, p_max=None,
                 p_balance=False, p_relative=False, p_self_relative=False,
                 p_normalize=False, p_mode="none",
                 p_vc=False, p_multi=False, p_col=3,
                 p_diff="norm_diff", p_genome=None, p_structure=False,
                 p_display="full",
                 f_size=None, f_out=None,
                 f_loop=None, f_dot=None, f_fasta=None,
                 f_index=None,
                 p_downsample=1.0, p_figsize=(4, 4), p_dot_size=5.0, p_point=False,
                 p_label_res=1000, p_dpi=100, p_ax=None, p_verbose=True):
    
    if p_genome:
        f_size, f_fasta = get_ref(p_genome)
    
    dict_kwargs = {"f_pairs": f_pairs, "p_chrom1": p_chrom1, "p_chrom2": p_chrom2,
                   "p_start1": p_start1, "p_start2": p_start2,
                   "p_end1": p_end1, "p_end2": p_end2,
                   "p_gene1": p_gene1, "p_gene2": p_gene2,
                   "p_resolution": p_resolution, "p_max": p_max,
                   "p_balance": p_balance, "p_relative": p_relative, "p_self_relative": p_self_relative,
                   "p_normalize": p_normalize, "p_mode": p_mode,
                   "p_vc": p_vc, "p_multi": p_multi, "p_col": p_col, "p_diff": p_diff, "p_genome": p_genome, "p_structure": p_structure,
                   "p_display": p_display, "f_size": f_size, "f_out": f_out, "f_loop": f_loop,
                   "f_dot": f_dot, "f_fasta": f_fasta, "f_index": f_index,
                   "p_downsample": p_downsample, "p_figsize": p_figsize,
                   "p_dot_size": p_dot_size, "p_point": p_point, "p_label_res": p_label_res,
                   "p_dpi": p_dpi, "p_ax": p_ax, "p_verbose": p_verbose}

    (loci, p_chrom1, p_start1, p_end1,
     p_chrom2, p_start2, p_end2,
     m_bins, n_bins, bins1, bins2) = get_loci(**dict_kwargs)
    
    sizes = get_size(f_size)    
    raw = get_contactmap(**dict_kwargs)
    p_mode = set_p_mode(**dict_kwargs)
    norm = normalize_contactmap(raw, p_mode=p_mode)
    n_max = set_n_max(norm, p_max, p_mode)
    n_rows, n_cols = set_dims(**dict_kwargs)
    fig, ax = create_axes(n_rows, n_cols, **dict_kwargs)
    
    if p_verbose:
        print(loci)        
        for datum in raw:
            print("Dim: %s | Total Counts: %d" % (datum.shape, datum.sum()//2))
        print(n_rows, n_cols)

    # =================
    # Differential map
    # =================
    results = compute_differential(norm, p_diff=p_diff)
    im = ax[0][0].imshow(results, cmap="coolwarm", vmin=-abs(float(n_max)),
                         vmax=abs(float(n_max)))
    for tick in ax[0][0].get_xticklabels():
        tick.set_rotation(45)
    
    # ============================
    # Adjust fig and add colorbar
    # ============================
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.02, 0.7])
    fig.colorbar(im, cax=cbar_ax)
        
    # =======================
    # Label axes and ticks
    # =======================
    ax = label_axes(ax, sizes, p_multi, p_label_res,
                    p_chrom1, p_start1, p_end1,
                    p_chrom2, p_start2, p_end2,
                    bins1, bins2)
    
    # ======================================
    # Set up structure overlay onto heatmap
    # ======================================
    ax = add_structure_overlay(ax, **dict_kwargs)

    # ============
    # Plot Loops
    # ============
    ax = add_loop_overlay(ax, **dict_kwargs)

    if f_out:
        fig.savefig(f_out)

    return fig, ax


def plot_contactmap_w_arcbands(f_pairs, f_cluster, p_chrom1=None, p_chrom2=None,
                               p_start1=None, p_start2=None,
                               p_end1=None, p_end2=None, p_gene1=None, p_gene2=None,
                               p_resolution=100, p_max=None,
                               p_balance=False, p_relative=False, p_self_relative=False,
                               p_normalize=False, p_mode="none",
                               p_vc=False, p_multi=False, p_col=3,
                               p_genome=None, p_structure=False,
                               p_display="full",
                               f_size=None, f_out=None,
                               f_loop=None, f_dot=None, f_fasta=None,
                               f_index=None,
                               p_downsample=1.0, p_figsize=(4, 4), p_dot_size=5.0, p_point=False,
                               p_label_res=1000, p_dpi=100, p_ax=None, p_verbose=True, p_order=None):
    
    if p_genome:
        f_size, f_fasta = get_ref(p_genome)
    
    dict_kwargs = {"f_pairs": f_pairs, "f_cluster": f_cluster, "p_chrom1": p_chrom1, "p_chrom2": p_chrom2,
                   "p_start1": p_start1, "p_start2": p_start2,
                   "p_end1": p_end1, "p_end2": p_end2,
                   "p_gene1": p_gene1, "p_gene2": p_gene2,
                   "p_resolution": p_resolution, "p_max": p_max,
                   "p_balance": p_balance, "p_relative": p_relative, "p_self_relative": p_self_relative,
                   "p_normalize": p_normalize, "p_mode": p_mode,
                   "p_vc": p_vc, "p_multi": p_multi, "p_col": p_col, "p_genome": p_genome,
                   "p_structure": p_structure,
                   "p_display": p_display, "f_size": f_size, "f_out": f_out, "f_loop": f_loop,
                   "f_dot": f_dot, "f_fasta": f_fasta, "f_index": f_index,
                   "p_downsample": p_downsample, "p_figsize": p_figsize,
                   "p_dot_size": p_dot_size, "p_point": p_point, "p_label_res": p_label_res,
                   "p_dpi": p_dpi, "p_ax": p_ax, "p_verbose": p_verbose, "p_order": p_order}

    (loci, p_chrom1, p_start1, p_end1,
     p_chrom2, p_start2, p_end2,
     m_bins, n_bins, bins1, bins2) = get_loci(**dict_kwargs)
    
    sizes = get_size(f_size)    
    raw = get_contactmap(**dict_kwargs)
    p_mode = set_p_mode(**dict_kwargs)
    norm = normalize_contactmap(raw, p_mode=p_mode)
    n_max = set_n_max(norm, p_max, p_mode)
    n_rows, n_cols = set_dims(**dict_kwargs)

    import matplotlib.gridspec as grd
    fig = plt.figure(figsize=(10, 10))
    gs = grd.GridSpec(2, 2,
                      height_ratios=[2, 8],
                      width_ratios=[8, 2],
                      wspace=0.1)

    ax = np.array([[plt.subplot(gs[2])]])  # contactmap
    ax_arcband = plt.subplot(gs[0])  # arcband
    ax_colorbar = plt.subplot(gs[3])  # colorbar
    
    
    if p_verbose:
        print(ax.shape)
        print(loci)        
        for datum in raw:
            print("Dim: %s | Total Counts: %d" % (datum.shape, datum.sum()//2))
        print(n_rows, n_cols)

    # ================
    # Plot Contactmap
    # ================
    for i in range(ax.shape[0]*ax.shape[1]):
        r, c = i//p_col, i%p_col
        indices = (r*p_col)+c
        if p_normalize:
            im = ax[r][c].imshow(norm[indices], cmap="coolwarm",
                                 vmin=-abs(float(n_max)), vmax=abs(float(n_max)))
        else:
            im = ax[r][c].imshow(norm[indices], cmap = "Reds",
                                 vmin=0, vmax=float(n_max), aspect="auto")

        for tick in ax[r][c].get_xticklabels():
            tick.set_rotation(45)

    # ======================================
    # Add colorbar associated w contact map
    # ======================================
    fig.colorbar(im, cax=ax_colorbar)
    
    # =======================
    # Label axes and ticks
    # =======================
    ax = label_axes(ax, sizes, p_multi, p_label_res,
                    p_chrom1, p_start1, p_end1,
                    p_chrom2, p_start2, p_end2,
                    bins1, bins2)

    # =======================
    # Set up ArcBand
    # =======================
    from .utils_arcband import add_arcband_overlay
    ax_arcband = add_arcband_overlay(f_cluster, sizes,
                                     ax_arcband, p_chrom1, p_order)
    
    # ======================================
    # Set up structure overlay onto heatmap
    # ======================================
    ax = add_structure_overlay(ax, **dict_kwargs)

    # ============
    # Plot Loops
    # ============
    ax = add_loop_overlay(ax, **dict_kwargs)

    if f_out:
        fig.savefig(f_out)

    return fig, ax
