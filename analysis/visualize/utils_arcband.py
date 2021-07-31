#!/usr/bin/env python3

import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches

import numpy as np
import argparse

proj_dir = ""

def calc_height(left, right):
    x1, y1 = left, 0.0
    x2, y2 = right, 0.0
    mxmy = mx, my = [(x1 + x2) / 2, (y1 + y2) / 2]
    r = np.sqrt((x1 - mx)**2 + (y1 - my)**2)
    width = 2 * r
    height = 2 * r
    start_angle = np.arctan2(y1 - my, x1 - mx) * 180 / np.pi
    end_angle = np.arctan2(my - y2, mx - x2) * 180 / np.pi
    return height

def get_verts(coords):
    s1, e1, s2, e2 = coords
    m1 = (s1 + e2)/2.0
    m2 = (s2 + e1)/2.0

    h1 = calc_height(s1, e2)
    h2 = calc_height(s2, e1)

    verts = [(s1, 0.0),
             (m1, h1),
             (e2, 0.0),
             (s2, 0.0),
             (m2, h2),
             (e1, 0.0),
             (0.0, 0.0)]
    return verts

def get_sizes(f_size):
    sizes = {}
    with open(f_size, "r") as f:
        for line in f:
            row = line.strip("\r\n").split("\t")
            chrom, size = row[0], int(row[1])
            sizes[chrom] = size
    return sizes


def add_arcband_overlay(f_cluster, sizes, ax, p_chrom, p_order):
    p_genomic_height = sizes[p_chrom]//2
    p_genomic_width = sizes[p_chrom]
    
    coords = []
    with open(f_cluster, "r") as f:
        for line in f:
            row = line.strip("\r\n").split("\t")
            if row[0] == p_chrom:
                s1, e1, s2, e2, lt = (int(row[1]), int(row[2]),
                                      int(row[4]), int(row[5]), row[6])

                if e1 <= s2:
                    coords.append((s1, e1, s2, e2, lt))
                else:
                    coords.append((s1, e1, e1, e2, lt))

    type_to_color = {"loops": "royalblue",
                     "lstripe": "darkorange",
                     "rstripe": "deeppink"}
    
    if p_order:
        if p_order == "by_width":
            coords = sorted(coords, key = lambda q: abs(q[0]-q[3]), reverse=True)
        elif p_order == "by_category":
            coords = sorted(coords, key = lambda q: q[4])
    
    for s1, e1, s2, e2, lt in coords:
        
        coord = (s1, e1, s2, e2)
        verts = get_verts(coord)

        codes = [
            Path.MOVETO,
            Path.CURVE3,
            Path.CURVE3,
            Path.LINETO,
            Path.CURVE3,
            Path.CURVE3,
            Path.CLOSEPOLY,
        ]

        path = Path(verts, codes)

        patch = patches.PathPatch(path, facecolor=type_to_color[lt],
                                  lw=0.5, alpha=0.5)
        ax.add_patch(patch)

    ax.set_xlim(0, p_genomic_width)
    ax.set_ylim(0, p_genomic_height)
    return ax


def plot_arcband(p_chrom, p_genome, p_order, f_size, f_cluster, f_output):

    if p_genome:
        if p_genome == "hg19":
            f_size = proj_dir + "data/reference/hg19/hg19_refseq_wo_version.sizes"
        elif p_genome == "mm10":
            f_size = proj_dir + "data/reference/mm10/mm10_refseq.sizes"
    
    sizes = get_sizes(f_size)

    fig, ax = plt.subplots(figsize=(8, 4))
    ax = add_arcband_overlay(ax)
    
    if f_output is not None:
        plt.savefig(f_output, dpi=300)

    return fig, ax

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-c", default=None, help="chromosome")
    parser.add_argument("-s", default=None, help="sizes file")
    parser.add_argument("-g", default=None, help="genome")
    parser.add_argument("-r", default=None, help="order arcbands by")
    parser.add_argument("-i", default=None, help="input cluster file")
    parser.add_argument("-o", default=None, help="output")
    
    args = parser.parse_args()

    p_chrom = args.c
    p_genome = args.g
    p_order = args.r
    f_size = args.s
    f_cluster = args.i
    f_output = args.o

    plot_arcband(p_chrom, p_genome, p_order, f_size, f_cluster, f_output)

