#!/usr/bin/env python3

import matplotlib.pyplot as plt
import sys
import gzip
import argparse

import numpy as np
from sklearn.manifold import TSNE
from sklearn.cluster import SpectralClustering
from sklearn.cluster import DBSCAN
from sklearn.neighbors import kneighbors_graph

import networkx as nx
from networkx.algorithms.community import greedy_modularity_communities

def overlap_function(i1, i2):

    i1_s1, i1_e1, i1_s2, i1_e2 = i1[1], i1[2], i1[3], i1[4]
    i2_s1, i2_e1, i2_s2, i2_e2 = i2[1], i2[2], i2[3], i2[4]

    def calc_overlap(a_s, a_e, b_s, b_e):
        if (b_s > a_e) or (a_s > b_e):
            return 0
        else:
            o_s = max(a_s, b_s)
            o_e = min(a_e, b_e)
            return abs(o_e-o_s)

    left = calc_overlap(i1_s1, i1_e1, i2_s1, i2_e1)
    right = calc_overlap(i1_s2, i1_e2, i2_s2, i2_e2)

    if (left == 0) or (right == 0):
        return 0, 0
    else:
        #return min(left, right)
        return (left, right)

def load_sizes(f_sizes):
    sizes = {}
    with open(f_sizes, "r") as f:
        for line in f:
            row = line.strip("\r\n").split("\t")
            sizes[row[0]] = int(row[1])
    return sizes

def load_eligible_transcript(f_pairs):
    eligible = {}
    with gzip.open(f_pairs, "rt") as f:
        for line in f:
            row = line.strip("\r\n").split("\t")
            if row[1] != row[3]:
                continue
            try:
                eligible[row[1]] += 1
            except KeyError:
                eligible[row[1]] = 1
    return eligible

# currently not functional
def community_detection():
    import community
    partition = community.best_partition(G)

    labels = np.array(list(partition.values()))

    for num in list(set(labels)):

        condition = (labels == num)
        c = np.where(condition)[0]

        if len(c) < 3:
            continue
        out = open("clusters_%s.bed" % num, "w")
        for j in list(c):
            out.write(info[j] + "\n")
        out.close()

def spectral_clustering():
    cluster = SpectralClustering(n_clusters=p_clusters,
                                 affinity="nearest_neighbors",
                                 random_state=386).fit(simimat)
    labels = cluster.labels_

def spectral_clustering_eigengap():
    connectivity = kneighbors_graph(data, n_neighbors=10)
    affinity_matrix = 0.5 * (connectivity + connectivity.T)
    cluster = SpectralClustering(n_clusters=p_clusters,
                                 affinity="precomputed", random_state=386,
                                 assign_labels="kmeans").fit(affinity_matrix)

    k = predict_k(affinity_matrix)
    p_clusters = k
    labels = cluster.labels_

def spectral_clustering_precomputed():
    cluster = SpectralClustering(n_clusters=20, affinity="precomputed",
                                 random_state=386)
    cluster.fit(simimat)
    labels = cluster.labels_

def dbscan():
    cluster = DBSCAN(eps=5, min_samples=3)
    labels = cluster.labels_

def hdbscan():
    cluster = hdbscan.HDBSCAN(min_cluster_size=10, metric="precomputed")
    labels = cluster.fit_predict(distmat)

def clustering():
    print(max(labels))

def call_clusters(data, counter, o_bed, o_output,
                  p_overlap, p_min_reads, p_method="greedy"):
    
    txt = []
    
    G = nx.Graph()
    m = len(data)

    fullmat = np.zeros((m, m))
    simimat = np.zeros((m, m))

    for i in range(len(data)):
        G.add_node(i)

    for i in range(len(data)):
        for j in range(i+1):
            overlap = overlap_function(data[i], data[j])
            simimat[i,j] = min(overlap)
            simimat[j,i] = min(overlap)

            fullmat[i,j] = min(data[i][5], data[i][6], data[j][5], data[j][6])
            fullmat[j,i] = min(data[i][5], data[i][6], data[j][5], data[j][6])


            #if overlap > p_overlap:
            if overlap[0] > p_overlap and overlap[1] > p_overlap:
                G.add_edge(i, j)

    #DEBUG print G.number_of_nodes()
    #DEBUG print G.number_of_edges()

    if p_method == "biconnected":
        clusters = nx.biconnected_components(G)
    elif p_method == "greedy":

        try:
            clusters = greedy_modularity_communities(G)
        except ZeroDivisionError:
            return counter

    labels = np.zeros(m, dtype=int)
    for num, c in enumerate(clusters):

        c = list(c)

        if len(c) < p_min_reads:
            continue

        labels[c] = num

        c1 = [ data[j][0] for j in c ]
        s1 = [ data[j][1] for j in c ]
        e1 = [ data[j][2] for j in c ]
        s2 = [ data[j][3] for j in c ]
        e2 = [ data[j][4] for j in c ]

        s1_m = np.percentile(s1, 50)
        e1_m = np.percentile(e1, 50)
        s2_m = np.percentile(s2, 50)
        e2_m = np.percentile(e2, 50)

        # ==============
        # write median
        # ==============
        #print c1[0], int((s1_m+e1_m)//2), int((s2_m+e2_m)//2)

        # ==============
        # write bed
        # ==============
        chrom = c1[0]
        cid = "%s|%s_%s|%s:%s:%s:%s:%s" % (counter, chrom, num,
                                           chrom, int(s1_m), int(e1_m),
                                           int(s2_m), int(e2_m))
        for j in c:
            start1, start2, length1, length2, strand = (data[j][1], data[j][3],
                                                        data[j][5], data[j][6],
                                                        data[j][8])

            bed_coords = [chrom, start1, start2+length2, cid, 1000,
                          strand, start1, start2+length2, 0, 2,
                          "%s,%s" % (length1, length2),
                          "%s,%s" % (0, start2-start1)]
            
            if o_bed is not None:
                o_bed.write("\t".join(map(str, bed_coords)) + "\n")

        # =============
        # write txt
        # =============
        out_coords = [c1[0], int(s1_m), int(e1_m), c1[0], int(s2_m), int(e2_m)]
        txt.append(out_coords)
        if o_output is not None:
            o_output.write("\t".join(map(str, out_coords)) + "\n")

        counter += 1

    return counter, txt

def visualize(simimat, fullmat, f_tsne=None, ax=None):

    distmat = (fullmat - simimat)

    # compute embedding
    model = TSNE(random_state=386, metric="precomputed")
    embeddings = model.fit_transform(distmat)

    # visualization
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(embeddings[:,0], embeddings[:,1], s=4, c=labels, cmap="Spectral")

    for i in range(max(labels)+1):
        condition = (labels == i)

        ax.annotate(i, embeddings[condition,:].mean(axis=0),
                    horizontalalignment="center",
                    verticalalignment="center", size=10,
                    weight="bold",
                    color="black")

    if f_tsne is not None:
        plt.savefig(f_tsne, dpi=300)

    return fig, ax

def cluster_chimeric_pairs(p_iid, p_genome,
                           p_long, p_span, p_overlap, p_min_reads,
                           p_min_interactions, p_max_interactions,
                           f_sizes, f_tsne, f_pairs, f_output, f_bed):
    
    results = []
    
    eligible = load_eligible_transcript(f_pairs)

    if f_output is not None:
        o_bed = open(f_bed, "w")
        o_output = open(f_output, "w")
    else:
        o_bed = None
        o_output = None
    
    current_iid = None
    counter = 0
    data = []
    with gzip.open(f_pairs, "rt") as f:

        for line in f:

            row = line.strip("\r\n").split("\t")

            if row[1] != row[3]:
                continue

            if row[1] != "all" and row[1] != p_iid:
                continue
            
            if current_iid == None:
                current_iid = row[1]

            elif current_iid != row[1]:
                if (eligible[current_iid] >= p_min_interactions
                    and eligible[current_iid] <= p_max_interactions):
                    counter, txt = call_clusters(data, counter, o_bed, o_output,
                                                 p_overlap, p_min_reads)
                    results.append(txt)

                data = []
                current_iid = row[1]


            # append these information

            c1, s1, l1, c2, s2, l2 = (row[1], int(row[2]), int(row[9]),
                                      row[3], int(row[4]), int(row[10]))
            iid, strand1, strand2 = (row[0], row[5], row[6])

            if abs(s1-(s2+l2)) < p_long:
                continue
            
            if abs(s2-(s1+l1)) < p_span:
                continue

            if strand1 == strand2:
                if strand1 == "-":
                    strand = "+"
                else:
                    strand = "-"

            data.append((c1, s1, s1+l1, s2, s2+l2, l1, l2, iid, strand))

    # ===========
    # last call
    # ===========
    if (eligible[current_iid] >= p_min_interactions
        and eligible[current_iid] <= p_max_interactions):
        counter, txt = call_clusters(data, counter, o_bed, o_output,
                                     p_overlap, p_min_reads)
        results.append(txt)
        
    if f_output is not None:
        o_bed.close()
        o_output.close()

    return results

def main(p_iid, p_genome,
         p_long, p_span, p_overlap, p_min_reads,
         p_min_interactions, p_max_interactions,
         f_sizes, f_tsne, f_pairs, f_output, f_bed):

    cluster_chimeric_pairs(p_iid, p_genome,
                           p_long, p_span, p_overlap, p_min_reads,
                           p_min_interactions, p_max_interactions,
                           f_sizes, f_tsne, f_pairs, f_output, f_bed)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-c", default="all", help="all|id")
    parser.add_argument("-b", default=None, help="bed")
    parser.add_argument("-o", default=None, help="clusters")
    parser.add_argument("-t", default=None, help="tsne")
    parser.add_argument("-s", default=None, help="sizes file")
    parser.add_argument("-g", default=None, help="genome")
    parser.add_argument("-l", default=200, type=int, help="long range threshold")
    parser.add_argument("-p", default=10, type=int, help="span threshold")
    parser.add_argument("-ol", default=10, type=int, help="overlapping threshold")
    parser.add_argument("-m", default=3, type=int, help="minimum no. of reads per cluster")
    parser.add_argument("-i", default=200, type=int, help="minimum no. of interactions per transcript")
    parser.add_argument("-j", default=20000, type=int, help="max no. of interactions per transcript")
    parser.add_argument("-pairs", default=None, help="paired txt file")
    args = parser.parse_args()

    p_iid = args.c
    p_genome = args.g
    p_long = args.l
    p_span = args.p
    p_overlap = args.ol
    p_min_reads = args.m
    p_min_interactions = args.i
    p_max_interactions = args.j

    f_sizes = args.s
    f_tsne = args.t
    f_pairs = args.pairs
    f_output = args.o
    f_bed = args.b

    main(p_iid=p_iid, p_genome=p_genome,
         p_long=p_long, p_span=p_span, p_overlap=p_overlap,
         p_min_reads=p_min_reads,
         p_min_interactions=p_min_interactions,
         p_max_interactions=p_max_interactions,
         f_sizes=f_sizes, f_tsne=f_tsne, f_pairs=f_pairs, f_output=f_output, f_bed=f_bed)
    
