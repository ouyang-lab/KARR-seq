#!/usr/bin/env python3

import sys
import re
from cigar import Cigar

def get_strand(flag):
    if (flag & 16) == 16:
        return "-"
    else:
        return "+"

def get_cigar_match(cigar):
    soft = 0
    matches_M = []
    matches_N = []
    for m, n in list(Cigar(cigar).items()):
        if n == "M":
            matches_M.append(m)
        elif n == "N":
            matches_N.append(m)
    
    if (len(matches_M) >= 2) and (len(matches_N)>=1):
        return matches_M[0], matches_M[1], matches_N[0]
    else:
        return None, None, None

def get_cigar_chimeric(cigar):
    matches = []
    for m, n in list(Cigar(cigar).items()):
        if n == "M":
            matches.append(m)
    if len(matches) >= 1:
        return matches[0]
    else:
        return None

    
def main(f_app, mapq_threshold, span_threshold):
    
    prev_id = None
    prev_list = []
    
    for line in sys.stdin:
        if line.startswith("@"):
            continue
        
        row = line.strip("\r\n").split("\t")
        
        iid = row[0]
        flag = int(row[1])
        chrom = row[2]
        start = int(row[3])  # SAM: 1-based
        qual = int(row[4])
        cigar = row[5]
        
        # get flag id
        strand = get_strand(flag)
        
        if f_app == "Aligned":
            
            if "N" in cigar:
                
                # get fragment length
                # get gap
                l1, l2, gap = get_cigar_match(cigar)
                if gap:

                    s1 = start
                    s2 = start+l1+gap

                    if qual >= mapq_threshold:

                        if abs(s1-s2) >= span_threshold:
                            
                            print("\t".join(map(str, [iid,
                                                      chrom, start,
                                                      chrom, start+l1+gap,
                                                      strand, strand,
                                                      qual, qual,
                                                      l1, l2, 
                                                      "L", "R", 0])))
            
        if f_app == "Chimeric":
            
            # get fragment length
            l = get_cigar_chimeric(cigar)
            
            if prev_id == None:
                prev_list.append((chrom, start, qual, iid, l, strand))
                prev_id = iid
            elif (prev_id == row[0]):
                prev_list.append((chrom, start, qual, iid, l, strand))
            else:
                # unload previous list
                if len(prev_list) == 2:

                    mapq1 = prev_list[0][2]
                    mapq2 = prev_list[1][2]
                    s1 = prev_list[0][1]
                    s2 = prev_list[1][1]
                    
                    if mapq1 >= mapq_threshold and mapq2 >= mapq_threshold:

                        if abs(s1-s2) >= span_threshold:
                            
                            if prev_list[0][1] < prev_list[1][1]:
                                print("\t".join(map(str, [prev_list[0][3],
                                                          prev_list[0][0], prev_list[0][1],
                                                          prev_list[1][0], prev_list[1][1],
                                                          prev_list[0][5], prev_list[1][5],
                                                          prev_list[0][2], prev_list[1][2],
                                                          prev_list[0][4], prev_list[1][4],
                                                          "R", "L", 0])))
                            else:
                                print("\t".join(map(str, [prev_list[0][3],
                                                          prev_list[1][0], prev_list[1][1],
                                                          prev_list[0][0], prev_list[0][1],
                                                          prev_list[1][5], prev_list[0][5],
                                                          prev_list[1][2], prev_list[0][2],
                                                          prev_list[1][4], prev_list[0][4],
                                                          "R", "L", 1])))

                # reset previous list
                prev_list = [(chrom, start, qual, iid, l, strand)]
                prev_id = iid


if __name__ == "__main__":

    f_app = sys.argv[1]  # Aligned || Chimeric
    mapq_threshold = int(sys.argv[2])
    span_threshold = int(sys.argv[3])
    
    main(f_app, mapq_threshold, span_threshold)
