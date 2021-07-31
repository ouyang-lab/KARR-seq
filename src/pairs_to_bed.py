#!/usr/bin/env python3

import sys

def main(p_mapq_threshold,
         p_min_dist_threshold, p_max_dist_threshold):
    
    for no, line in enumerate(sys.stdin):
        
        row = line.strip("\r\n").split("\t")
        
        if (row[1] == row[3] and
            int(row[7]) >= p_mapq_threshold and
            int(row[8]) >= p_mapq_threshold):  # only display intra and specified MAPQ
            
            iid = row[0]
            chrom1, s1, strand1, o1 = row[1], int(row[2]), row[5], int(row[9])
            chrom2, s2, strand2, o2 = row[3], int(row[4]), row[6], int(row[10])

            s1 = s1 - 1  # bed is 0-based, sam is 1-based
            
            span = abs(s1-(s2+o2))
            
            if strand1 == strand2:
                if strand1 == "-":
                    strand = "+"
                else:
                    strand = "-"
            else:
                strand = "."

            if span >= p_min_dist_threshold and span <= p_max_dist_threshold:
                print("\t".join(map(str, [chrom1, s1, s2+o2, iid, 1000,
                                          strand, s1, s2+o2, 0, 2,
                                          "%s,%s" % (o1, o2),
                                          "%s,%s" % (0, s2-s1)])))
            
if __name__ == "__main__":

    p_mapq_threshold = int(sys.argv[1])
    p_min_dist_threshold = 83
    p_max_dist_threshold = 100000
    
    main(p_mapq_threshold,
         p_min_dist_threshold,
         p_max_dist_threshold)
