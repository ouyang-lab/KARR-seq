#!/usr/bin/env python3

import sys

def main(f_app):
    current_pair = None
    current_row = None
    for no, line in enumerate(sys.stdin):
        if f_app == "dedup":
            row = line.strip("\r\n").split("\t")
            pair = (row[1], int(row[2]), row[3], int(row[4]))
            if no == 0:
                current_pair = pair
                current_row = row
            else:
                if pair == current_pair:
                    # skip
                    pass
                else:
                    print("\t".join(current_row))
                    current_pair = pair
                    current_row = row
        else:
            print(line.strip("\r\n"))
            
    if f_app == "dedup":
        print("\t".join(current_row))

if __name__ == "__main__":
    f_app = sys.argv[1]
    main(f_app)
    
    
