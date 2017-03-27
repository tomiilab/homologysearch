import numpy as np
from path import Path
from collections import defaultdict
import re
import os
import sys

def parse_block(ifh):
    B = []
    B_sub = []
    prev_s = None
    prev_e = None
    for l in ifh:
        if l.startswith("###"):
            pos, s, e = filter(lambda a:a, re.split('\D+', l[3:].strip()))
            pos = int(pos)
            e = int(e)
            s = int(s)
            if pos == 0 and B_sub:
                B.append(B_sub)
                B_sub = []
            if prev_s != s and prev_e != e:
                prev_e = e
                prev_s = s
            B_sub.append(e-s+1)
    B.append(B_sub)
    return B

def main(p, max_iteration=5):
    p = Path(p)

    if not p.isdir():
        sys.stderr.write("%s is not a dir\n" % (p))
        sys.exit(0)

    result_files = []
    for f in Path(p).files():
        if f.endswith('.seq'):
            result_files.append(f)
 
    if not len(result_files):
        sys.stderr.write("Unable to find files in %s\n" % (p))
        sys.exit(0)

    print("parsing %s " % p)
    print("num of files: %s " % len(result_files))
    outfile = p.rstrip('/')

    oname = outfile+'.block_sizes.tsv'
    oname2 = outfile+'.block_sizes.stats.tsv'

    o = open(oname,'w')
    o2 = open(oname2,'w')

    for f in result_files:
        basename,ext = os.path.splitext(f.split('/')[-1])
        B = parse_block(open(f))

        o.write(basename)
        o2.write(basename)

        for i, B_sub in enumerate(B):
            o.write('\t')
            o.write(','.join(map(str, B_sub)))

            B_sub = np.array(B_sub)
            o2.write('\t')
            o2.write(str(np.count_nonzero(B_sub==1)/float(len(B_sub))))
            o2.write(',')
            o2.write(str(np.count_nonzero(B_sub<5)/float(len(B_sub))))
        o.write('\n')
        o2.write('\n')

    o.close()
    o2.close()
    print("saved at %s" % (oname))
    print("saved at %s" % (oname2))

if __name__ == '__main__':
    import argparse
    par = argparse.ArgumentParser()
    par.add_argument('path')
    args = par.parse_args()
    main(args.path)
    
