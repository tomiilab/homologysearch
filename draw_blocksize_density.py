import numpy as np
from path import Path
from collections import defaultdict
import os


def main(p, max_iteration=5):
    print("parsing %s " % p)

    blocks = []
    for _ in range(max_iteration-1):
        blocks.append([])
 
    iter2num = defaultdict(int)
    with open(p) as f:
        for l in f:
            cs = l.strip().split()
            name = cs[0]
            B = cs[1:]
   
            for i in range(len(B)):
                iter2num[i+2] += 1

            if len(B) >= max_iteration-1:
                for i in range(max_iteration-1):
                    blocks[i].extend([int(e) for e in B[i].split(',')])


    iters = iter2num.keys()
    iters.sort()
    for i in iters:
        print("iteration %s: %s hits" % (i, iter2num[i]))

    weighted_blocks = []
    labels = []
    for i, b in enumerate(blocks):
        b = np.array(b)
        weighted_blocks.append(np.ones_like(b)/float(len(b)))
        labels.append('iteration %s' % (i+2))
    return blocks, weighted_blocks, labels


if __name__ == '__main__':
    import argparse
    from matplotlib import use
    use('Agg')
    from matplotlib import pyplot as plt
    from pyplot_ext import *


    par = argparse.ArgumentParser()
    par.add_argument('max_iteration', type=int, default=5)
    par.add_argument('-max_block_size', type=int, default=100)
    par.add_argument('-num_bins', type=int, default=20)
    par.add_argument('-ylim', type=float, default=0.40)
    args = par.parse_args()
    print("show until iteration %s" % (args.max_iteration))
    coor = Coor((5, 3),(0.8, 0.2),(.7,.3),(.5,.5),(2,1))
    fig1 = plt.figure(**coor.figsize())
    fig1.subplots_adjust(**coor.adjust())

    ps = [
        'repo/raw/scop20_training/psiblast_block_BL62_scop20_training_8/0.002/uniref50.block_sizes.tsv',
        'repo/raw/scop20_validation/psiblast_block_BL62_scop20_validation_8/0.002/uniref50.block_sizes.tsv',
        'repo/raw/cath20-scop/psiblast_block_BL62_cath20-scop_8/0.002/uniref50.block_sizes.tsv',
    ]


    titles = ['A','B','C']
    for i, p in enumerate(ps):
        ax = fig1.add_subplot(len(ps),1,i+1)
        blocks, weighted_blocks, labels = main(p, max_iteration=args.max_iteration)
        ax.hist(blocks, bins = args.num_bins, range=(0, args.max_block_size), weights=weighted_blocks, histtype='bar', label=labels, edgecolor='none')
        ax.set_xticks(np.arange(0,args.max_block_size, args.max_block_size/float(args.num_bins)))
        if not i:
            ax.legend(fontsize=9)
        ax.set_xlabel('block size')
        ax.set_ylabel('density')
        ax.set_title(titles[i])
        ax.set_ylim([0,args.ylim])

    outfig = './fig/blocksize.%s.%s.%s.png' % (args.max_iteration, args.max_block_size, args.num_bins)
    print("saved at %s" % (outfig))
    fig1.savefig(outfig, dpi=300)

     
