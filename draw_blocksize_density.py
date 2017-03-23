from matplotlib import use
use('Agg')
from matplotlib import pyplot as plt
from pyplot_ext import *
import numpy as np
from path import Path
import re
import os

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

def main(p):
    outfile = p.rstrip('/')
    o = open(outfile+'.block_sizes.tsv','w')
    o2 = open(outfile+'.block_sizes.stats.tsv','w')

    blocks = [[],[],[],[]]
    for f in Path(p).files():
        if f.endswith('seq'):
            basename,ext = os.path.splitext(f.split('/')[-1])
            B = parse_block(open(f))
            if len(B) < 4:
                continue
            o.write(basename)
            o2.write(basename)

            for i, B_sub in enumerate(B):
                blocks[i].extend(B_sub)
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

    weighted_blocks = []
    labels = []
    for i, b in enumerate(blocks):
        b = np.array(b)
        weighted_blocks.append(np.ones_like(b)/float(len(b)))
        labels.append('iteration %s' % (i+2))
    return blocks, weighted_blocks, labels

coor = Coor((5, 3),(0.8, 0.2),(.7,.3),(.5,.5),(2,1))
fig1 = plt.figure(**coor.figsize())
print(coor.figsize())
fig1.subplots_adjust(**coor.adjust())

ps = [
    'repo/raw/scop20_training/psiblast_block_BL62_scop20_training_5/0.002/uniref50/',
    'repo/raw/scop20_validation/psiblast_block_BL62_scop20_validation_5/0.002/uniref50/',
    'repo/raw/cath20-scop/psiblast_block_BL62_cath20-scop_5/0.002/uniref50/',
]
titles = ['A','B','C']
for i, p in enumerate(ps):
    ax = fig1.add_subplot(len(ps),1,i+1)
    blocks, weighted_blocks, labels = main(p)
    ax.hist(blocks, bins = 20, range=(0,100), weights=weighted_blocks, histtype='bar', normed=True, label=labels)
    ax.set_xticks(np.arange(0,100,5))
    if not i:
        ax.legend()
    ax.set_xlabel('block size')
    ax.set_ylabel('density')
    ax.set_title(titles[i])
    ax.set_ylim([0,0.08])


print('./fig/blocksize.png')
fig1.savefig('./fig/blocksize.png')

 
