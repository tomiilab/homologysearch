import ROC
from collections import defaultdict
from itertools import chain
from scipy import stats
from matplotlib import use
use('Agg')
from matplotlib import pyplot as plt
from matplotlib import rcParams, rc
import numpy as np
from pyplot_ext import *
from path import Path
fontsize = 8
rcParams.update({'font.size': fontsize})


class Benchmarks(list):
    def __init__(self, benchmark_name, args, fdr, standard, ax=None, max_evalue=10.0, max_FP_count=5, rocn_max_evalue=None):
        self.fdr = fdr
        self.standard = standard
        self.max_evalue = max_evalue
        self.max_FP_count = max_FP_count
        self.args = args

        self.MEAN_ROCN = []
        self.labels = []

        self.proteins = ROC.Proteins('./repo/fastadb/'+benchmark_name, standard, JG_path='./data/jg')

        for arg in self.args:

            hits = self.eval(benchmark_name, *arg[:5], db=arg[5], max_evalue=max_evalue)

            style = arg[-1]

            bench = ROC.Benchmark(hits)
            bench.iter = arg[3]
            bench.cal_ROC(fdr=fdr)
            if ax is not None:
                bench.draw(ax, fdr=fdr, style=style)
            
            if rocn_max_evalue is None:
                print("evalue cut for ROCN: %s" % bench.evalue_fdr)
                bench.cal_ROCN(max_evalue=bench.evalue_fdr, max_FP_count=max_FP_count)
            else:
                print("evalue cut for ROCN: %s" % rocn_max_evalue)
                bench.cal_ROCN(max_evalue=rocn_max_evalue, max_FP_count=max_FP_count)

            rocns = [e.ROCN for e in bench.ROCN_benchmarks]
            rocns = rocns + [0.0]*(self.proteins.valid_seqnum-len(rocns))
            print('ROC{}: {}'.format(max_FP_count, np.mean(rocns)))

            self.MEAN_ROCN.append(np.mean(rocns))
            self.append(bench)
            print('#End\n')

    def eval(self, bench, method, matrix, subject, n_iter, e_iter, db=None, max_evalue=10.0):
        #dir containing result
        dirpath = './repo/raw/'+bench+'/'+method+'_'+matrix+'_'+subject
        if 'msa' in method: 
            dirpath += '_5'
        elif n_iter > 1:
            dirpath += '_'+str(n_iter)

        if e_iter != None:
            dirpath += '/'+str(e_iter)
        if db != None:
            dirpath += '/'+db

        print(dirpath)       
        print('num of iterations: %s' % n_iter) 

        i = 0
        group = []
        
        P = Path(dirpath)
        for f in P.walkfiles(pattern="*.tsv"):
            hits = self.proteins.parse_hits(f, max_evalue=max_evalue)
            if hits != None:
                group.append(hits)
            i += 1

        print(method+' file num: '+str(i))
        print(method+' query num: '+str(len(group)))
        return group


    def compare_ROCN(self):

        name2iter2block1 = defaultdict(dict)
        name2iter2block4 = defaultdict(dict)

        fs = [
            './repo/raw/scop20_training/psiblast_block_BL62_scop20_training_8/0.002/uniref50.block_sizes.stats.tsv',
            './repo/raw/scop20_validation/psiblast_block_BL62_scop20_validation_8/0.002/uniref50.block_sizes.stats.tsv',
            './repo/raw/cath20-scop/psiblast_block_BL62_cath20-scop_8/0.002/uniref50.block_sizes.stats.tsv'
        ]

        for f in fs:
            for l in open(f):
                cs = l.strip().split()
                for i, e in enumerate(cs[1:]):
                    b1,b4 =e.split(',')
                    name2iter2block1[cs[0]][i+2] = float(b1)
                    name2iter2block4[cs[0]][i+2] = float(b4)


        name2seqlen = {}
        name2n_TP = {}
        name2rocn1 = {}
        name2rocn2 = {}

        bench1_iter = int(self[0].iter)
        print(bench1_iter) 

        names = []
        for bench in self[0].ROCN_benchmarks:
            name2seqlen[bench[0].id] = bench[0].seqlen
            name2rocn1[bench[0].id] = bench.ROCN
            name2n_TP[bench[0].id] = bench[0].size

        for bench in self[1].ROCN_benchmarks:
            name2seqlen[bench[0].id] = bench[0].seqlen
            name2rocn2[bench[0].id] = bench.ROCN
            name2n_TP[bench[0].id] = bench[0].size

        print(len(name2seqlen))
        print(len(name2iter2block1))

        X = []
        Y = []
        B1 = []
        B4 = []
    
        for name, size in self.proteins.id2size.items():
            if size <= 1:
                continue
            seqlen = name2seqlen.get(name)
            rocn1 = name2rocn1.get(name)
            rocn2 = name2rocn2.get(name)
            if rocn1 is None:
                rocn1 = 0
            if rocn2 is None:
                rocn2 = 0
            d = rocn2-rocn1

            b1 = name2iter2block1[name].get(bench1_iter)
            b4 = name2iter2block4[name].get(bench1_iter)
            if b1 is not None and b4 is not None:
                X.append(seqlen)
                Y.append(d)
                B1.append(b1)
                B4.append(b4)
                print('#%s\t%s\t%s\t%s\t%s\t%s\t%s' % (name, seqlen, name2n_TP[name], rocn1, rocn2, b1, b4))

        print('num of hits: %s' % (len(X)))
        print('sum of d: %s' % np.sum(Y))
        print('pearsonr len vs d: %s' % stats.pearsonr(X,Y)[0])
        print('pearsonr block1 vs d: %s' % stats.pearsonr(B1,Y)[0])
        print('pearsonr block4 vs d: %s' % stats.pearsonr(B4,Y)[0])


def draw_SF(args, figfile, standards, xlim=500, ylim=900, dpi=300, max_evalue=10.0, fdr = 0.1, rocn_max_evalue=None):
    
    n_panel = len(standards)

    coor = Coor((2.7, 2.7),(0.6, 0.1),(.5,.2),(.3,.3),(n_panel,1))
    fig1 = plt.figure(**coor.figsize())
    fig1.subplots_adjust(**coor.adjust())
    print(coor.figsize())

    B = []
    labels = [d[-1]['label'] for d in args['args']]

    AX = []
    for i, standard in enumerate(standards):

        AX.append(fig1.add_subplot(n_panel, 1, i+1))
        
        AX[i].set_xlabel('weighted false positive count')
        AX[i].set_ylabel('weighted true positive count')
        AX[i].set_ylim([0,ylim])
        AX[i].set_xlim([0,xlim])
        AX[i].grid()
        if n_panel > 1:
            AX[i].set_title('(%s)' % chr(65+i))

        benchmarks = Benchmarks(args['benchmark'], args['args'], fdr, standard, ax=AX[i], max_evalue=max_evalue, rocn_max_evalue=rocn_max_evalue)
        B.append(benchmarks)
        if i == 0:
            if len(labels) < 9:
                AX[i].legend(loc='lower right',fontsize=7)
            else:
                AX[i].legend(loc='lower right', ncol=2, fontsize=7)


    figfile_roc = figfile+'.png'
    
    print(figfile_roc)
    fig1.savefig(figfile_roc, dpi=dpi)

    return benchmarks

def calc_rocn_stats(args, standard, max_evalue=10.0, fdr = 0.1, rocn_max_evalue=None):
    benchmarks = Benchmarks(args['benchmark'], args['args'], fdr, standard, max_evalue=max_evalue, rocn_max_evalue=rocn_max_evalue)
    benchmarks.compare_ROCN()

def draw_mean_rocn(ax, M, labels, legends, cs):
    width = 0.35
    rects = []
    legs = []
    for i, m in enumerate(M):
        r = ax.bar(np.arange(len(m))+(width*(i+.5)), m, width, color=cs[i])
        rects.append(r)
    ax.legend(rects,legends, fontsize=8)
    ax.set_ylabel('mean ROC5')
    ax.set_xlim([0,len(m)])
    ax.set_xticks(np.arange(len(m))+width)
    ax.set_xticklabels(labels, rotation=70)


if __name__ == '__main__':
    import argparse, json
    par = argparse.ArgumentParser()
    par.add_argument('action', choices=['draw_roc', 'roc5_stats'])
    par.add_argument('args')
    par.add_argument('-xlim', type=int, default=500)
    par.add_argument('-ylim', type=int, default=900)
    par.add_argument('-max_evalue', type=float, default=1)
    par.add_argument('-rocn_max_evalue', type=float)
    par.add_argument('-standards', nargs='+', choices=["fold", "superfamily", "JG"])
    par.add_argument('--compare', action="store_true")
    args = par.parse_args()

    ARGS = json.load(open(args.args))
    if args.action == "draw_roc":
        name = ARGS['benchmark']
        draw_SF(ARGS, './fig/ROC_'+name+'_'+'_'.join(args.standards), args.standards, xlim=args.xlim, ylim=args.ylim, dpi=600, max_evalue=args.max_evalue, rocn_max_evalue=args.rocn_max_evalue)

    elif args.action == "roc5_stats":
        calc_rocn_stats(ARGS, args.standards[0], max_evalue=args.max_evalue, rocn_max_evalue=args.rocn_max_evalue)

