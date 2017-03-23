import argparse, pathlim, sys
import numpy as np
from seqlim import MSeq
import parse_aln as pa
from pyplot_ext import *

from matplotlib import use
use('Agg')
from matplotlib import pyplot as plt
from matplotlib import rcParams, rc
fontsize = 8
rcParams.update({'font.size': fontsize})


class Align_match:
    '''
    calculate accuracy and sensitivity
    '''
    def __init__(self, tags, bench, ref_aln_dir, ref_aln_ext='.dal'):
        '''
        tags : tags of seqs
        bench : path to the benchmark seqs
        '''
        self.error=False

        self.tags = tags
        
        ref_aln_seqs = self.parse_ref_aln(open(ref_aln_dir+'/'+'.'.join(self.tags)+ref_aln_ext).read())

        print '\noriginal ref aln'
        print '\n'.join(ref_aln_seqs)

        ref_aln = np.array([[e] for e in ref_aln_seqs]).view('S1')
    
        true_seq0 = MSeq.parse(open(bench+'/'+self.tags[0]+'.seq'))[0].seq.upper()
        true_seq1 = MSeq.parse(open(bench+'/'+self.tags[1]+'.seq'))[0].seq.upper()
        
        
        r = self.find_missing(ref_aln, [true_seq0, true_seq1])
        if not r[0]:
            self.error=True
        else:
            self.seqarr_true = r[1]
            print '\nmarked ref aln'
            for e in self.seqarr_true.view('S'+str(self.seqarr_true.shape[1])):
                print e[0]    


    def find_missing(self, aln, true_seqs):
        j = 0
        i = 0
        while 1:
            try:
                r = aln[0][j]
            except:
                break
            
            try:
                r0 = true_seqs[0][i]
            except:
                return None, None
            
            if r in ['-','X']:
                j += 1
            elif r == r0:
                j += 1
                i += 1
            else:
                aln = np.insert(aln, j , ['#','-'], axis=1)
                i += 1
                j += 1

        while 1:
            try:
                r = aln[1][j]
            except:
                break
            
            try:
                r0 = true_seqs[1][i]
            except:
                return None, None
            
            if r in ['-','X']:
                j += 1
            elif r == r0:
                j += 1
                i += 1
            else:
                aln = np.insert(aln, j , ['-','#'], axis=1)
                i += 1
                j += 1
        return True, aln


    def parse_ref_aln(self, seqtxt):
        '''
        parse sequences from TM-align format
        '''
        lines = seqtxt.splitlines()
        if 'TM-align' in seqtxt: 
            return [lines[-4].upper(),lines[-2].upper()]
        else:
            return [e.split()[1].upper().replace('X','-') for e in lines[1:]]

    def load_aln(self, aln, query_start=1, target_start=1):
        '''
        load alignments to be test with their starting positions
        '''
        self.query_start = query_start
        self.target_start = target_start
        if self.query_start== 1 and self.target_start==1:
            self.starts = None
        else:
            self.starts = np.array([self.query_start, self.target_start])
        self.seqarr = np.array([[e] for e in aln]).view('S1')

        print '\nalign'
        print '\n'.join(aln)

    def seqarray2indice(self, res_ref_str, starts=None):
        '''
        convert seqarray to indice_array which contains gap-ignored residue postions.
        '''
        res_ref = np.ones(res_ref_str.shape, dtype='int')
        #- and . are regarded as gap.
        gap_arg = np.logical_or(res_ref_str=='-',res_ref_str=='.')
        res_ref[gap_arg] = 0
        
        if np.any(starts) == None:
            res_indice = np.cumsum(res_ref, axis=-1)
        else:
            res_indice = np.cumsum(res_ref, axis=-1) + (starts - 1)[:,None]
        #gap positions are labelled 0
        res_indice[gap_arg] = 0
        return res_indice

    def match_algs(self):
        '''
        label correctly and incorrectly aligned columns 1 and 0, respectively.
        label ignored columns -1
        '''
        #convert seq array to indice array that marks individual positions ignoring gaps
        ind_test = self.seqarray2indice(self.seqarr, starts=self.starts)
        ind_true = self.seqarray2indice(self.seqarr_true)
        ind_true = ind_true[:,(self.seqarr_true=='#').sum(axis=0)==0]
        ind_true = ind_true[:,(ind_true==0).sum(axis=0)==0]
        #ind_test = ind_test[:,(ind_test==0).sum(axis=0)==0]
        gap_p = ind_test == 0    

        #num of seqs    
        seq_num = ind_true.shape[0]
        col_match = np.zeros((ind_test.shape[1]))
        for i in range(ind_true.shape[1]):
            arg = (ind_test - ind_true[:,i][:,None]) == 0
            arg[gap_p] = False
            col_match[np.sum(arg, axis=0) == 2] = 1
        col_match[gap_p.sum(axis=0)!=0] = -1
        
        print '\nmatched res'
        temparr = self.seqarr.copy()
        temparr[:,col_match!=1] = 'x'
        for e in temparr:
            print ''.join(e.tolist())
        
        return col_match    
    
    def ref_identity(self):
        I = np.logical_or(self.seqarr_true == '-', self.seqarr_true == '#').sum(axis=0) == 0

        L= I.sum()
        arg = self.seqarr_true[0][I] == self.seqarr_true[1][I]
        
        ma = np.count_nonzero(arg)
        return ma/float(L), ma, L

def draw_roc(ax, D1, D2, label, c, ls):
    ax.plot(np.cumsum(1.0-D2),np.cumsum(D2), label=label.upper(),color=c,linewidth=ls)
    return ax

def draw(args, figfile):

    coor = Coor((2.8, 2.8),(.5,.1),(.5,.5),(.3,.3),(1,1))
    print coor.figsize()
    fig = plt.figure(**coor.figsize())
    fig.subplots_adjust(**coor.adjust())


    ax1 = fig.add_subplot(111)
    ax1.set_xlabel('sum(1 - precision)')
    ax1.set_ylabel('sum(sensitivity)')
    ax1.set_xlim([0,200])
    ax1.set_ylim([0,100])

    for arg in args:
        label,c,ls = arg[-3:]
        draw_roc(ax1, arg[0], arg[1], label,c,ls)

    ax1.legend(loc='lower right',fontsize=8)
    ax1.grid()
    fig.savefig(figfile+'.png', dpi=600)

def cal(fs, bench, ref_aln_dir, ref_aln_ext):
    M = []
    N = []
    for i,f in enumerate(fs):
    
        print '\n'+f
        r = pa.parse(open(f).read())
        #print r
        if not r:
            print f
            continue
        query, target =f.split('/')[-1].split('.')[:2]
        am = Align_match([query, target], bench, ref_aln_dir, ref_aln_ext=ref_aln_ext)
        if am.error:
            continue
        iden, ma, tot = am.ref_identity() # tot is total length of residue-residue alignment in ref alignment
        matchs = []
        for each in r[0]:
            q, t, e = each[:3]
            qid, qseq, qs, qe = q
            tid, tseq, ts, te = t
            am.load_aln([qseq,tseq], query_start=qs, target_start=ts)
            col_match = am.match_algs()
            subm = np.array([np.count_nonzero(col_match==1),np.count_nonzero(col_match==0)])
            print('evalue:%s' % e)
            print('match:%s non-match:%s' % tuple(subm))
            if subm.sum() == 0:
                continue
            else:
                prec = subm[0]/float(subm.sum())
                sens = subm[0]/float(tot)
                print('Prec:%s Sens:%s' % (prec, sens))
            matchs.append(subm)
            N.append([e,prec,sens])
            break

        #print matchs
        matchs = np.array(matchs,dtype=float).sum(axis=0)
        if matchs.sum() == 0:
            continue
        else:
            prec = matchs[0]/matchs.sum()
            sens = matchs[0]/float(tot)
            M.append([iden, prec, sens])

    M = np.array(M)
    N = np.array(N)
    return M, N


def draw_all(args):
    A = []
    sys.stdout.write('Summary of first 100 hits\n')
    for arg in args:
        N, label, ls, color = arg
        arg_e = np.argsort(N[:,0])
        Prec = N[:,1][arg_e]
        Sens = N[:,2][arg_e]
        A.append([Prec,Sens,label,color,ls])
    

        Prec_sum = Prec[:100].sum()
        Sens_sum = Sens[:100].sum()
        valid = Prec[:100] != 0.0
        Cover_sum = (Sens[:100][valid]/Prec[:100][valid]).sum()
        sys.stdout.write('%s total_hit_n=%s sum(Pres)=%s sum(Sens)=%s sum(Cover)=%s\n' % (label, N.shape[0], Prec_sum, Sens_sum,Cover_sum))
    draw(A,'./fig/prec_sens.png')    



if __name__ == '__main__':
    ref_aln_dir = './miqs_yamada_benchmark_dataset/dali_reference'
    ref_aln_ext = '.dal'
    bench = './repo/cath20-scop'

    d = [
        ['./repo/pw/ssearch_MIQS_cath20-scop','SSEARCH MIQS',1,'grey'],
        ['./repo/pw/ssearch_BL62_cath20-scop','SSEARCH BL62',1,'black'],
        ['./repo/pw/blastp_BL62_cath20-scop','BLASTP BL62',1,'g'],
        ['./repo/pw/last10000_MIQS_cath20-scop','LAST MIQS',1,'r'],
        ['./repo/pw/last10000_BL62_cath20-scop','LAST BL62',1,'b'],
        ['./repo/pw/csblast_NA_cath20-scop','CS-BLAST',1,'orange'],
    ]

    A = []
    for each, label, ls, col in d:
        M,N = cal(pathlim.Files(each,exts=['.seq']).get_paths(), bench, ref_aln_dir, ref_aln_ext)        
        A.append([N,label,ls,col])
    draw_all(A)
    
