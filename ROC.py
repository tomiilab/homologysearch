import os

import numpy as np
import parse_aln as pa

from collections import Counter
from seqlim import MSeq
from path import Path
    
class Protein:
    def __init__(self, id_, family):
        self.id = id_
        self.family = family
        self.superfamily, self.fold = self.parse_family(self.family) 

    @classmethod
    def parse_family(cls, g):
        if g[0].isdigit():
            sfam = '.'.join(g[:4])
            fold = '.'.join(g[:3])
        else:
            sfam = '.'.join(g[:3])
            fold = '.'.join(g[:2])
        return sfam, fold


class Proteins(dict):
    def __init__(self, file_path, standard, JG_path=None):
        self.file_path = file_path
        self.name = os.path.splitext(os.path.basename(self.file_path))[0]
        self.JG_path = JG_path

        self.initial_seqnum = 0
        sfam = []
        fold = [] 
        for o in MSeq.parse(open(self.file_path)):
            basename, cat = o.tag.split()[:2]
            family = cat.split('.')
            protein = Protein(basename, family)
            protein.seqlen = len(o.seq)
            self[protein.id] = protein
            self.initial_seqnum += 1
            sfam.append(protein.superfamily)
            fold.append(protein.fold)

        self.standard = standard
        self.valid_seqnum = 0
        self.id2size = {}
        
        print("standard: %s" % self.standard)
        
        if self.standard == "superfamily":
            superfamily2size = Counter(sfam)
            for id_, each in self.items():
                n = superfamily2size[each.superfamily]
                each.size = n
                self.id2size[each.id] = n
                if n > 1:
                    self.valid_seqnum += 1
        elif self.standard == "fold":
            fold2size = Counter(fold)
            for id_, each in self.items():
                n = fold2size[each.fold]
                each.size = n
                self.id2size[each.id] = n
                if n > 1:
                    self.valid_seqnum += 1
        elif self.standard == "JG":
            for l in open(self.JG_path+'/count/'+self.name):
                cs = l.split()
                n = int(cs[3])
                self.id2size[cs[0]] = n
                if n > 1:
                    self.valid_seqnum += 1
            
            for id_, each in self.items():
                each.size = self.id2size[id_]

            
            self.pair2JG_judge = {}
            for l in open(self.JG_path+'/pair/'+self.name):
                cs = l.split()
                pair = [cs[1], cs[2]]
                pair.sort()
                pair_str = ','.join(pair)
                self.pair2JG_judge[pair_str] = cs[0]

        print("benchmark: %s" % self.name)        
        print("valid seqnum: %s" % self.valid_seqnum)        
            

    def parse_hits(self, f, max_evalue=10.0):
        f = Path(f)
        query_name = f.basename().split('.seq')[0]
        query = self[query_name]

        if query.size == 1:
            return None
        
        
        lines = open(f).read().splitlines()
        
        if not lines:
            print 'No content: '+f
            return None


        aln_num = 0
        hits = Hits()

        tids= set()
        for each in lines[::-1]:
            tid, e = each.split()
            e = float(e)
            tid = tid.split()[0].replace('lcl|','')
            if query.id == tid:
                continue
            if e > max_evalue:
                continue
            if tid in tids:
                continue
            else:
                tids.add(tid)
            target = self[tid]
            hits.append(Hit(query, target, e, standard=self.standard, db=getattr(self, "pair2JG_judge", None)))

        hits.seqlen = query.seqlen
        hits.id = query.id
        
        hits.superfamily = query.superfamily
        hits.fold = query.fold
        hits.family = query.family
        hits.size = query.size

        return hits



class Hit:
    #match search result, ignore false negatives according to given conditions 
    def __init__(self, query, target, e, standard, db=None):
        self.evalue = e
        self.query = query 
        self.target = target
        self.standard = standard

        if self.standard == "JG":
            pair = [query.id, target.id]
            pair.sort()
            judge = db.get(','.join(pair))
            if judge == 'A':
                self.bin = True
            elif judge == 'B':
                self.bin = None
            else:
                self.bin = False

        elif self.standard == "superfamily":
            if self.query.superfamily == self.target.superfamily:
                self.bin = True
            elif self.query.fold != self.target.fold:
                self.bin = False
            else:
                self.bin = None

        elif self.standard == "fold":
            if self.query.fold == self.target.fold:
                self.bin = True
            else:
                self.bin = False


class Hits(list):
    pass



def first_passage_index(arr, cutoff):
    ind = np.where(arr > cutoff)
    if np.any(ind) is not None and len(ind[0]):
        return ind[0][0]

class Benchmark(list):
    def __init__(self, L):
        self.extend(L)
        self.chain_and_arrayify()

    def chain_and_arrayify(self):
        self.bin = []
        self.evalue = []
        self.weight = []
        self.superfamily = []
        self.fold = []
        self.size = []

        for hits in self:
            size = hits.__dict__['size']
            for e in hits:
                bin = e.__dict__['bin']
                if bin is None:
                    continue
                self.bin.append(bin)
                self.evalue.append(e.evalue)
                self.weight.append(1.0/(size - 1))
                self.size.append(size)
                self.superfamily.append(e.query.superfamily)
                self.fold.append(e.query.fold)

        self.evalue = np.array(self.evalue)
        arg = self.evalue.argsort()

        self.evalue = self.evalue[arg]
        self.bin = np.array(self.bin)[arg].astype('bool') #True or False
        self.weight = np.array(self.weight)[arg]
        self.superfamily = np.array(self.superfamily)[arg]
        self.fold = np.array(self.fold)[arg]
        self.size = np.array(self.size)[arg]
        

    def _cal_ROC(self, bin, evalue=None, weight=None, fdr=None, max_evalue=None, max_FP_count=None):
        if np.any(bin) is None:
            return None

        limit = len(bin)
        if not limit:
            return None
        if max_evalue is not None:
            assert evalue is not None

            idx = first_passage_index(evalue, max_evalue)
            if idx == 0:
                return None
            if idx is not None:
                limit = idx

        if weight is None:
            TP = np.ones(limit, dtype='float')
            FP = np.ones(limit, dtype='float')
        else:
            TP = weight[:limit].copy()
            FP = weight[:limit].copy()

        TP[~bin[:limit]] = 0.0
        FP[bin[:limit]] = 0.0
        TP = np.cumsum(TP)
        FP = np.cumsum(FP)

        if max_FP_count is not None:
            limit = first_passage_index(FP, max_FP_count)

        TP = TP[:limit]
        FP = FP[:limit]
        
        evalue_fdr = None
        TP_fdr = None
        if fdr is not None:
            FDR = FP/(TP + FP)
            FDR_diff = FDR - fdr
            arg_close_left = first_passage_index(FDR, fdr)
            if arg_close_left is not None:
                evalue_fdr = evalue[arg_close_left]
                TP_fdr = TP[arg_close_left]
            else:
                evalue_fdr = evalue[limit-1]
                TP_fdr = TP[limit-1]

        return {
            'TP':TP,
            'FP':FP,
            'evalue_fdr':evalue_fdr,
            'TP_fdr':TP_fdr,
        }

    def cal_ROC(self, superfamily=None, fold=None, weighted=True, fdr=0.05):
        
        wfp = 0.0
        wtp = 0.0
        X = []
        Y = []

        e = 0.0
        t = True
        h = True
        self.evalue_fdr = None

        if superfamily:
            arg = self.superfamily == superfamily
            bin = self.bin[arg]
            w = self.weight[arg]
            evalues = self.evalue[arg]
        elif fold:
            arg = self.fold == fold
            bin = self.bin[arg]
            w = self.weight[arg]
            evalues = self.evalue[arg]
        else:
            bin = self.bin
            w = self.weight
            evalues = self.evalue
 
 
        if weighted:
            r = self._cal_ROC(bin, evalue=evalues, weight=w, fdr=fdr)
        else:
            r = self._cal_ROC(bin, evalue=evalues, weight=None, fdr=fdr)
        
        self.FP = r['FP']
        self.TP = r['TP']
        self.evalue_fdr = r['evalue_fdr']
        self.TP_fdr = r['TP_fdr']
        a = evalues < self.TP_fdr

        print(np.count_nonzero(bin[a]))
        print(self.FP[-1])
        print(self.TP[-1])
        print("evalue fdr: %s " % r['evalue_fdr'])
        print("TP fdr: %s " % r['TP_fdr'])

        '''
        passwtp = True

        for i, each in enumerate(bin):
            if each:
                if w is None:
                    wtp += 1
                else:
                    wtp += w[i]
            else:
                if w is None:
                    wfp += 1
                else:
                    wfp += w[i]

            Y.append(wtp)
            X.append(wfp)
            if t and evalues[i] > 0.1:
                print('E-value=0.1 wtp=%s' % wtp)
                t =False
            if passwtp and wfp > 40:
                print('wfp=40 e=%s' % evalues[i])
                passwtp =False
            if self.evalue_fdr is None and wfp/(wtp+wfp) > fdr:
                self.evalue_fdr =  evalues[i]
                print('FDR=%.3f %.3f %.1f' % (fdr, evalues[i], wtp))
        '''

        if self.evalue_fdr == None:
            self.evalue_fdr = evalues[i]

    
    def cal_ROCN(self, max_FP_count=5, max_evalue=1.0):
        rocs = []
        self.ROCN_benchmarks = []
        for E in self:
            total_TP_count = E.size
            if total_TP_count < 2:
                print('ignore')
                continue

            bench2 = Benchmark([E])
            r = bench2._cal_ROC(bench2.bin, evalue=bench2.evalue, max_evalue=max_evalue, max_FP_count=max_FP_count)
            if r is not None:
                TP_sum = r['TP'][~bench2.bin[:len(r['TP'])]].sum() + r['TP'][-1]*(max_FP_count-r['FP'][-1])
                bench2.ROCN = TP_sum / float(max_FP_count*(total_TP_count-1)) #total_TP_count includes self
                if bench2.ROCN > 1:
                    print(TP_sum, total_TP_count)
            else:
                TP_sum = 0
                bench2.ROCN = 0
            self.ROCN_benchmarks.append(bench2)


    def draw(self, ax, fdr=0.05, style={}):
        
        ax.plot(self.FP, self.TP, **style)
        ax.plot(np.arange(0, 4000, (1.0/fdr)-1),c='black')





