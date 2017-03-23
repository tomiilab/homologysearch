import re
import csv
import numpy as np
from seqlim import MSeq
from collections import defaultdict
from path import Path

#bl_pa = re.compile(r"Score \= {1,4}(\w.+) bits \(\d+\),  Expect \= +([\.\deE\-\+]+)[^\n]*\n Identities \= \d+\/(\d+) \((\w+)\%\), Positives \= .+ \(\S+\), Gaps \= .+\(\S+\)\n\n(.+)$",re.DOTALL)
bl_pa = re.compile(r"Score \= +([\d\.]+) bits \(\d+\), +Expect \= +([\.\deE\-\+]+)[^\n]*\n Identities +\= +\d+\/(\d+) \((\w+)\%\), +Positives +\= .+ \(\S+\)[^\n]*\n\n(.+)$",re.DOTALL)

footer_bl_pa = re.compile(r"\nLambda|\n  Database\:")
footer_ssbl_pa = re.compile(r"\n\n\n\n")
tag_bl_pa = re.compile(r"\n *Length *\=")


def classify(seqtxt):
    line = seqtxt.split('\n',1)[0]
    if seqtxt.startswith('BLAST') or seqtxt.startswith('PSIBLAST'):
        return 'blast'
    elif seqtxt.startswith('# LAST'):
        return 'maf'
    elif 'ssearch' in line or '\n# ssearch' in seqtxt:
        return 'ssearch blast'

def filter(seqtxt):
    cl = classify(seqtxt)
    if cl == 'blast':
        return fil_blast(seqtxt)
    elif cl == 'maf':
        return fil_maf(seqtxt)
    elif cl == 'ssearch blast':
        return fil_blast(seqtxt)


def fil_blast(seqtxt):
    o = ''
    accept = True
    for l in seqtxt.splitlines():
        if l and l[0] == '>':
            tag = l[1:].strip()
            if tag.startswith('UniRef') or tag.startswith('pfam21') or tag.startswith('up_'):
                accept = False
            else:
                accept = True
        if accept:
            o += l+'\n'
        elif l.startswith('Results from round'):
            o += l+'\n'
            accept = True

    return o

def fil_maf(seqtxt):
    o = ''
    a =''
    s = []
    for l in seqtxt.splitlines():
        if not l:
            continue
        if l[:2] == 'a ':
            a = l
            s = []
        elif l[:2] == 's ':
            s.append(l)
        elif l[0] == '#':
            o += l+'\n'
            
        if len(s) == 2:
            if s[0].startswith('s UniRef') or s[0].startswith('s pfam21') or s[0].startswith('s up_'):
                s = []
                continue
            o += '\n'.join([a]+s)+'\n\n'
    return o


def parse(seqtxt, max_aln_num=float('inf')):
    cl = classify(seqtxt)
    if cl == 'blast':
        return parse_blast(seqtxt, regex=bl_pa, max_aln_num=max_aln_num)
    elif cl == 'maf':
        return parse_maf(seqtxt, max_aln_num=max_aln_num)
    elif cl == 'ssearch blast':
        return parse_blast(seqtxt, no_hit_mess="!! No sequences with E()",
            sig_mess ='The best scores are:',
            footer_start=footer_ssbl_pa, regex= bl_pa, max_aln_num=max_aln_num
            )

def parse_blast(seqtxt, no_hit_mess="***** No hits found *****",
                sig_mess ='Sequences producing significant alignments:',
                footer_start=footer_bl_pa, regex= bl_pa, max_aln_num=float('inf')
            ):

    r = []
    if no_hit_mess in seqtxt:
        return r
    round_mess = 'Results from round'
    if round_mess in seqtxt:
        seqtxt = seqtxt.split(round_mess)[-1]
    if sig_mess in seqtxt:
        seqtxt = seqtxt.split(sig_mess)[1]
    block = re.split(footer_start, seqtxt)[0]

    if '\n\n>' in block:
        blocks = block.strip().split("\n\n>")[1:]
    else:
        blocks = [block]
    qid = 'query'
    aln_num = 0
    for block in blocks:
        if aln_num > max_aln_num:
            break
        if not block:
            continue
        subblocks = block.split('\n\n ')
        tid = ' '.join(re.split(tag_bl_pa,subblocks[0])[:-1])
        tid = tid.replace('\n',' ')
        sub = []
        for subblock in subblocks[1:]:
            if not subblock:
                continue
            g = re.search(regex, subblock.strip())
            evalue = g.group(2)
            if evalue[0] == 'e':
                evalue = '1'+evalue
            evalue = float(evalue)
            align_len = int(g.group(3))
            iden = float(g.group(4))
            pairs = g.group(5).strip('\n').split('\n\n')
            qseq = ''
            tseq = ''
            count = 0
            pair_num = len(pairs)
            pos = ''
            err = False
            for pair in pairs:
                lines = pair.splitlines()
                if len(lines) < 3:
                    print('error')
                    print(pair)
                    err =True
                    break
                cs = lines[0].split()
                if len(cs) > 2:
                    if count == 0:
                        qs = int(cs[1])
                    if count == pair_num-1:
                        qe = int(cs[-1])
                    qseq += cs[2]
                elif len(cs) == 2:
                    qseq += cs[1]
                    
                cs = lines[2].split()
                if len(cs) > 2:
                    if count == 0:
                        ts = int(cs[1])
                    if count == pair_num-1:
                        te = int(cs[-1])
                    tseq += cs[2]
                elif len(cs) == 2:
                    tseq += cs[1]
            
                #seqidx = lines[2].index(tseq)
                #pos += lines[1][seqidx:seqidx+len(cs[2])]
    
                count += 1
            if not err:
                sub.append(((qid, qseq, qs, qe),(tid, tseq, ts, te), evalue))
        if sub:
            r.append(sub)
            aln_num += 1
    return r

def parse_maf(seqtxt, max_aln_num=float('inf')):
    r = []
    is_target = True
    tid2alns = defaultdict(list)
    tids = []
    aln_num = 0
    for b in seqtxt.split('\na ')[1:]:
        if aln_num > max_aln_num:
            break
        cs = b.strip().split()
        if len(cs) < 17:
            continue
        e = float(cs[2].split('=')[1])
        score = int(cs[0].split('=')[1])

        tid = cs[4]
        ts = int(cs[5])+1
        te = ts -1 +int(cs[6])
        tseq = cs[9]
        
        qid = cs[11]
        qs = int(cs[12])+1
        qe = qs -1 +int(cs[13])
        qseq = cs[16]
        tid2alns[tid].append([(qid,qseq,qs,qe),(tid,tseq,ts,te), e])

        if not tid in tids:
            tids.append(tid)
            aln_num += 1

    for tid in tids:
        r.append(tid2alns[tid])
    return r

def make_msa(oritag, oriseq, seqtxt, insert_ori=True):
    r = parse(seqtxt)
    seqob = MSeq()
    seqs = []
    seq_len = len(oriseq)
    if insert_ori:
        seq = oriseq.upper().replace('J','X').replace('Z','X').replace('B','X')
        seqob.add(oritag, seq)
        seqs.append(seq)
    if not r:
        return seqob
    for i,sub in enumerate(r):
        arr = np.zeros(seq_len,dtype='S1')
   
        for j, each in enumerate(sub):
            q, t, e = each[:3]
            if j >0:
                print(e)
            qid, qseq, sq, eq = q
            tid, tseq, st, et = t
            SEQS = np.array([[qseq],[tseq]]).view('S1')
            gaparg = SEQS == '-'

            #if overlap pass
            if np.any(arr[sq-1:eq] != ''):
                continue
            arr[sq-1:eq] = SEQS[1][~(gaparg[0])]
        arr[arr=='']='-'
        seq = arr.view('S'+str(seq_len))[0].upper().replace('J','X').replace('Z','X').replace('B','X')
        if not seq in seqs:
            seqob.add(tid.strip()+'|'+str(e), seq)
            seqs.append(seq)
    return seqob

def make_msa_fil_gap(oritag, oriseq, seqtxt, insert_ori=True, fil_all=False):
    r = parse(seqtxt)
    seqob = MSeq()
    seqs = []
    seq_len = len(oriseq)
    oriseq_ = np.array([oriseq.upper().replace('J','X').replace('Z','X').replace('B','X')]).view('S1')
    if insert_ori:
        seq = oriseq.upper().replace('J','X').replace('Z','X').replace('B','X')
        seqob.add(oritag, seq)
        seqs.append(seq)
    if not r:
        return seqob
    for i,sub in enumerate(r):
        arr = np.zeros(seq_len,dtype='S1')
   
        for j, each in enumerate(sub):
            q, t, e = each[:3]
            if j >0:
                print(e)
            qid, qseq, sq, eq = q
            tid, tseq, st, et = t
            SEQS = np.array([[qseq],[tseq]]).view('S1')
            gaparg = SEQS == '-'

            #if overlap pass
            if np.any(arr[sq-1:eq] != ''):
                continue
            if fil_all:
                arr[:] = oriseq_[:]
            else:
                arr[sq-1:eq] = oriseq_[sq-1:eq]
            arr[sq-1:eq][~gaparg[1][~gaparg[0]]] = SEQS[1][~(gaparg[0])][~gaparg[1][~gaparg[0]]]
        
        arr[arr=='']='-'
        seq = arr.view('S'+str(seq_len))[0].upper().replace('J','X').replace('Z','X').replace('B','X')
        if not seq in seqs:
            seqob.add(tid.strip()+'|'+str(e), seq)
            seqs.append(seq)
    return seqob



def summarize_hits(f):
    f = Path(f)
    cs = f.basename().split('.')
    query = cs[0]

    outfile = f.dirname() / (f.basename()+'.tsv')
    alns = parse(open(f).read())
    
    if not alns:
        print 'No content: '+f
        return None
    print(outfile)
    with open(outfile, 'w') as csvfile:
        writer = csv.writer(csvfile, delimiter='\t')

       
        for each in alns:
            q, t, e = each[0][:3]

            tid, tseq, ts, te = t
            qid, qseq, qs, qe = q

            qid = qid.split()[0].replace('lcl|','')
            tid = tid.split()[0].replace('lcl|','')

            if tid.startswith('UniRef'):
                continue
            if tid.startswith('pfam'):
                continue
            if tid.startswith('up_'):
                continue
           
            writer.writerow((tid, e)) 

def summarize_ssearch_summary(f):
    f = Path(f)
    cs = f.basename().split('.')
    query = cs[0]

    outfile = f.dirname() / (f.basename()+'.tsv')

    r = []
    started = False    

    with open(f) as ifh:
        for l in ifh:
            if l.startswith("The best scores are:"):
                started = True
                continue
            if started:
                cs = l.strip().split()
                if cs:
                    e = cs[-1]
                    tid = cs[0]
                    tid = tid.replace('lcl|','')
                    if tid.startswith('UniRef'):
                        continue
                    if tid.startswith('pfam'):
                        continue
                    if tid.startswith('up_'):
                        continue
                    r.append([tid, e])
                else:
                    break

 
    if not r:
        print 'No content: '+f
        return None

    print(outfile)
    with open(outfile, 'w') as csvfile:
        writer = csv.writer(csvfile, delimiter='\t')
        for row in r:
            writer.writerow(row) 



if __name__ == '__main__':
    import argparse
    par = argparse.ArgumentParser()
    par.add_argument('action')
    par.add_argument('path')
    args = par.parse_args()

    if args.action == "summarize_hits":
        for f in Path(args.path).walkfiles(pattern="*.seq"):
            summarize_hits(f)
        for f in Path(args.path).walkfiles(pattern="*.hits"):
            summarize_hits(f)
        for f in Path(args.path).walkfiles(pattern="*.res"):
            summarize_ssearch_summary(f)



