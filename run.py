import subprocess, time
import argparse, subprocess, os, job, re, time

from psiblast import PSIBLAST
from psiblast_msa import PSIBLAST as PSIBLASTmsa
from psiblast_msa_fil import PSIBLAST as PSIBLASTmsafil
from psiblast_msa_filall import PSIBLAST as PSIBLASTmsafilall
from psiblast_block import PSIBLAST as PSIBLASTblock

from pypsi import *
from path import Path
from seqlim import MSeq
import parse_aln as pa


class PSIBLASToda5(PSIBLAST):
    psiblast_path = "/home/klim/local/build_all/ReleaseMT_5/bin/psiblast"

class PSIBLASToda9(PSIBLAST):
    psiblast_path = "/home/klim/local/build_all/ReleaseMT_9/bin/psiblast"

class PSIBLASToda13(PSIBLAST):
    psiblast_path = "/home/klim/local/build_all/ReleaseMT_13/bin/psiblast"

class PSIBLASToda17(PSIBLAST):
    psiblast_path = "/home/klim/local/build_all/ReleaseMT_17/bin/psiblast"

class PSIBLASToda21(PSIBLAST):
    psiblast_path = "/home/klim/local/build_all/ReleaseMT_21/bin/psiblast"

class PSIBLASToda25(PSIBLAST):
    psiblast_path = "/home/klim/local/build_all/ReleaseMT_25/bin/psiblast"

class PSIBLASToda29(PSIBLAST):
    psiblast_path = "/home/klim/local/build_all/ReleaseMT_29/bin/psiblast"

class PSIBLASToda41(PSIBLAST):
    psiblast_path = "/home/klim/local/build_all/ReleaseMT_41/bin/psiblast"

class PSIBLASToda9msa(PSIBLASTmsa):
    psiblast_path = "/home/klim/local/build_all/ReleaseMT_9/bin/psiblast"

class PSIBLASToda13msa(PSIBLASTmsa):
    psiblast_path = "/home/klim/local/build_all/ReleaseMT_13/bin/psiblast"

class PSIBLASToda17msa(PSIBLASTmsa):
    psiblast_path = "/home/klim/local/build_all/ReleaseMT_17/bin/psiblast"

class PSIBLASToda21msa(PSIBLASTmsa):
    psiblast_path = "/home/klim/local/build_all/ReleaseMT_21/bin/psiblast"

class PSIBLASToda25msa(PSIBLASTmsa):
    psiblast_path = "/home/klim/local/build_all/ReleaseMT_25/bin/psiblast"

class PSIBLASToda25msafil(PSIBLASTmsafil):
    psiblast_path = "/home/klim/local/build_all/ReleaseMT_25/bin/psiblast"

class PSIBLASTfullmsa(PSIBLASTmsa):
    psiblast_path = "/home/klim/local/ncbi-blast-2.5.0+psiblast_fullhenikoff/bin/psiblast"


famat2blmat = {
    'BL50':'BLOSUM50',
    'BL62':'BLOSUM62',
}

def ss_mat(mat):
    if mat == 'MIQS':
        return './data/'+mat+'.txt'
    else:
        return mat

def bl_mat(mat):
    g = famat2blmat.get(mat)
    if g:
        return g
    else:
        return 'BLOSUM62'


class PSIBLAST_run:
    def __init__(self, fs, method, matrix, query, subject, db, n_iter, evalue_pssm):
        self.out_dir = './repo/raw/'+query+'/'+'_'.join([method, matrix, subject])
        self.temp_dir = './temp/'+'_'.join([method, matrix, query, subject])

        if n_iter > 1:
            self.out_dir+='_'+str(n_iter)
            self.temp_dir+='_'+str(n_iter)
        if evalue_pssm:
            self.out_dir += '/'+str(evalue_pssm)
            self.temp_dir += '/'+str(evalue_pssm)
        if db:
            db_name = db.split('/')[-1]
            self.out_dir += '/'+db_name
            self.temp_dir += '/'+db_name
        
        Path(self.out_dir).makedirs_p()
        Path(self.temp_dir).makedirs_p()

        self.files = fs
        self.method = method
        self.matrix = matrix
        self.subject = subject
        self.db = db
        if not '/' in self.subject:
            self.subject = './repo/blastdb/'+self.subject
        if not '/' in self.db:
            self.db = './repo/blastdb/'+self.db

        self.n_iter = n_iter
        self.evalue_pssm = float(evalue_pssm)

    def _run(self, func):
        d = {
            'mat':self.matrix,
            'evalue':10,
        }
        if self.n_iter is not None:
            d['n_iter'] = self.n_iter
        if self.db is not None:
            d['pssm_db'] = self.db
        if self.evalue_pssm is not None:
            d['evalue_pssm'] = self.evalue_pssm
        if '_msa' in self.method or '_arelia' in self.method:
            d['continuous'] = True

        st =os.times()
        for f in self.files:
            basename = f.split('/')[-1]
            out_file = self.out_dir+'/'+basename
            caller = func()
            r = caller.go(
                f,
                self.subject,
                out_file,
                **d
            )
        et =os.times()
        exe_time = et[2]+et[3]-st[2]-st[3]
        return exe_time


        #parsed = pa.filter(r[0])
        #print(parsed)
        #open(out_file,'w').write(parsed)

    def run(self):
        if self.method == 'psiblast':
            self._run(PSIBLAST)
        elif self.method == 'psiblast_msa':
            self._run(PSIBLASTmsa)
        elif self.method == 'psiblast_block':
            self._run(PSIBLASTblock)
        elif self.method == 'psiblast_msa_fil':
            self._run(PSIBLASTmsafil)
        elif self.method == 'psiblastfull_msa':
            self._run(PSIBLASTfullmsa)
        elif self.method == 'psiblast_msa_filall':
            self._run(PSIBLASTmsafilall)
        elif self.method == 'psiblastoda5':
            self._run(PSIBLASToda5)
        elif self.method == 'psiblastoda9':
            self._run(PSIBLASToda9)
        elif self.method == 'psiblastoda13':
            self._run(PSIBLASToda13)
        elif self.method == 'psiblastoda17':
            self._run(PSIBLASToda17)
        elif self.method == 'psiblastoda21':
            self._run(PSIBLASToda21)
        elif self.method == 'psiblastoda25':
            self._run(PSIBLASToda25)
        elif self.method == 'psiblastoda29':
            self._run(PSIBLASToda29)
        elif self.method == 'psiblastoda41':
            self._run(PSIBLASToda41)
        elif self.method == 'psiblastoda21_msa':
            self._run(PSIBLASToda21msa)
        elif self.method == 'psiblastoda25_msa':
            self._run(PSIBLASToda25msa)
        elif self.method == 'psiblastoda25_msa_fil':
            self._run(PSIBLASToda25msafil)


def run_search(fs, method, matrix, query, subject, db, n_iter, evalue_pssm):
    
    out_dir = './repo/raw/'+query+'/'+'_'.join([method, matrix, subject])
    temp_dir = './temp/'+'_'.join([method, matrix, query, subject])

    if n_iter > 1:
        out_dir+='_'+str(n_iter)
        temp_dir+='_'+str(n_iter)
    if evalue_pssm:
        out_dir += '/'+str(evalue_pssm)
        temp_dir += '/'+str(evalue_pssm)
    if db:
        out_dir += '/'+db
        temp_dir += '/'+db

        
    Path(out_dir).makedirs_p()
    Path(temp_dir).makedirs_p()

    
    st = os.times()

    for f in fs:
        basename = f.split('/')[-1]
        out_file = out_dir+'/'+basename

        if method[:5] == 'lastm':
            s = MSeq.parse(open(f))
            m = 10**float(method[5:])
            proc = PYPSI.go_last(f, './repo/lastmdb/'+subject, query_len=len(s[0].seq), matrix=matrix, m=int(m), evalue=2000)
    
        elif method[:4] == 'last':
            s = MSeq.parse(open(f))
            m = 10**float(method[4:])
            proc = PYPSI.go_last(f, './repo/lastdb/'+subject, query_len=len(s[0].seq), matrix=matrix, m=int(m), evalue=2000)
            
        if method == 'blastp':
            proc = subprocess.Popen(['blastp','-query',f,'-db','./repo/blastdb/'+subject,'-evalue','1000','-num_alignments','100000000'], stdout=subprocess.PIPE).communicate()
        
        if method == 'blastpgp':
            proc = subprocess.Popen(['blastpgp','-i',f,'-d','./repo/blastpgpdb/'+subject,'-e','20'], stdout=subprocess.PIPE).communicate()
       
        if method == 'ssearch':
            cmd = ['ssearch36','-m','B','-E','1000','-s',ss_mat(matrix), f,'./repo/fastadb/'+subject]
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE).communicate()

        if method == 'cs-blast':
            new_f2 = './repo/blastpgpdb/'+subject
            cmd = ['csblast', '-i' ,f, '-d',new_f2,'-e','1000','-b','10000000','-D', '/home/klim/local/csblast-2.2.3_linux64/data/K4000.crf', '--blast-path','/home/klim/local/blast-2.2.26/bin/']
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE).communicate()

        if proc:
            open(out_file,'w').write(pa.filter(proc[0]))


    et =os.times()
    exe_time = et[2]+et[3]-st[2]-st[3]
    return exe_time

if __name__ == '__main__':
    par = argparse.ArgumentParser()
    par.add_argument('method')
    par.add_argument('matrix')
    par.add_argument('query')
    par.add_argument('subject')
    par.add_argument('-db', default=None)
    par.add_argument('-n_iter', default=1, type=int)
    par.add_argument('-evalue_pssm',type=float,default=None)
    args = par.parse_args()

    query = './repo/'+args.query
    if os.path.isdir(query):
        files = Path(query).walkfiles(pattern='*.seq')
    else:
        files = [query]

    if "psiblast" in args.method:
        run = PSIBLAST_run(files, args.method, args.matrix, args.query, args.subject, args.db, args.n_iter, args.evalue_pssm)
        exe_time = run.run()
    else:
        exe_time = run_search(files, args.method, args.matrix, args.query, args.subject, args.db, args.n_iter, args.evalue_pssm)

    open('./repo/'+args.method+'_'+args.matrix+'_'+args.query+'_'+args.subject+'_'+str(args.db)+'_'+str(args.n_iter)+'_'+str(args.evalue_pssm)+'.txt','w').write(str(exe_time))
            

