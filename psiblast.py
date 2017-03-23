from parse_aln import *
import subprocess, os
from seqlim import Seq
from pypsi import PYPSI
import filecmp

class PSIBLAST(PYPSI):
    psiblast_path = "/home/klim/local/build_all/ReleaseMT_orig/bin/psiblast"
    
    def go(self, query, subject, outfile, cutoff=0.3, pssm_db=None, n_iter=1, evalue=10, evalue_pssm=0.1, pssm=None, mat='BL62'):
        self.query = query
        self.subject = subject
        self.pssm = pssm
        self.evalue = evalue
        self.evalue_pssm = evalue_pssm
        self.mat = mat

        if n_iter == 1:
            output = self.go_psiblast(self.query, self.subject, evalue=self.evalue, matrix=self.mat)

        elif n_iter > 1:
            if pssm_db:
                self.pssm_db = pssm_db
            else:
                self.pssm_db = self.subject

            pssm_file = outfile+'.'+str(n_iter)+'.pssm'

            if os.path.isfile(pssm_file):
                os.remove(pssm_file)

            self.go_psiblast(self.query, self.pssm_db, evalue=self.evalue, evalue_pssm=evalue_pssm, matrix=self.mat, n_iter=n_iter, out_pssm=pssm_file)
            #final iteration
            if os.path.isfile(pssm_file):
                output = self.go_psiblast(pssm_file, self.subject, evalue=self.evalue, matrix=self.mat, is_query_pssm=True)
            else:
                output = self.go_psiblast(self.query, self.subject, evalue=self.evalue, matrix=self.mat)

        open(outfile,'w').write(output[0])
        return output


if __name__ == '__main__':
    s = PSIBLAST()
    s.go('./repo/scop20_validation/d1gcya1.seq','./repo/blastdb/scop20_validation','./temp/temp.seq', n_iter=5, evalue=1.0, evalue_pssm=100.0, pssm_db='./repo/blastdb/uniref50', mat='BL62')
            


