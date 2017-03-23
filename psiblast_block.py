import os
from parse_aln import *
from pypsi import PYPSI

class PSIBLAST(PYPSI):
    psiblast_path = "/home/klim/local/ReleaseMT_print_block/bin/psiblast"
    
    def go(self, query, subject, outfile, cutoff=0.3, pssm_db=None, n_iter=1, evalue=10, evalue_pssm=0.1, pssm=None, mat='BL62'):
        self.query = query
        self.subject = subject
        self.pssm = pssm
        self.evalue = evalue
        self.evalue_pssm = evalue_pssm
        self.mat = mat

        pssm_output = None

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
            pssm_output = self.go_psiblast(self.query, self.pssm_db, evalue=evalue, evalue_pssm=self.evalue_pssm, matrix=self.mat, n_iter=n_iter, out_pssm=pssm_file)
            #final iteration
            if os.path.isfile(pssm_file):
                output = self.go_psiblast(pssm_file, self.subject, evalue=self.evalue, matrix=self.mat, is_query_pssm=True)
            else:
                output = self.go_psiblast(self.query, self.subject, evalue=self.evalue, matrix=self.mat)

       
        open(outfile,'w').write(pssm_output[0])
        open(outfile,'w').write('\n')
        open(outfile,'w').write(pssm_output[1])

        return pssm_output
 

