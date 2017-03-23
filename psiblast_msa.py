from parse_aln import *
import subprocess
from seqlim import Seq
from pypsi import PYPSI
import filecmp

class PSIBLAST(PYPSI):
    psiblast_path = "/home/klim/local/build_all/ReleaseMT_orig/bin/psiblast"

    def go(self, query, subject, outfile, pssm_db=None, n_iter=1, evalue=10, evalue_pssm=0.1, pssm=None, mat='BL62', continuous=False):
        self.query = query
        self.subject = subject
        self.pssm = pssm
        self.evalue = evalue
        self.evalue_pssm = evalue_pssm
        self.mat = mat

        cur_iter = 1 
        if n_iter == 1:
            output = self.go_psiblast(self.query, self.subject, evalue=self.evalue, matrix=self.mat)

        elif n_iter > 1:
            if pssm_db:
                self.pssm_db = pssm_db
            else:
                self.pssm_db = self.subject
            
            ori_seqob = Seq.parse(open(self.query))
            self.query_seq = ori_seqob[0].seq
            self.query_tag = ori_seqob[0].tag

            #iteration 1 againt pssm db
            output = self.go_psiblast(self.query, self.pssm_db, evalue=self.evalue_pssm, matrix=self.mat)
            if continuous:
                cur_output = self.go_psiblast(self.query, self.subject, evalue=self.evalue, matrix=self.mat)
                open(outfile+'.'+str(cur_iter)+'.hits','w').write(cur_output[0])
    
            prvs_msa_file = None
            msa_file = None
            for i in range(1, n_iter, 1):
                cur_iter += 1

                #merge
                msa = make_msa(self.query_tag, self.query_seq, output[0], insert_ori=True)

                msa_len = len(msa)
                if msa_len == 1:
                    break
                if msa_file:
                    prvs_msa_file = msa_file
                msa_file = outfile+'.'+str(i)+'.fasta'
                msa.write(open(msa_file,'w'))

                if prvs_msa_file and filecmp.cmp(msa_file, prvs_msa_file):
                    print('converged')
                    break

                if i < n_iter-1:

                    output = self.go_psiblast(msa_file, self.pssm_db, evalue=self.evalue_pssm, matrix=self.mat, is_query_msa=True)
                    if continuous:
                        cur_output = self.go_psiblast(msa_file, self.subject, evalue=self.evalue, matrix=self.mat, is_query_msa=True)
                        open(outfile+'.'+str(cur_iter)+'.hits','w').write(cur_output[0])

            #final iteration
            if msa_file:    
                output = self.go_psiblast(msa_file, self.subject, evalue=self.evalue, matrix=self.mat, is_query_msa=True)
            else:
                output = self.go_psiblast(self.query, self.subject, evalue=self.evalue, matrix=self.mat)

        for i in range(cur_iter, n_iter+1, 1):
            open(outfile+'.'+str(i)+'.hits','w').write(output[0])

        return output

if __name__ == '__main__':
    s = PSIBLAST()
    s.go('./repo/scop20_2738/d1brta_.seq','./repo/blastdb/scop20_validation','./temp/temp.seq', n_iter=3, evalue=1.0, evalue_pssm=0.002, pssm_db='./repo/blastdb/scop20_validation', mat='BL62', continuous=True)
            


