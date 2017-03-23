import subprocess

class PYPSI:
    famat2blmat = {
        'BL50':'BLOSUM50',
        'BL62':'BLOSUM62',
    }

    psiblast_path = 'psiblast'

    def __init__(self, ssearchcmd='ssearch36'):
        self.outputs = []
        self.ssearchcmd = ssearchcmd

    def blast_mat(self, mat):
        g = self.famat2blmat.get(mat)
        if g:
            return g
        else:
            return 'BL0SUM62'

    @classmethod
    def go_last(self, query, subject, evalue=10, m=10, query_len=300, pssm=None, matrix='BL62', is_query_pssm=False):
        D = float(query_len)/float(evalue)
        cmd = ['lastal', subject, query,'-D', str(D), '-p', matrix, '-m', str(m)]
        cmd = ['lastal', subject, query,'-e', '20', '-p', matrix, '-m', str(m)]
        if is_query_pssm:
            cmd.extend(['-Q', '5'])
        print ' '.join(cmd)    
        return subprocess.Popen(cmd, stdout=subprocess.PIPE).communicate()

    def go_ssearch(self, query, subject, evalue=10, pssm=None, mat='BL50'):
        cmd = [self.ssearchcmd, query, subject, '-s', mat, '-E', str(evalue), '-m', 'B']
        if pssm:
            cmd.extend(['-P', pssm+' 2'])
        print ' '.join(cmd)    
        return subprocess.Popen(cmd, stdout=subprocess.PIPE).communicate()

    @classmethod
    def go_blastp(self, query, subject, evalue=10, matrix='BL62',n_aln=10000000, is_subject_fasta=False, raise_error=True, quiet=False):
        cmd = ['blastp', '-evalue', str(evalue),'-num_alignments',str(n_aln)]
        cmd.extend(['-query',query])

        if is_subject_fasta:
            cmd.extend(['-subject',subject])
        else:
            cmd.extend(['-db',subject])
        
        if not quiet:
            print(' '.join(cmd))
        output = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
        if raise_error and not output[0]:
            raise OSError(output[1])

        return output

    @classmethod
    def go_psiblast(self, query, subject, evalue=10, evalue_pssm=None, out_pssm=None, out_ascii_pssm=None, n_iter=None,is_subject_fasta=False, is_query_msa=False, is_query_pssm=False, matrix='BL62',n_aln=10000000, raise_error=True, quiet=False):
        cmd = [self.psiblast_path, '-evalue', str(evalue),'-num_alignments',str(n_aln)]
        if is_query_msa:
            cmd.extend(['-in_msa',query])
        elif is_query_pssm:
            cmd.extend(['-in_pssm',query])
        else:
            cmd.extend(['-query',query])

        if is_subject_fasta:
            cmd.extend(['-subject',subject])
        else:
            cmd.extend(['-db',subject])
        if out_pssm:
            cmd.extend(['-out_pssm', out_pssm])
        if out_ascii_pssm:
            cmd.extend(['-out_ascii_pssm', out_ascii_pssm])
        if n_iter:
            cmd.extend(['-num_iterations', str(n_iter)])
        if evalue_pssm:
            cmd.extend(['-inclusion_ethresh', str(evalue_pssm)])
        
        if not quiet:
            print(' '.join(cmd))
        output = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
        if raise_error and not output[0]:
            raise OSError(output[1])

        return output

