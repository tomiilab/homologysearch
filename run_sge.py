from path import Path
import argparse, subprocess, os, job, re, time
from run import *

if __name__ == '__main__':
    class MYJOB(job.UGE):
    
        def main_job(self, *args): 
            par = argparse.ArgumentParser()
            par.add_argument('method')
            par.add_argument('matrix')
            par.add_argument('query')
            par.add_argument('subject')
            par.add_argument('subjobname')
            par.add_argument('-db', default=None)
            par.add_argument('-n_iter', default=1, type=int)
            par.add_argument('-evalue_pssm',type=float,default=None)
            par.add_argument('-nqueues', metavar='number_of_queues', default=5, type=int)
            par.add_argument('-node', nargs='+', default=[])
            par.add_argument('-max_concurrent', type=int)
            args = par.parse_args(args)

            self.args = args
            self.nqueues = args.nqueues
            self.subjobname = args.subjobname
            qs = []
            for node in args.node:
                qs.append('all.q@'+node)
            if qs:
                self.set_q(','.join(qs))
            

            self.max_concurrent = args.max_concurrent

            ids = []
            files = Path('./repo/'+args.query).walkfiles(pattern='*.seq')
            new_files = []
            for f in files:
                b = os.path.basename(f)
                if not b in ids:
                    new_files.append(f)
            files =  new_files

            print(len(files))
            return files, args.method, args.matrix, args.query, args.subject, args.db, args.n_iter, args.evalue_pssm

        def sub_job(self, files, method, matrix, query, subject, db, n_iter, evalue_pssm):
            if 'psiblast' in method:
                run = PSIBLAST_run(files, method, matrix, query, subject, db, n_iter, evalue_pssm)
                exe_time = run.run()
            else:
                exe_time = run_search(files, method, matrix, query, subject, db, n_iter, evalue_pssm)
            self.return_args([exe_time])
            
    qs = []
    
    for hostname in job.qhost_names()[2:-1]:
        qs.append('all.q@'+hostname)

    print(qs)
    myjob = MYJOB(UGE_options=['-q',','.join(qs)])
    myjob.go()

