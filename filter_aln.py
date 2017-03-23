import argparse
import parse_aln as pa
from path import Path

parser = argparse.ArgumentParser()
parser.add_argument('dirpath', metavar='PATH', help='MSA file or dir containing MSA files')
args = parser.parse_args()

i = 0
print args.dirpath

for f in Path(args.dirpath).walkfiles(pattern= '*.seq'):
	seqtxt = open(f).read()
	cl = pa.classify(seqtxt)
	if cl in ['blast','ssearch blast']:
		o = pa.fil_blast(seqtxt)
		open(f,'w').write(o)
		i += 1
	elif cl == 'maf':
		o = pa.fil_maf(seqtxt)
		open(f,'w').write(o)
		i += 1

print(i)
