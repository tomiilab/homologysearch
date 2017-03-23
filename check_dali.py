import seqlim
from path import Path

'''
compare dali alignments and reference seqs.
to check missing residues in dali alignments
'''

um = []
all = []
files = Path('./miqs_yamada_benchmark_dataset/dali_reference/').walkfiles(pattern='*.dal')
for f in files:
	tags = f.split('/')[-1].split('.')[:2]
	true_seqs = [e.split()[1].upper().replace('-','') for e in open(f).read().splitlines()[1:]]
	seq1 = seqlim.MSeq.parse(open('./repo/cath20-scop/'+tags[0]+'.seq'))[0].seq.upper()
	seq2 = seqlim.MSeq.parse(open('./repo/cath20-scop/'+tags[1]+'.seq'))[0].seq.upper()
	
	true_seqs[0] = true_seqs[0].replace('X','')
	true_seqs[1] = true_seqs[1].replace('X','')

	if seq1 != true_seqs[0]:
		print tags[0]
		print seq1
		print true_seqs[0]
		um.append(tags[0])

	if seq2 != true_seqs[1]:
		print tags[1]
		print seq2
		print true_seqs[1]
		um.append(tags[1])
	all.extend(tags)

print len(list(set(um)))
print len(list(set(all)))
