## Convert blastp, psi-blast, last, ssearch and cs-blast results to tsv files.
`python filter_aln.py ./repo/raw/YOURDIR`

## Oda et al. Fig 2
`python draw_blocksize_density.py`

## Oda et al. Fig 3A
`python calc_roc.py draw_roc ./args/scop20_training.iter5.blocksizes.json -standards JG -ylim 1100`

## Oda et al. Fig 3B
`python calc_roc.py draw_roc ./args/scop20_validation.iter2_3_5.json -standards JG -ylim 1100`

## Oda et al. Fig 3C
`python calc_roc.py draw_roc ./args/cath20-scop.iter5.json -standards superfamily -ylim 600`

## Oda et al. Fig 5
`python calc_roc.py draw_roc ./args/scop20_validation.msa.iter3_5.json -standards JG -ylim 1100`

## Oda et al. Fig 4A raw data
`python calc_roc.py roc5_stats ./args/scop20_training.iter2.json -standards JG -rocn_max_evalue 1.0`

## Oda et al. Fig 4B raw data
`python calc_roc.py roc5_stats ./args/scop20_validation.iter2.json -standards JG -rocn_max_evalue 1.0`

## Oda et al. Fig 4C raw data
`python calc_roc.py roc5_stats ./args/cath20-scop.iter2.json -standards superfamily -rocn_max_evalue 1.0`

