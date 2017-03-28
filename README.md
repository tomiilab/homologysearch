## Convert blastp, psi-blast, last, ssearch and cs-blast results to tsv files.
`python filter_aln.py ./repo/raw/YOURDIR`

## Oda et al. Fig 2
parse block infomation
`python parse_blocksize.py ./repo/raw/scop20_training/psiblast_block_BL62_scop20_training_8/0.002/uniref50/`
`python parse_blocksize.py ./repo/raw/scop20_validation/psiblast_block_BL62_scop20_validation_8/0.002/uniref50/`
`python parse_blocksize.py ./repo/raw/cath20-scop/psiblast_block_BL62_cath20-scop_8/0.002/uniref50/`

draw block size density until the 5th itration
`python draw_blocksize_density.py 5`

draw block size density until the 8th itration
`python draw_blocksize_density.py 8`

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

