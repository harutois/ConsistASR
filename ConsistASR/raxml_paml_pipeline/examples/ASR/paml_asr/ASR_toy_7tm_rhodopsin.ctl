
seqfile = ../toy_7tm_rhodopsin.fasta
treefile = ../toy_7tm_rhodopsin.raxml.bestTree
outfile = ASR_toy_7tm_rhodopsin.out

noisy = 9
verbose = 1
runmode = 0

seqtype = 2     * 2: amino acid
CodonFreq = 2

model = 2       * Empirical (e.g. LG)
aaRatefile = ../lg.dat
Mgene = 0

fix_alpha = 0
alpha = 0.5
ncatG = 4

fix_kappa = 1
kappa = 2
fix_omega = 1
omega = 1

RateAncestor = 1   * output ancestral sequences
cleandata = 0
