from matplotlib import use
use('Agg')
from matplotlib import pyplot as plt
from matplotlib import rcParams, rc
import numpy as np
from pyplot_ext import *

fontsize = 7
rcParams.update({'font.size': fontsize})

algos = ['LAST5','LAST6','LAST5_small','LAST6_small','BLASTP','CS-BLAST','SSEARCH']


def parse(text):
	M = np.zeros((len(algos),2))
	for l in text.splitlines():
		cs = l.split()
		if len(cs) == 1:
			continue
		if len(cs) == 3:
			algo, time1, ram1 = cs
			print algo, cs
			if algo in algos:
				M[algos.index(algo)]= [float(time1), float(ram1)]
	return M


def draw(m, outpath):
	print m
	coor = Coor((2.7, 3.2),(.5,.1),(.6,.2),(.1,.1),(1,1))
	print coor.figsize()
	fig = plt.figure(**coor.figsize())
	fig.subplots_adjust(**coor.adjust())

	ax1 = fig.add_subplot(111)

	width = 0.35
	obs = []
	names = []
	for i,e in enumerate(m):
		ob = ax1.plot([e[1]], [e[0]], markersize=8, linewidth=0, markeredgecolor='none',marker='o', color='deepskyblue')
		ob = ax1.text(e[1]+1, e[0], algos[i])


	ax1.set_ylabel('computation time (seconds)')
	ax1.set_xlim([0, 27])
	ax1.set_ylim([10, 5000])
	ax1.set_yscale('log')
	ax1.set_xlabel('max memory usage (GB)')
	fig.savefig(outpath,dpi=600)


if __name__ == '__main__':
	M = parse(open('./archive/speed.txt').read())
	draw(M, './fig/speed.png')

