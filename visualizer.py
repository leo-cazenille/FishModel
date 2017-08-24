#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy.misc import imread



if __name__ == "__main__":
	nbAgents = int(sys.argv[2])
	fig, ax = plt.subplots()
	ax.set_xlim(0.0, 1.0)
	ax.set_ylim(0.0, 1.0)
	colormap = plt.cm.jet
	ax.set_color_cycle([colormap(i) for i in np.linspace(0, 0.9, nbAgents)])

	arenaMatrixColor_ = imread(sys.argv[3])
	ax.imshow(arenaMatrixColor_, extent=(0.0, 1.0, 1.0, 0.0), interpolation="none", aspect='auto')

	points = [ax.plot([], [], marker='o', linestyle='None')[0] for i in range(nbAgents)]

	data = np.loadtxt(sys.argv[1], skiprows=1)
	foo = np.hstack(data[:, 0::3])
	avgDetected = len(foo[~np.isnan(foo)]) / float(data.shape[0])
	nbFrame = 0.0
	textNbFrame = plt.text(0.1, 0.9, "#" + str(nbFrame),
			horizontalalignment='center',
			verticalalignment='center',
			transform = ax.transAxes)
	print("AvgDetected=" + str(avgDetected))
	for t in data:
		for i in range(nbAgents):
			points[i].set_data(t[1 + i * 3], t[1 + i * 3 + 1])
		nbFrame = t[0]
		textNbFrame.set_text("#" + "%.3f" % nbFrame)
		plt.pause(0.001)

