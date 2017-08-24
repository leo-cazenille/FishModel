#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy.misc import imread


def createPlot():
	fig, ax = plt.subplots()
	points, = ax.plot([], [], marker='o', linestyle='None', color='red')
	pointsTail, = ax.plot([], [], marker='o', linestyle='None', color='green')
	pointsLeft, = ax.plot([], [], marker='o', linestyle='None', color='blue')
	pointsRight, = ax.plot([], [], marker='o', linestyle='None', color='yellow')
	lines = [ax.plot([], [], lw=2)[0] for x in range(agentsNb)]
	return fig, ax, (points, pointsTail, pointsLeft, pointsRight), lines


def plotStep(fig, ax, sim, p, lines):
	(points, pointsTail, pointsLeft, pointsRight) = p

	linesX = zip(sim.agentsHeadsPosX, sim.agentsTailPosX)
	linesY = zip(sim.agentsHeadsPosY, sim.agentsTailPosY)

	ax.set_xlim(0.0, sim.arena.sizeX)
	ax.set_ylim(0.0, sim.arena.sizeY)

	points.set_data(sim.agentsHeadsPosX, sim.agentsHeadsPosY)
	pointsTail.set_data(sim.agentsTailPosX, sim.agentsTailPosY)
	pointsLeft.set_data(sim.agentsLeftSidePosX, sim.agentsLeftSidePosY)
	pointsRight.set_data(sim.agentsRightSidePosX, sim.agentsRightSidePosY)
	#pointsTail.set_data(sim.agentsTailPosX, sim.agentsTailPosY)
	for lindex in range(len(lines)):
		lines[lindex].set_data(linesX[lindex], linesY[lindex])



if __name__ == "__main__":
	from optparse import OptionParser
	from scipy.misc import imread
	usage = "%prog [command] [options]"
	parser = OptionParser(usage=usage)

	parser.add_option("-i", "--inputFilename", dest = "inputFilename", default = "",
			help = "Path of input data file")
	parser.add_option("-n", "--agentsNb", dest = "agentsNb", default = 5,
			help = "Number of agents in simulation")
	parser.add_option("-a", "--arenaImage", dest = "arenaImage", default = "maze.bmp",
			help = "Filename of arena image")
	parser.add_option("-D", "--duration", dest = "duration", default = 60000,
			help = "Duration")

	(options, args) = parser.parse_args()

	arenaMatrix_ = imread(options.arenaImage, flatten=True)
	arenaMatrix = np.empty(shape=arenaMatrix_.shape, dtype=bool)
	arenaMatrix[arenaMatrix_ == 0] = False
	arenaMatrix[arenaMatrix_ != 0] = True

# XXX TODO
#	arena = Arena(arenaMatrix)
	agentsNb = int(options.agentsNb) or 5

	arenaMatrixColor_ = imread(options.arenaImage)
	fig, ax, p, lines = createPlot()
# XXX TODO
#	ax.imshow(arenaMatrixColor_, extent=(0.0, arena.sizeX, arena.sizeY, 0.0), interpolation="none", aspect='auto')


	positions = []
	duration = int(options.duration) or 60000
	for t in range(duration):
		sim.step()
		if options.showAnimation:
			plotStep(fig, ax, sim, p, lines)
			plt.pause(0.001)
		currentPositions = np.hstack(zip(sim.agentsHeadsPosX, sim.agentsHeadsPosY))
		positions.append(currentPositions)
		if t % 1000 == 0:
			print("# ", t, " timestep")
	positions = np.array(positions)




#if __name__ == "__main__":
#	nbAgents = int(sys.argv[2])
#	fig, ax = plt.subplots()
#	ax.set_xlim(0.0, 1.0)
#	ax.set_ylim(0.0, 1.0)
#	colormap = plt.cm.jet
#	ax.set_color_cycle([colormap(i) for i in np.linspace(0, 0.9, nbAgents)])
#
#	arenaMatrixColor_ = imread(sys.argv[3])
#	ax.imshow(arenaMatrixColor_, extent=(0.0, 1.0, 1.0, 0.0), interpolation="none", aspect='auto')
#
#	points = [ax.plot([], [], marker='o', linestyle='None')[0] for i in range(nbAgents)]
#
#	data = np.loadtxt(sys.argv[1], skiprows=1)
#	foo = np.hstack(data[:, 0::3])
#	avgDetected = len(foo[~np.isnan(foo)]) / float(data.shape[0])
#	nbFrame = 0.0
#	textNbFrame = plt.text(0.1, 0.9, "#" + str(nbFrame),
#			horizontalalignment='center',
#			verticalalignment='center',
#			transform = ax.transAxes)
#	print("AvgDetected=" + str(avgDetected))
#	for t in data:
#		for i in range(nbAgents):
#			points[i].set_data(t[1 + i * 3], t[1 + i * 3 + 1])
#		nbFrame = t[0]
#		textNbFrame.set_text("#" + "%.3f" % nbFrame)
#		plt.pause(0.001)

