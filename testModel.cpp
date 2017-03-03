/**
 * @file
 * @author Leo Cazenille <leo.cazenille@gmail.com>
 *
 *
 */

#include "model.hpp"
#include "bm.hpp"
#include "factory.hpp"
#include "random.h"

#include <stdlib.h>
#include <unistd.h>
#include <fstream>
#include <iostream>

#include <boost/algorithm/string.hpp>

using namespace Fishmodel;
using namespace std;

// For getopt
extern char *optarg;
extern int  optind, opterr;

void usage() {
	std::cout << "model:" << std::endl;
	std::cout << "\t-h" << std::endl;
	std::cout << "\t-n nbSteps" << std::endl;
	std::cout << "\t-F nbFishes" << std::endl;
	std::cout << "\t-i arenaImageFilename" << std::endl;
	std::cout << "\t-s seed" << std::endl;
	std::cout << "\t-d dt" << std::endl;
	std::cout << "\t-I inputPositions" << std::endl;
	std::cout << std::endl;
}


int main(int argc, char** argv) {
	size_t nbSteps = 1000;
	size_t nbFishes = 10;
	size_t nbZones = 1;
	std::string imgFilename = "arena.png";
	std::string inputFilename = "";
	unsigned seed = 0;
	double dt = 1. / 3.;

	int c;
	static char optstring[] = "n:F:i:s:d:I:z:";
	opterr=0;
	while ((c=getopt(argc, argv, optstring)) != -1) {
		switch(c) {
			case 'n':
				nbSteps = atoi(optarg);
				break;

			case 'F':
				nbFishes = atoi(optarg);
				break;

			case 'i':
				imgFilename = optarg;
				break;

			case 's':
				seed = atoi(optarg);
				break;

			case 'd':
				dt = 1. / static_cast<double>(atoi(optarg));
				break;

			case 'I':
				inputFilename = optarg;
				break;

			case 'z':
				nbZones = atoi(optarg);
				break;

			case 'h':
			default:
				usage();
				return 0;
				break;
		}
	}

	// Init Random Number Generator
	rne().seed(seed);

	// Launch model simulation
	Arena arena(imgFilename);
	//FishSimulation sim(arena, nbFishes);

	//for(size_t i = 0; i < nbSteps; ++i) {
	//	double t = static_cast<double>(i) / 15.0;
	//	std::cout << t << " ";
	//	std::cout << sim.printCurrentPositions() << std::endl;
	//	sim.step();
	//}

	SimulationFactory factory(arena);
//	factory.nbFishes = nbFishes - 3;
	factory.nbFishes = nbFishes;
	factory.nbRobots = 0;
	factory.nbVirtuals = 0;
	factory.nbZones = nbZones;

	//factory.behaviorFishes = "BM"; //"TrajectoryFollowing";
	//factory.behaviorRobots = "BM";
	//factory.behaviorVirtuals = "BM";
	factory.behaviorFishes = "ZoneDependantBM";
	factory.behaviorRobots = "ZoneDependantBM";
	factory.behaviorVirtuals = "ZoneDependantBM";

	// Find input trajectories, if wanted
	if(inputFilename != "") {
		ifstream in(inputFilename);
		if(in.is_open()) {
			std::string line;
			std::getline(in, line);
			while(std::getline(in, line)) {
				std::vector<std::string> strs;
				boost::split(strs, line, boost::is_any_of(" \t"));
				strs.erase(strs.begin());
				if(factory.trajectories.size() == 0) {
					factory.trajectories.assign(strs.size() / 2, std::vector<Coord3D_t>());
				}
				for(size_t i = 0; i < strs.size() / 2; ++i) {
					//factory.trajectories[i].push_back({ std::stof(strs[i*3]), std::stof(strs[i*3+1]), std::stof(strs[i*3+2]) });
					factory.trajectories[i].push_back({ std::stof(strs[i*2]), std::stof(strs[i*2+1]), 0.0 });
				}
			}
		}
	}

	auto sim = factory.create();
	sim->dt = dt;

	std::cout << std::fixed; // << std::setprecision(3);
	for(size_t i = 0; i < nbSteps; ++i) {
		double t = static_cast<double>(i) / 15.0;
		std::cout << t << "\t";
		std::cout << sim->printCurrentPositions() << std::endl;
		//std::cout << std::fixed << sim->printCurrentSpeeds() << std::endl;
		sim->step();
	}

	return 0;
}


// MODELINE	"{{{1
// vim:noexpandtab:softtabstop=4:shiftwidth=4:fileencoding=utf-8
// vim:foldmethod=marker
