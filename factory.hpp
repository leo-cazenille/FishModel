/**
 * @file
 * @author Leo Cazenille <leo.cazenille@gmail.com>
 *
 *
 */

#ifndef FACTORY_H
#define FACTORY_H

#include <vector>
#include <string>

#include "model.hpp"
#include "bm.hpp"
#include "zones.hpp"


namespace Fishmodel {


class SimulationFactory {
public:
	Arena& arena;
	size_t nbFishes = 7;
	size_t nbRobots = 2;
	size_t nbVirtuals = 1;
	std::string behaviorFishes = "Model3D";
	std::string behaviorRobots = "Model3D";
	std::string behaviorVirtuals = "Model3D";

	std::vector<std::vector<Coord3D_t>> trajectories;
	size_t currentTrajectoryIndex = 0;

	size_t nbZones = 1;

protected:
	Simulation* _sim = nullptr;

	template <class A>
	std::pair<Agent*, Behavior*> _createAgent(std::string const& behaviorType);

public:
	SimulationFactory(Arena& _arena) : arena(_arena) {}
	std::unique_ptr<Simulation> create();
};


}
#endif

// MODELINE	"{{{1
// vim:noexpandtab:softtabstop=4:shiftwidth=4:fileencoding=utf-8
// vim:foldmethod=marker
