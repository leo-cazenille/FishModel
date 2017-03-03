/**
 * @file
 * @author Leo Cazenille <leo.cazenille@gmail.com>
 *
 *
 */

#include <boost/math/special_functions/bessel.hpp>
#include <iostream>
#include <limits>

#include "factory.hpp"
#include "random.h"

using namespace Fishmodel;
using namespace std;


/////////////////////////////////////////////////////////
////////////////// SimulationFactory //////////////////// {{{1
/////////////////////////////////////////////////////////

template <class A>
std::pair<Agent*, Behavior*> SimulationFactory::_createAgent(std::string const& behaviorType) {
	Agent* a = new A(arena);
	Behavior* b = nullptr;
	if(behaviorType == "NoBehavior") {
		b = new NoBehavior(*_sim, a);
	} else if(behaviorType == "StraightAhead") {
		b = new StraightAheadBehavior(*_sim, a);
	} else if(behaviorType == "Random") {
		b = new RandomBehavior(*_sim, a);
	} else if(behaviorType == "RandomWalk") {
		b = new RandomWalkBehavior(*_sim, a);
	} else if(behaviorType == "Couzin") {
		b = new CouzinBehavior(*_sim, a);
	//} else if(behaviorType == "Model") {
	//	b = new ModelBehavior(*_sim, a);
	//} else if(behaviorType == "Model3D") {
	//	b = new Model3DBehavior(*_sim, a);
	} else if(behaviorType == "BM") {
		b = new BM(*_sim, a);
	} else if(behaviorType == "TrajectoryFollowing") {
		b = new TrajectoryFollowingBehaviour(*_sim, a);
	} else if(behaviorType == "ZoneDependantBM") {
		std::vector<std::shared_ptr<Behavior>> zones;
		zones.emplace_back(new NoBehavior(*_sim, a));
		for(size_t i = 0; i < nbZones; ++i) {
			zones.emplace_back(new ZonedBM(*_sim, a));
		}
		auto* zdb = new ZoneDependantBehavior(*_sim, a, _sim->arena, zones);
		auto affinities = std::vector<real_t>(nbZones, 1.0);
		for(size_t i = 0; i < nbZones; ++i) {
			reinterpret_cast<ZonedBM*>(zones[i+1].get())->zdb(zdb);
			reinterpret_cast<ZonedBM*>(zones[i+1].get())->zonesAffinity(affinities);
		}
		b = zdb;
	} else {
		b = new Behavior(*_sim, a);
	}
	return {a, b};
}

std::unique_ptr<Simulation> SimulationFactory::create() {
	_sim = new Simulation(arena);
	std::unique_ptr<Simulation> sim(_sim);
	currentTrajectoryIndex = 0;

	// Create Agents
	for(size_t i = 0; i < nbFishes; ++i) {
		auto ab = _createAgent<FishAgent>(behaviorFishes);
		sim->agents.push_back({std::unique_ptr<Agent>(ab.first), std::unique_ptr<Behavior>(ab.second)});
		sim->fishes.push_back(ab);
	}
	for(size_t i = 0; i < nbRobots; ++i) {
		auto ab = _createAgent<VirtualAgent>(behaviorRobots);
		sim->agents.push_back({std::unique_ptr<Agent>(ab.first), std::unique_ptr<Behavior>(ab.second)});
		sim->robots.push_back(ab);
	}
	for(size_t i = 0; i < nbVirtuals; ++i) {
		auto ab = _createAgent<VirtualAgent>(behaviorVirtuals);
		sim->agents.push_back({std::unique_ptr<Agent>(ab.first), std::unique_ptr<Behavior>(ab.second)});
		sim->virtuals.push_back(ab);
	}
	return std::move(sim);
}



// MODELINE	"{{{1
// vim:noexpandtab:softtabstop=4:shiftwidth=4:fileencoding=utf-8
// vim:foldmethod=marker
