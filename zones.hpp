/**
 * @file
 * @author Leo Cazenille <leo.cazenille@gmail.com>
 *
 *
 */

#ifndef ZONES_H
#define ZONES_H

#include <utility>
#include <vector>
#include <string>
#include <memory>

#include <iostream>

#include "model.hpp"


namespace Fishmodel {

//template <typename behavior_t = Behavior>
class ZoneDependantBehavior : public Behavior { // XXX Dependant vs dependent ??
protected:
	Arena* _zonesMatrix;
	std::vector<std::shared_ptr<Behavior>> _zones; // XXX Use a map<size_t, Behavior&> instead ??

public:
	ZoneDependantBehavior(Simulation& simulation, Agent* agent = nullptr) : Behavior(simulation, agent) {}

	ZoneDependantBehavior(Simulation& simulation, Agent* agent, Arena& zonesMatrix, std::vector<std::shared_ptr<Behavior>>& zones) :
			Behavior(simulation, agent),
			_zonesMatrix(&zonesMatrix),
			_zones(zones) {}


	virtual void reinit();
	virtual void step();

	inline size_t currentZoneId() { return static_cast<size_t>(_zonesMatrix->at(_agent->headPos) / 10); } 
	//inline size_t addZone(Behavior& z) { _zones.push_back(std::make_shared<Behavior>(&z)); return _zones.size() - 1; }
	inline Behavior* behavior(size_t index) { return _zones[index].get(); }
	inline Behavior* currentBehavior() { return behavior(currentZoneId()); } 
	inline void zoneMatrix(decltype(_zonesMatrix) zoneMatrix) { _zonesMatrix = zoneMatrix; }

};


}
#endif

// MODELINE	"{{{1
// vim:noexpandtab:softtabstop=4:shiftwidth=4:fileencoding=utf-8
// vim:foldmethod=marker
