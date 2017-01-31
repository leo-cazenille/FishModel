/**
 * @file
 * @author Leo Cazenille <leo.cazenille@gmail.com>
 *
 *
 */

#include <iostream>

#include "zones.hpp"

using namespace Fishmodel;
using namespace std;


/////////////////////////////////////////////////////////
/////////////// ZoneDependantBehavior /////////////////// {{{1
/////////////////////////////////////////////////////////

void ZoneDependantBehavior::reinit() {
	for(auto& z: _zones) {
		z->reinit();
	}
}

void ZoneDependantBehavior::step() {
	auto& currentZone = _zones[this->currentZoneId()];
	//std::cout << "DEBUG1: " << this->currentZoneId() << std::endl;
	currentZone->step();
}


// MODELINE	"{{{1
// vim:noexpandtab:softtabstop=4:shiftwidth=4:fileencoding=utf-8
// vim:foldmethod=marker
