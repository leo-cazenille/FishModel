/**
 * @file
 * @author Leo Cazenille <leo.cazenille@gmail.com>
 *
 *
 */

#include <boost/math/special_functions/bessel.hpp>
#include <iostream>
#include <limits>

#include "model.hpp"
#include "random.h"

using namespace Fishmodel;
using namespace std;
using namespace CATS;

/////////////////////////////////////////////////////////
/////////////////// Utility functions /////////////////// {{{1
/////////////////////////////////////////////////////////
namespace {

template<typename val_t = real_t>
inline val_t constexpr normAngle(val_t const angle) {
	return ::fmod(angle + 2 * pi(), 2 * pi());
}

}



/////////////////////////////////////////////////////////
///////////////////////// Arena ///////////////////////// {{{1
/////////////////////////////////////////////////////////

void Arena::resetDistanceMatrix() {
	//distanceMatrix = Mat_t(distanceMatrixDim, distanceMatrixDim, CV_8S, cv::Scalar::all(-1));
	distanceMatrix = Mat_t(distanceMatrixDim, distanceMatrixDim, CV_8S);
	nearestObstacleCoordsMatrix = std::vector<std::vector<IndexCoord_t>>(distanceMatrixDim, std::vector<IndexCoord_t>(distanceMatrixDim, {0, 0}));
	tangentMatrix = Mat_t(distanceMatrixDim, distanceMatrixDim, CV_32F);

	// Identify location of walls
	size_t remainingCells = distanceMatrixDim * distanceMatrixDim;
	for(size_t i = 0; i < distanceMatrixDim; ++i) {
		for(size_t j = 0; j < distanceMatrixDim; ++j) {
			Coord_t const c = {
				static_cast<real_t>(i) / static_cast<real_t>(distanceMatrixDim) * size.first,
				static_cast<real_t>(j) / static_cast<real_t>(distanceMatrixDim) * size.second,
			};
			if(!isInArena(c)) {
				distanceMatrix.at<char>(i, j) = 0;
				--remainingCells;
			} else {
				distanceMatrix.at<char>(i, j) = -1;
			}
		}
	}

	// Propagate distance
	bool madeChanges;
	do {
		madeChanges = false;
		for(size_t i = 0; i < distanceMatrixDim; ++i) {
			for(size_t j = 0; j < distanceMatrixDim; ++j) {
				if(distanceMatrix.at<char>(i, j) == 0 || distanceMatrix.at<char>(i, j) == 1)
					continue;
				size_t minDist;
				if(distanceMatrix.at<char>(i, j) != -1)
					minDist = distanceMatrix.at<char>(i, j) - 1;
				else
					minDist = distanceMatrixDim * distanceMatrixDim;
				if(i > 0 && j > 0 && distanceMatrix.at<char>(i - 1, j - 1) != -1)
					if(distanceMatrix.at<char>(i - 1, j - 1) < minDist) minDist = distanceMatrix.at<char>(i - 1, j - 1);
				if(i > 0 && distanceMatrix.at<char>(i - 1, j) != -1)
					if(distanceMatrix.at<char>(i - 1, j) < minDist) minDist = distanceMatrix.at<char>(i - 1, j);
				if(i > 0 && j < distanceMatrixDim - 1 && distanceMatrix.at<char>(i - 1, j + 1) != -1)
					if(distanceMatrix.at<char>(i - 1, j + 1) < minDist) minDist = distanceMatrix.at<char>(i - 1, j + 1);
				if(j > 0 && distanceMatrix.at<char>(i, j - 1) != -1)
					if(distanceMatrix.at<char>(i, j - 1) < minDist) minDist = distanceMatrix.at<char>(i, j - 1);
				if(j < distanceMatrixDim - 1 && distanceMatrix.at<char>(i, j + 1) != -1)
					if(distanceMatrix.at<char>(i, j + 1) < minDist) minDist = distanceMatrix.at<char>(i, j + 1);
				if(i < distanceMatrixDim - 1 && j > 0 && distanceMatrix.at<char>(i + 1, j - 1) != -1)
					if(distanceMatrix.at<char>(i + 1, j - 1) < minDist) minDist = distanceMatrix.at<char>(i + 1, j - 1);
				if(i < distanceMatrixDim - 1 && distanceMatrix.at<char>(i + 1, j) != -1)
					if(distanceMatrix.at<char>(i + 1, j) < minDist) minDist = distanceMatrix.at<char>(i + 1, j);
				if(i < distanceMatrixDim - 1 && j < distanceMatrixDim - 1 && distanceMatrix.at<char>(i + 1, j + 1) != -1)
					if(distanceMatrix.at<char>(i + 1, j + 1) < minDist) minDist = distanceMatrix.at<char>(i + 1, j + 1);
				++minDist;
				if(minDist < distanceMatrix.at<char>(i, j)) {
					distanceMatrix.at<char>(i, j) = minDist;
					madeChanges = true;
				}
			}
		}
	} while(madeChanges);

	//std::cout << "DEBUG2: " << distanceMatrix << std::endl;
}



/////////////////////////////////////////////////////////
///////////////////////// Agent ///////////////////////// {{{1
/////////////////////////////////////////////////////////

void Agent::reinit() {
	lognormal_distribution<real_t> speedDistrib(meanSpeed, varSpeed);
	direction = unifAngle(rne());
	speed = speedDistrib(rne());

	angleFishSide = acos((length/3.0) / sqrt((length/3.0) * (length/3.0) + (width/2.0) * (width/2.0)));
	normBody = sqrt((length / 3.0) * (length / 3.0) + (width / 2.0) * (width / 2.0));

	// Initialize agent position
	uniform_real_distribution<real_t> unifSizeX(0.0, arena.size.first);
	uniform_real_distribution<real_t> unifSizeY(0.0, arena.size.second);
	for(;;) {
		Coord_t pos(unifSizeX(rne()), unifSizeY(rne()));
		Coord_t newTailPos(
				pos.first - cos(direction) * length,
				pos.second - sin(direction) * length);
		if(arena.isInArena(pos) && arena.isInArena(newTailPos)) {
			headPos = pos;
			tailPos = newTailPos;
			break;
		}
	}
}

void Agent::updateAgentPosition(real_t dt) {
	Coord_t newHeadPos = {
		headPos.first + dt * speed * cos(direction),
		headPos.second + dt * speed * sin(direction)
	};

	if(newHeadPos.first > arena.size.first) {
		newHeadPos.first = arena.size.first - (newHeadPos.first - arena.size.first);
		direction = pi() - direction;
	} else if(newHeadPos.first < 0.0) {
		newHeadPos.first = -newHeadPos.first;
		direction = pi() - direction;
	}

	if(newHeadPos.second > arena.size.second) {
		newHeadPos.second = arena.size.second - (newHeadPos.second - arena.size.second);
		direction = - direction;
	} else if(newHeadPos.second < 0.0) {
		newHeadPos.second = -newHeadPos.second;
		direction = - direction;
	}

	if(!arena.isInArena(newHeadPos)) {
		direction = - direction;
		newHeadPos = {
			headPos.first + dt * speed * cos(direction),
			headPos.second + dt * speed * sin(direction)
		};
	}

	if(arena.isInArena(newHeadPos)) {
		headPos = newHeadPos;
	}

	tailPos = {
		headPos.first - cos(direction) * length,
		headPos.second - sin(direction) * length
	};
}

std::string Agent::printCurrentPositions() const {
	stringstream ss;
	ss << headPos.first << " " << headPos.second;
	ss << " " << direction;
	return ss.str();
}

std::string Agent::printCurrentSpeeds() const {
	stringstream ss;
	ss << std::fixed << speed;
	return ss.str();
}



/////////////////////////////////////////////////////////
///////////////// StraightAheadBehavior ///////////////// {{{1
/////////////////////////////////////////////////////////

void StraightAheadBehavior::step() {
	_agent->updateAgentPosition(_simulation.dt);
}


/////////////////////////////////////////////////////////
///////////////////// RandomBehavior //////////////////// {{{1
/////////////////////////////////////////////////////////

void RandomBehavior::reinit() {
	// TODO
}

void RandomBehavior::step() {
	// TODO
	_agent->updateAgentPosition(_simulation.dt);
}

/////////////////////////////////////////////////////////
/////////////////// RandomWalkBehavior ////////////////// {{{1
/////////////////////////////////////////////////////////

void RandomWalkBehavior::reinit() {
	// TODO
}

void RandomWalkBehavior::step() {
	// TODO
	_agent->updateAgentPosition(_simulation.dt);
}

/////////////////////////////////////////////////////////
///////////// TrajectoryFollowingBehaviour ////////////// {{{1
/////////////////////////////////////////////////////////

void TrajectoryFollowingBehaviour::reinit() {
	_currentPointIndex = 0;
	step();
}

// TODO
void TrajectoryFollowingBehaviour::step() {
	if(_currentPointIndex < _trajectory.size()) {
		_agent->headPos.first = _trajectory[_currentPointIndex][0];
		_agent->headPos.second = _trajectory[_currentPointIndex][1];
		_agent->direction = _trajectory[_currentPointIndex][2];
	} else {
		_agent->headPos.first = 0.;
		_agent->headPos.second = 0.;
		_agent->direction = 0.;
	}
	_agent->tailPos = {
		_agent->headPos.first - cos(_agent->direction) * _agent->length,
		_agent->headPos.second - sin(_agent->direction) * _agent->length
	};
	++_currentPointIndex;
}



/////////////////////////////////////////////////////////
/////////////////// CouzinBehavior /////////////////////// {{{1
/////////////////////////////////////////////////////////

void CouzinBehavior::reinit() {
	//_directionNoiseDistribution = std::uniform_real_distribution<double>(-directionNoiseStdDev/2.0, directionNoiseStdDev/2.0);
	_directionNoiseDistribution = std::normal_distribution<double>(0., 0.1);
}

real_t CouzinBehavior::_computeAgentSpeed() {
	lognormal_distribution<real_t> speedDistrib(log(0.12), 0.045); ///(_agent->meanSpeed, _agent->varSpeed);
	return speedDistrib(rne());
	//return _agent->meanSpeed; // XXX
	//return 0.06; // XXX
}

real_t CouzinBehavior::_computeAgentDeltaDirection() {
	real_t deltaDirection = 0.0;
	real_t normalization = 0.0;
	real_t deltaDirection1x = 0.0;
	real_t deltaDirection1y = 0.0;
	real_t deltaDirection2x = 0.0;
	real_t deltaDirection2y = 0.0;
	size_t nbRepulsed = 0;

	for(auto& _a: _simulation.agents) {
		Agent* a = _a.first.get();
		if(a == _agent)
			continue;

		real_t const relativeAgentsHeadsPosX = a->headPos.first - _agent->headPos.first;
		real_t const relativeAgentsHeadsPosY = a->headPos.second - _agent->headPos.second;
		real_t const relativeAgentsTailPosX = a->tailPos.first - _agent->tailPos.first;
		real_t const relativeAgentsTailPosY = a->tailPos.second - _agent->tailPos.second;
		real_t const headRetinaX = relativeAgentsHeadsPosX * cos(_agent->direction) + relativeAgentsHeadsPosY * sin(_agent->direction);
		real_t const headRetinaY = -relativeAgentsHeadsPosX * sin(_agent->direction) + relativeAgentsHeadsPosY * cos(_agent->direction);
		real_t const tailRetinaX = relativeAgentsTailPosX * cos(_agent->direction) + relativeAgentsTailPosY * sin(_agent->direction);
		real_t const tailRetinaY = -relativeAgentsTailPosX * sin(_agent->direction) + relativeAgentsTailPosY * cos(_agent->direction);
		real_t const headAngle = normAngle(atan2(headRetinaY, headRetinaX));
		real_t const distAgent = sqrt(headRetinaX * headRetinaX + headRetinaY * headRetinaY);
		real_t const tailAngle = normAngle(atan2(tailRetinaY, tailRetinaX));
		//real_t const agentCenter = normAngle(atan2(cos(headAngle) + cos(tailAngle), sin(headAngle) + sin(tailAngle)));
		real_t const agentCenter = headAngle;
		real_t const diffDirections = normAngle(a->direction - _agent->direction);

		// Find if the agent is in the attraction zone
		if(		(distAgent < attractionDistance) &&
				((agentCenter < _agent->fov * pi() / 180.0) || (agentCenter > (360. - _agent->fov) * pi() / 180.0))) {
			normalization += 1.0;
			if(distAgent < repulsionDistance) {
				++nbRepulsed;
				//deltaDirection1 += normAngle(agentCenter + pi());
				deltaDirection1x += cos(normAngle(agentCenter + pi()));
				deltaDirection1y += sin(normAngle(agentCenter + pi()));
			} else {
				//deltaDirection2 += diffDirections / 2.0;
				deltaDirection2x += cos(diffDirections);
				deltaDirection2y += sin(diffDirections);
			}
		}
	}

	if(nbRepulsed > 0) {
		//deltaDirection = deltaDirection1;
		deltaDirection = atan2(deltaDirection1y, deltaDirection1x);
	} else {
		//deltaDirection = deltaDirection2;
		deltaDirection = atan2(deltaDirection2y, deltaDirection2x);
	}

	if(deltaDirection > pi())
		deltaDirection -= 2. * pi();

	//if(normalization > 0.0) {
	//	deltaDirection /= normalization;
	//}

	return deltaDirection;
}

void CouzinBehavior::step() {
	//for(;;) {
		real_t deltaDirection = _computeAgentDeltaDirection();
		//std::cout << "DEBUG deltaDirection=" << deltaDirection << std::endl;

		//deltaDirection = normAngle(deltaDirection);
		if(deltaDirection > maxRotationAngle * _simulation.dt)
			deltaDirection = maxRotationAngle * _simulation.dt;
		else if(deltaDirection < - maxRotationAngle * _simulation.dt)
			deltaDirection = - maxRotationAngle * _simulation.dt;

		real_t newSpeed = _computeAgentSpeed();
		real_t directionNoise = _directionNoiseDistribution(rne());
		real_t newDirection = _agent->direction + deltaDirection + directionNoise;
		newDirection = normAngle(newDirection);

		Coord_t newAgentsHeadPos = {
			_agent->headPos.first + _simulation.dt * newSpeed * cos(newDirection),
			_agent->headPos.second + _simulation.dt * newSpeed * sin(newDirection)
		};

		if(_simulation.arena.isInArena(newAgentsHeadPos)) {
			_agent->direction = newDirection;
			_agent->speed = newSpeed;
			_agent->headPos = newAgentsHeadPos;
			return; //break;
		}

		newDirection = normAngle(newDirection + pi());
		newSpeed = 0.0; // /= 5.0;
		newAgentsHeadPos = {
			_agent->headPos.first + _simulation.dt * newSpeed * cos(newDirection),
			_agent->headPos.second + _simulation.dt * newSpeed * sin(newDirection)
		};
		if(_simulation.arena.isInArena(newAgentsHeadPos)) {
			_agent->direction = newDirection;
			_agent->speed = newSpeed;
			_agent->headPos = newAgentsHeadPos;
		}
	//}
}



/////////////////////////////////////////////////////////
////////////////////// Simulation /////////////////////// {{{1
/////////////////////////////////////////////////////////

void Simulation::step() {
#ifdef ENABLE_OPENMP
	#pragma omp parallel for
#endif
	for(size_t i = 0; i < agents.size(); ++i) {
		agents[i].second->step();
	}
	++currentTimeStep;
}

std::string Simulation::printCurrentPositions() const {
	stringstream ss;
	for(size_t i = 0; i < agents.size(); ++i) {
		ss << agents[i].first->printCurrentPositions();
		if(i != agents.size() - 1)
			ss << "\t";
	}
	return ss.str();
}

std::string Simulation::printCurrentSpeeds() const {
	stringstream ss;
	for(size_t i = 0; i < agents.size(); ++i) {
		ss << std::fixed << agents[i].first->printCurrentSpeeds();
		if(i != agents.size() - 1)
			ss << "\t";
	}
	return ss.str();
}




// MODELINE	"{{{1
// vim:noexpandtab:softtabstop=4:shiftwidth=4:fileencoding=utf-8
// vim:foldmethod=marker
