/**
 * @file
 * @author Leo Cazenille <leo.cazenille@gmail.com>
 *
 *
 */

#ifndef MODEL_H
#define MODEL_H

#include <utility>
#include <vector>
#include <string>
#include <memory>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>

#include <iostream>


namespace Fishmodel {

typedef double real_t;

typedef std::pair<real_t, real_t> Coord_t;
typedef std::array<real_t, 3> Coord3D_t;
typedef cv::Mat Mat_t;
typedef std::pair<size_t, size_t> IndexCoord_t;


std::size_t constexpr evalAnglesNb = 360;


template<typename val_t = real_t>
inline val_t constexpr pi() { return std::atan(1.)*4.; }


/** Arena used in simulation, with matrix-based support of obstacles */
class Arena {
public:
	Mat_t arenaMatrix;
	Mat_t distanceMatrix;
	std::vector<std::vector<IndexCoord_t>> nearestObstacleCoordsMatrix;
	Mat_t tangentMatrix;
	Coord_t size = {1.0, 1.0};
	Coord_t pixelSize;
	size_t distanceMatrixDim = 1000;

public:
	Arena(Mat_t& _arenaMatrix, Coord_t _size = {1.0, 1.0}, size_t _distanceMatrixDim = 100) : arenaMatrix(_arenaMatrix), size(_size), distanceMatrixDim(_distanceMatrixDim) {
		pixelSize = discreteToRealCoords(this->size);
		resetDistanceMatrix();
	}
	Arena(std::string filename, Coord_t _size = {1.0, 1.0}, size_t _distanceMatrixDim = 100) : size(_size), distanceMatrixDim(_distanceMatrixDim) {
		arenaMatrix = cv::imread(filename, CV_LOAD_IMAGE_GRAYSCALE); // XXX verify
		pixelSize = discreteToRealCoords(this->size);
		resetDistanceMatrix();
	}

	inline Coord_t discreteToRealCoords(Coord_t const c) const {
		return Coord_t(
			c.first * this->size.first / this->arenaMatrix.cols,
			c.second * this->size.second / this->arenaMatrix.rows
		);
	}

	inline Coord_t realToDiscreteCoords(Coord_t const c) const {
		return Coord_t(
			c.first / this->size.first * this->arenaMatrix.cols,
			c.second / this->size.second * this->arenaMatrix.rows
		);
	}

	inline bool isInArenaDiscrete(Coord_t const c) const {
		if(c.first < 0.0 || c.second < 0.0 || c.first >= this->arenaMatrix.cols || c.second >= this->arenaMatrix.rows) {
			return false;
		} else {
			//return this->arenaMatrix.at<uchar>(this->arenaMatrix.rows - c.first, this->arenaMatrix.cols - c.second) != 0;
			return this->arenaMatrix.at<uchar>(c.second, c.first) != 0;
		}
	}

	inline bool isInArena(Coord_t const c) const {
		return this->isInArenaDiscrete(this->realToDiscreteCoords(c));
	}

	//template <typename ret_t=uchar>
	//inline ret_t at(Coord_t const c) const {
	inline uchar at(Coord_t const c) const {
		Coord_t const discreteCoords = this->realToDiscreteCoords(c);
		return this->arenaMatrix.at<uchar>(discreteCoords.second, discreteCoords.first);
	}

	inline double getDistanceToWall(Coord_t const c) const {
		return (size.first / distanceMatrix.rows) * distanceMatrix.at<uchar>(
			static_cast<size_t>(c.first / size.first * distanceMatrix.rows),
			static_cast<size_t>(c.second / size.second * distanceMatrix.cols));
	}

	void resetDistanceMatrix();
};


struct Agent {
	bool present = true;

	real_t length = 0.02; //0.035;
	real_t width = 0.02; //0.01;
	real_t height = 0.01;
	std::size_t fov = 135;

	real_t meanSpeed = -2.65;
	real_t varSpeed = 0.51;

	//real_t maxTurningRate = pi() / 2.;
	real_t maxTurningRate = 2. * pi();

	real_t wallDistancePerception = 0.05;
	real_t agentDistancePerception = 3.0;
	real_t agentMinAnglePerception = 0.0;

	Arena& arena;

	real_t direction;
	real_t speed;
	Coord_t headPos;
	Coord_t tailPos;

public:
	real_t angleFishSide;
	real_t normBody;

public:
	Agent(Arena& _arena) : arena(_arena) {
		reinit();
	}
	virtual ~Agent() {}

	virtual void reinit();
	virtual void updateAgentPosition(real_t dt);
	virtual std::string printCurrentPositions() const;
	virtual std::string printCurrentSpeeds() const;
};


struct FishAgent : public Agent {
	FishAgent(Arena& _arena) : Agent(_arena) {}
};

struct VirtualAgent : public Agent {
	VirtualAgent(Arena& _arena) : Agent(_arena) {}
};


struct Simulation;

class Behavior {
protected:
	Simulation& _simulation;
	Agent* _agent;

public:
	Behavior(Simulation& simulation, Agent* agent = nullptr) : _simulation(simulation), _agent(agent) {
		reinit();
	}
	virtual ~Behavior() {}

	virtual void reinit() {}
	virtual void step() {}
};

struct NoBehavior : public Behavior {
	virtual void reinit() {}
	virtual void step() {}
	NoBehavior(Simulation& simulation, Agent* agent = nullptr) : Behavior(simulation, agent) {}
};

struct StraightAheadBehavior : public Behavior {
	virtual void step();
	StraightAheadBehavior(Simulation& simulation, Agent* agent = nullptr) : Behavior(simulation, agent) {}
};

struct RandomBehavior : public Behavior {
	virtual void reinit();
	virtual void step();
	RandomBehavior(Simulation& simulation, Agent* agent = nullptr) : Behavior(simulation, agent) {}
};

struct RandomWalkBehavior : public Behavior {
	virtual void reinit();
	virtual void step();
	RandomWalkBehavior(Simulation& simulation, Agent* agent = nullptr) : Behavior(simulation, agent) {}
};


class TrajectoryFollowingBehaviour : public Behavior {
protected:
	std::vector<Coord3D_t> _trajectory;
	size_t _currentPointIndex = 0;

public:
	virtual void reinit();
	virtual void step();
	TrajectoryFollowingBehaviour(Simulation& simulation, Agent* agent = nullptr) : Behavior(simulation, agent) { reinit(); }
	void trajectory(decltype(_trajectory)& trajectory) { _trajectory = trajectory; reinit(); }
};


struct CouzinBehavior : public Behavior {
public:
	real_t repulsionDistance = 0.02;
	//real_t alignmentDistance = 0.1;
	real_t attractionDistance = 1.50;
	real_t directionNoiseStdDev = 0.0;
	real_t maxRotationAngle = 20.0;

protected:
	//std::uniform_real_distribution<double> _directionNoiseDistribution;
	std::normal_distribution<double> _directionNoiseDistribution;

	virtual real_t _computeAgentSpeed();
	virtual real_t _computeAgentDeltaDirection();

public:
	virtual void reinit();
	virtual void step();

	CouzinBehavior(Simulation& simulation, Agent* agent = nullptr) :
		Behavior(simulation, agent) {
		reinit();
	}
};




typedef std::pair<std::unique_ptr<Agent>, std::unique_ptr<Behavior>> AgentBehaviorStorage_t;
typedef std::pair<Agent*, Behavior*> AgentBehavior_t;

struct Simulation {
public:
	std::vector<AgentBehaviorStorage_t> agents;
	Arena& arena;

public:
	std::vector<AgentBehavior_t> fishes;
	std::vector<AgentBehavior_t> robots;
	std::vector<AgentBehavior_t> virtuals;

public:
	std::size_t currentTimeStep = 0;
	real_t dt = 1.0 / 3.0;

	Simulation(Arena& _arena) : arena(_arena) {}
	virtual ~Simulation() {}

	void step();

	std::string printCurrentPositions() const;
	std::string printCurrentSpeeds() const;
};


}
#endif

// MODELINE	"{{{1
// vim:noexpandtab:softtabstop=4:shiftwidth=4:fileencoding=utf-8
// vim:foldmethod=marker
