/**
 * @file
 * @author Leo Cazenille <leo.cazenille@gmail.com>
 *
 *
 */

#include <boost/math/special_functions/bessel.hpp>
#include <iostream>
#include <limits>
#include <exception>

#include "bm.hpp"
#include "random.h"

using namespace Fishmodel;
using namespace std;

/////////////////////////////////////////////////////////
/////////////////// Utility functions /////////////////// {{{1
/////////////////////////////////////////////////////////
namespace {

template<typename val_t = real_t>
inline val_t constexpr normAngle(val_t const angle) {
	return ::fmod(angle + 2 * pi(), 2 * pi());
}

template<typename val_t = real_t>
inline std::pair<val_t, val_t> cart2sph2(val_t const x, val_t const y, val_t const z) {
	val_t const azimuth = ::atan2(y, x);
	val_t const elevation = ::atan2(z, ::sqrt(x * x + y * y));
	return {azimuth, elevation};
}
template<typename val_t = Coord3D_t, typename ret_t = real_t>
inline std::pair<ret_t, ret_t> cart2sph2(val_t const coords) {
	auto const x = coords[0];
	auto const y = coords[1];
	auto const z = coords[2];
	return cart2sph2(x, y, z);
}

template<typename val_t = std::vector<double>, typename ret_t = double>
inline ret_t sum(val_t const l) {
	ret_t s{0.0};
	for(auto foo: l) s += foo;
	return s;
}

template<typename ret_t = real_t>
inline ret_t surfaceSphereFish(std::array<Coord3D_t, 4> const c) {
	std::array<ret_t, 3> d1;
	std::array<ret_t, 3> d2;
	std::array<std::pair<ret_t, ret_t>, std::tuple_size<decltype(c)>::value> sph;
	for(size_t i = 0; i < sph.size(); ++i) sph[i] = cart2sph2(c[i]);
	std::array<size_t, 4> const points0 = {0, 1, 2, 0};
	std::array<size_t, 4> const points1 = {0, 2, 3, 0};

	for(size_t i = 0; i < d1.size(); ++i) {
		d1[i] = acos(sin(sph[points0[i]].second) * sin(sph[points0[i+1]].second) +
					cos(sph[points0[i]].first - sph[points0[i+1]].first) * cos(sph[points0[i]].second) * cos(sph[points0[i+1]].second));
		d2[i] = acos(sin(sph[points1[i]].second) * sin(sph[points1[i+1]].second) +
					cos(sph[points1[i]].first - sph[points1[i+1]].first) * cos(sph[points1[i]].second) * cos(sph[points1[i+1]].second));
	}

	ret_t const s1 = sum(d1) / 2.0;
	ret_t const s2 = sum(d2) / 2.0;
	ret_t const e1 = 4. * atan(sqrt( tan(s1 / 2.) * tan((s1 - d1[0]) / 2.) * tan((s1 - d1[1]) / 2.) * tan((s1 - d1[2]) / 2.0) ));
	ret_t const e2 = 4. * atan(sqrt( tan(s2 / 2.) * tan((s2 - d2[0]) / 2.) * tan((s2 - d2[1]) / 2.) * tan((s2 - d2[2]) / 2.0) ));
	ret_t const result = e1 + e2;
	return result;
}

template<typename val_t, typename vec_t>
inline void normalizeVec(vec_t& v, val_t const coeff) {
	for(size_t k = 0; k < v.size(); ++k) {
		v[k] /= coeff;
	}
}

// NOTE: using constexpr template metaprogramming to populate evaluedProb is bugged with gcc version < 5.0, so we put all values by hand
std::array<real_t, 360> constexpr evaluedProb = {0.0174533, 0.0349066, 0.0523599, 0.0698132, 0.0872665, 0.10472, 0.122173, 0.139626, 0.15708, 0.174533, 0.191986, 0.20944, 0.226893, 0.244346, 0.261799, 0.279253, 0.296706, 0.314159, 0.331613, 0.349066, 0.366519, 0.383972, 0.401426, 0.418879, 0.436332, 0.453786, 0.471239, 0.488692, 0.506145, 0.523599, 0.541052, 0.558505, 0.575959, 0.593412, 0.610865, 0.628319, 0.645772, 0.663225, 0.680678, 0.698132, 0.715585, 0.733038, 0.750492, 0.767945, 0.785398, 0.802851, 0.820305, 0.837758, 0.855211, 0.872665, 0.890118, 0.907571, 0.925025, 0.942478, 0.959931, 0.977384, 0.994838, 1.01229, 1.02974, 1.0472, 1.06465, 1.0821, 1.09956, 1.11701, 1.13446, 1.15192, 1.16937, 1.18682, 1.20428, 1.22173, 1.23918, 1.25664, 1.27409, 1.29154, 1.309, 1.32645, 1.3439, 1.36136, 1.37881, 1.39626, 1.41372, 1.43117, 1.44862, 1.46608, 1.48353, 1.50098, 1.51844, 1.53589, 1.55334, 1.5708, 1.58825, 1.6057, 1.62316, 1.64061, 1.65806, 1.67552, 1.69297, 1.71042, 1.72788, 1.74533, 1.76278, 1.78024, 1.79769, 1.81514, 1.8326, 1.85005, 1.8675, 1.88496, 1.90241, 1.91986, 1.93732, 1.95477, 1.97222, 1.98968, 2.00713, 2.02458, 2.04204, 2.05949, 2.07694, 2.0944, 2.11185, 2.1293, 2.14675, 2.16421, 2.18166, 2.19911, 2.21657, 2.23402, 2.25147, 2.26893, 2.28638, 2.30383, 2.32129, 2.33874, 2.35619, 2.37365, 2.3911, 2.40855, 2.42601, 2.44346, 2.46091, 2.47837, 2.49582, 2.51327, 2.53073, 2.54818, 2.56563, 2.58309, 2.60054, 2.61799, 2.63545, 2.6529, 2.67035, 2.68781, 2.70526, 2.72271, 2.74017, 2.75762, 2.77507, 2.79253, 2.80998, 2.82743, 2.84489, 2.86234, 2.87979, 2.89725, 2.9147, 2.93215, 2.94961, 2.96706, 2.98451, 3.00197, 3.01942, 3.03687, 3.05433, 3.07178, 3.08923, 3.10669, 3.12414, 3.14159, 3.15905, 3.1765, 3.19395, 3.21141, 3.22886, 3.24631, 3.26377, 3.28122, 3.29867, 3.31613, 3.33358, 3.35103, 3.36849, 3.38594, 3.40339, 3.42085, 3.4383, 3.45575, 3.47321, 3.49066, 3.50811, 3.52557, 3.54302, 3.56047, 3.57792, 3.59538, 3.61283, 3.63028, 3.64774, 3.66519, 3.68264, 3.7001, 3.71755, 3.735, 3.75246, 3.76991, 3.78736, 3.80482, 3.82227, 3.83972, 3.85718, 3.87463, 3.89208, 3.90954, 3.92699, 3.94444, 3.9619, 3.97935, 3.9968, 4.01426, 4.03171, 4.04916, 4.06662, 4.08407, 4.10152, 4.11898, 4.13643, 4.15388, 4.17134, 4.18879, 4.20624, 4.2237, 4.24115, 4.2586, 4.27606, 4.29351, 4.31096, 4.32842, 4.34587, 4.36332, 4.38078, 4.39823, 4.41568, 4.43314, 4.45059, 4.46804, 4.4855, 4.50295, 4.5204, 4.53786, 4.55531, 4.57276, 4.59022, 4.60767, 4.62512, 4.64258, 4.66003, 4.67748, 4.69494, 4.71239, 4.72984, 4.7473, 4.76475, 4.7822, 4.79966, 4.81711, 4.83456, 4.85202, 4.86947, 4.88692, 4.90438, 4.92183, 4.93928, 4.95674, 4.97419, 4.99164, 5.00909, 5.02655, 5.044, 5.06145, 5.07891, 5.09636, 5.11381, 5.13127, 5.14872, 5.16617, 5.18363, 5.20108, 5.21853, 5.23599, 5.25344, 5.27089, 5.28835, 5.3058, 5.32325, 5.34071, 5.35816, 5.37561, 5.39307, 5.41052, 5.42797, 5.44543, 5.46288, 5.48033, 5.49779, 5.51524, 5.53269, 5.55015, 5.5676, 5.58505, 5.60251, 5.61996, 5.63741, 5.65487, 5.67232, 5.68977, 5.70723, 5.72468, 5.74213, 5.75959, 5.77704, 5.79449, 5.81195, 5.8294, 5.84685, 5.86431, 5.88176, 5.89921, 5.91667, 5.93412, 5.95157, 5.96903, 5.98648, 6.00393, 6.02139, 6.03884, 6.05629, 6.07375, 6.0912, 6.10865, 6.12611, 6.14356, 6.16101, 6.17847, 6.19592, 6.21337, 6.23083, 6.24828, 6.26573, 6.28319
};

std::array<std::array<real_t, 4>, 6> constexpr indicesCoords =
	{ std::array<real_t,4>{{0, 4, 1, 5}}, {0, 4, 2, 5}, {0, 4, 3, 5}, {2, 4, 3, 5}, {2, 4, 1, 5}, {3, 4, 1, 5} };

}




/////////////////////////////////////////////////////////
/////////////////// ModelBehavior /////////////////////// {{{1
/////////////////////////////////////////////////////////
//
//void ModelBehavior::_computeFishPDF() {
//	_fishPDF = vector<real_t>(evalAnglesNb, 0.0);
//	real_t fishPDFNormalization = 0.0;
//
//	for(auto& _a: _simulation.agents) {
//		Agent* a = _a.first.get();
//		if(a == _agent)
//			continue;
//
//		real_t const relativeAgentsHeadsPosX = a->headPos.first - _agent->headPos.first;
//		real_t const relativeAgentsHeadsPosY = a->headPos.second - _agent->headPos.second;
//		real_t const relativeAgentsTailPosX = a->tailPos.first - _agent->tailPos.first;
//		real_t const relativeAgentsTailPosY = a->tailPos.second - _agent->tailPos.second;
//		real_t const headRetinaX = relativeAgentsHeadsPosX * cos(_agent->direction) + relativeAgentsHeadsPosY * sin(_agent->direction);
//		real_t const headRetinaY = -relativeAgentsHeadsPosX * sin(_agent->direction) + relativeAgentsHeadsPosY * cos(_agent->direction);
//		real_t const tailRetinaX = relativeAgentsTailPosX * cos(_agent->direction) + relativeAgentsTailPosY * sin(_agent->direction);
//		real_t const tailRetinaY = -relativeAgentsTailPosX * sin(_agent->direction) + relativeAgentsTailPosY * cos(_agent->direction);
//		real_t const headAngle = normAngle(atan2(headRetinaY, headRetinaX));
//		real_t const distAgent = sqrt(headRetinaX * headRetinaX + headRetinaY * headRetinaY);
//		real_t const tailAngle = normAngle(atan2(tailRetinaY, tailRetinaX));
//		real_t const agentCenter = normAngle(atan2(cos(headAngle) + cos(tailAngle), sin(headAngle) + sin(tailAngle)));
//		real_t const diffDirections = normAngle(a->direction - _agent->direction);
//		real_t agentLength = abs(normAngle(headAngle) - normAngle(tailAngle));
//		if(agentLength > pi())
//			agentLength = 2.0 * pi() - agentLength;
//
//		// Find if the agent is within current agent 'a' FOV
//		if(		(distAgent < _agent->agentDistancePerception) &&
//				((agentCenter < _agent->fov * pi() / 180.0) || (agentCenter > (360. - _agent->fov) * pi() / 180.0))) {
//			fishPDFNormalization += agentLength;
//			// Compute PDF
//			for(size_t k = 0; k < evalAnglesNb; ++k) {
//				_fishPDF[k] += _fishPDFBessel * exp(kappaFishes * cos(evaluedProb[k] - agentCenter));
//				//_fishPDF[k] += _fishPDFBessel * exp(kappaFishes * cos(evaluedProb[k] - diffDirections));
//			}
//		}
//	}
//
//	// Normalize PDF
//	if(fishPDFNormalization > 0.0) {
//		_alphas = neighborsInfluence;
//		normalizeVec(_fishPDF, fishPDFNormalization);
//	} else {
//		_alphas = 0.0;
//	}
//}
//
//
//void ModelBehavior::_computeWallPDF() {
//	_wallPDF = vector<real_t>(evalAnglesNb, 0.0);
//	real_t wallPDFNormalization = 0.0;
//
//	std::vector<real_t> const impactLim({_agent->headPos.second, _agent->headPos.first, _agent->headPos.second, _agent->headPos.first});
//	std::vector<real_t> const wallXPosTmp({0.0 - _agent->headPos.first, _agent->headPos.first - _agent->headPos.first,
//			_simulation.arena.size.first - _agent->headPos.first, _agent->headPos.first - _agent->headPos.first});
//	std::vector<real_t> const wallYPosTmp({_agent->headPos.second - _agent->headPos.second, _simulation.arena.size.second - _agent->headPos.second,
//			_agent->headPos.second - _agent->headPos.second, 0.0 - _agent->headPos.first});
//
//	for(size_t k = 0; k < 4; ++k) {
//		real_t const wallRetinaX = wallXPosTmp[k] * cos(_agent->direction) + wallYPosTmp[k] * sin(_agent->direction);
//		real_t const wallRetinaY = -wallXPosTmp[k] * sin(_agent->direction) + wallYPosTmp[k] * cos(_agent->direction);
//		real_t const wallCenterAngleTmp = atan2(wallRetinaY, wallRetinaX);
//		real_t const distWallTmp = sqrt(wallRetinaX * wallRetinaX + wallRetinaY * wallRetinaY);
//		real_t const wallLength = distWallTmp * 100.0;
//		if(distWallTmp > _agent->wallDistancePerception) {
//			real_t const wallCenterAngle = normAngle(wallCenterAngleTmp);
//			real_t const distWall = distWallTmp;
//			real_t const impactPos = impactLim[k];
//			real_t const wallVal = 1.0;
//			if((distWall <= _agent->wallDistancePerception && impactPos >= 0 && impactPos <= _simulation.arena.size.first) &&
//					(wallCenterAngle < _agent->fov * pi() / 180.0 || wallCenterAngle > 360.0 - _agent->fov * pi() / 180.0 )) {
//				++wallPDFNormalization;
//				for(size_t j = 0; j < evalAnglesNb; ++j) {
//					_wallPDF[j] += wallVal / 2.0 * pi() * boost::math::cyl_bessel_i(0.0, wallLength) * exp(wallLength * cos(evaluedProb[j] - wallCenterAngle));
//				}
//			}
//		} else {
//			{
//				real_t const wallCenterAngle = normAngle(wallCenterAngleTmp + abs(acos(distWallTmp / _agent->wallDistancePerception)));
//				real_t const impactPos = impactLim[k] + pow(-1.0, ceil(k / 2.0)) * abs(sin(acos(distWallTmp / _agent->wallDistancePerception)) * _agent->wallDistancePerception);
//				real_t const wallVal = 0.5;
//				if((impactPos >= 0 && impactPos <= _simulation.arena.size.first) &&
//						(wallCenterAngle < _agent->fov * pi() / 180.0 || wallCenterAngle > 360.0 - _agent->fov * pi() / 180.0 )) {
//					++wallPDFNormalization;
//					for(size_t j = 0; j < evalAnglesNb; ++j) {
//						_wallPDF[j] += wallVal / 2.0 * pi() * boost::math::cyl_bessel_i(0.0, wallLength) * exp(wallLength * cos(evaluedProb[j] - wallCenterAngle));
//					}
//				}
//			}
//			{
//				real_t const wallCenterAngle = normAngle(wallCenterAngleTmp - abs(acos(distWallTmp / _agent->wallDistancePerception)));
//				real_t const impactPos = impactLim[k] + pow(-1.0, ceil(k / 2.0) + 1.0) * abs(sin(acos(distWallTmp / _agent->wallDistancePerception)) * _agent->wallDistancePerception);
//				real_t const wallVal = 0.5;
//				if((impactPos >= 0 && impactPos <= _simulation.arena.size.first) &&
//						(wallCenterAngle < _agent->fov * pi() / 180.0 || wallCenterAngle > 360.0 - _agent->fov * pi() / 180.0 )) {
//					++wallPDFNormalization;
//					for(size_t j = 0; j < evalAnglesNb; ++j) {
//						_wallPDF[j] += wallVal / 2.0 * pi() * boost::math::cyl_bessel_i(0.0, wallLength) * exp(wallLength * cos(evaluedProb[j] - wallCenterAngle));
//					}
//				}
//			}
//		}
//	}
//
//	// Normalize PDF
//	if(wallPDFNormalization > 0.0) {
//		_gammas = wallsInfluence;
//		normalizeVec(_wallPDF, wallPDFNormalization);
//	} else {
//		_gammas = 0.0;
//	}
//}
//
//
//real_t ModelBehavior::_computeAgentSpeed() {
//	//static lognormal_distribution<real_t> speedDistrib(_meanSpeed, _varSpeed);
//	lognormal_distribution<real_t> speedDistrib(_agent->meanSpeed, _agent->varSpeed);
//	return speedDistrib(rne());
//}
//
//
//void ModelBehavior::reinit() {
//	_fishPDFBessel = 1.0 / (2.0 * pi() * boost::math::cyl_bessel_i(0.0, kappaFishes));
//
//	for(size_t k = 0; k < evalAnglesNb; ++k) {
//		_neutPDF[k] = (1.0 / (2.0 * pi() * boost::math::cyl_bessel_i(0.0, kappaNeut))) * exp(kappaNeut * cos(evaluedProb[k]));
//		_neutPDF2[k] = (inertia / (2.0 * pi() * boost::math::cyl_bessel_i(0.0, kappaNeut))) * exp(kappaNeut * cos(evaluedProb[k]));
//	}
//}
//
//
//void ModelBehavior::step() {
//	// Compute PDFs
//	_computeFishPDF();
//	_computeWallPDF();
//
//	// Compute final PDF and CDF
//	real_t cdfsum = 0.0;
//	if(_alphas + _gammas > 0.0) {
//		for(size_t k = 0; k < evalAnglesNb; ++k) {
//			//real_t const val = (_inertia / (2.0 * pi() * boost::math::cyl_bessel_i(0.0, 2.3)) * exp(2.3 * cos(evaluedProb[k])) + _agent->perceptionInfluence *
//			//		(_alphas * _fishPDF[k] + _gammas * _wallPDF[k]) / (_alphas + _gammas))
//			//	/ (_inertia + _agent->perceptionInfluence);
//			real_t const val = (_neutPDF2[k] + perceptionInfluence * (_alphas * _fishPDF[k] + _gammas * _wallPDF[k]) / (_alphas + _gammas))
//				/ (inertia + perceptionInfluence);
//			_PDF[k] = val;
//			_CDF[k] = cdfsum + val;
//			cdfsum += val;
//		}
//	} else {
//		for(size_t k = 0; k < evalAnglesNb; ++k) {
//			real_t const val = _neutPDF[k];
//			_PDF[k] = val;
//			_CDF[k] = cdfsum + val;
//			cdfsum += val;
//		}
//	}
//
//	// Normalize CDF
//	normalizeVec(_CDF, cdfsum);
//
//	//for(size_t run = 0; run < 50; ++run) {
//	for(;;) {
//		real_t newSpeed = _computeAgentSpeed();
//		real_t z = unif01(rne());
//		int dir = -1;
//		for(int k = 0; k < _CDF.size(); ++k) {
//			if(z <= _CDF[k]) {
//				dir = k;
//				break;
//			}
//		}
//
//		real_t newDirection = _agent->direction;
//		if(dir >= 0) {
//			newDirection += evaluedProb[dir];
//		}
//		newDirection = normAngle(newDirection);
//		Coord_t newAgentsHeadPos = {
//			_agent->headPos.first + _simulation.dt * newSpeed * cos(newDirection),
//			_agent->headPos.second + _simulation.dt * newSpeed * sin(newDirection)
//		};
//
//		if(_simulation.arena.isInArena(newAgentsHeadPos)) {
//			_agent->direction = newDirection;
//			_agent->speed = newSpeed;
//			_agent->headPos = newAgentsHeadPos;
//			break;
//		} else {
//			//if(dir >= 0) {
//			//	_CDF[dir] = 0.0
//			//}
//		}
//	}
//
//	_agent->updateAgentPosition(_simulation.dt);
//}
//


/////////////////////////////////////////////////////////
////////////////// Model3DBehavior/////////////////////// {{{1
/////////////////////////////////////////////////////////
//
//void Model3DBehavior::_computeFishPDF() {
//	_fishPDF = vector<real_t>(evalAnglesNb, 0.0);
//	//real_t fishPDFNormalization = 0.0;
//	_totalAreaFishes = 0.0;
//
//	for(auto& _a: _simulation.agents) {
//		Agent* a = _a.first.get();
//		if(a == _agent || !a->present)
//			continue;
//
//		real_t const relativeAgentsHeadsPosX = a->headPos.first - _agent->headPos.first;
//		real_t const relativeAgentsHeadsPosY = a->headPos.second - _agent->headPos.second;
//		std::array<Coord3D_t, 6> tmpFishPosition;
//		tmpFishPosition[0] = {
//			relativeAgentsHeadsPosX,
//			relativeAgentsHeadsPosY,
//			0.0};
//
//		tmpFishPosition[1] = {
//			relativeAgentsHeadsPosX - cos(a->direction) * _agent->length,
//			relativeAgentsHeadsPosY - sin(a->direction) * _agent->length,
//			0.0};
//
//		tmpFishPosition[2] = {
//			relativeAgentsHeadsPosX + cos(a->direction + pi() + _agent->angleFishSide) * _agent->normBody,
//			relativeAgentsHeadsPosY + sin(a->direction + pi() + _agent->angleFishSide) * _agent->normBody,
//			0.0};
//
//		tmpFishPosition[3] = {
//			relativeAgentsHeadsPosX + cos(a->direction + pi() - _agent->angleFishSide) * _agent->normBody,
//			relativeAgentsHeadsPosY + sin(a->direction + pi() - _agent->angleFishSide) * _agent->normBody,
//			0.0};
//
//		tmpFishPosition[4] = {
//			relativeAgentsHeadsPosX - cos(a->direction) * _agent->length / 3.0,
//			relativeAgentsHeadsPosY - sin(a->direction) * _agent->length / 3.0,
//			_agent->height / 2.0};
//
//		tmpFishPosition[5] = {
//			relativeAgentsHeadsPosX - cos(a->direction) * _agent->length / 3.0,
//			relativeAgentsHeadsPosY - sin(a->direction) * _agent->length / 3.0,
//			-_agent->height / 2.0};
//
//		real_t const relativeDirection = normAngle(a->direction - _agent->direction);
//
//		std::array<real_t, 6> angleRetina;
//		for(size_t k = 0; k < 6; ++k) {
//			real_t const positionXRetina = tmpFishPosition[k][0] * cos(_agent->direction) + tmpFishPosition[k][1] * sin(_agent->direction);
//			real_t const positionYRetina = -tmpFishPosition[k][0] * sin(_agent->direction) + tmpFishPosition[k][1] * cos(_agent->direction);
//			angleRetina[k] = atan2(positionYRetina, positionXRetina);
//		}
//
//		real_t const headRetinaX = relativeAgentsHeadsPosX * cos(_agent->direction) + relativeAgentsHeadsPosY * sin(_agent->direction);
//		real_t const headRetinaY = -relativeAgentsHeadsPosX * sin(_agent->direction) + relativeAgentsHeadsPosY * cos(_agent->direction);
//		real_t const headAngle = normAngle(atan2(headRetinaY, headRetinaX));
//		real_t const distAgent = sqrt(headRetinaX * headRetinaX + headRetinaY * headRetinaY);
//		real_t const agentCenter = normAngle(atan2(cos(angleRetina[0]) + cos(angleRetina[1]), sin(angleRetina[0]) + sin(angleRetina[1])) );
//		real_t const diffDirections = normAngle(a->direction - _agent->direction);
//
//		std::array<real_t, 6> agentSize = {
//			std::abs(angleRetina[0] - angleRetina[1]),
//			std::abs(angleRetina[0] - angleRetina[2]),
//			std::abs(angleRetina[0] - angleRetina[3]),
//			std::abs(angleRetina[3] - angleRetina[2]),
//			std::abs(angleRetina[2] - angleRetina[1]),
//			std::abs(angleRetina[3] - angleRetina[1]) };
//		size_t indiceLargestSize = 0;
//		real_t largestSize = 0.0;
//		for(size_t k = 0; k < 6; ++k) {
//			if(agentSize[k] > pi())
//				agentSize[k] = std::abs(2. * pi() - agentSize[k]);
//			if(agentSize[k] > largestSize) {
//				largestSize = agentSize[k];
//				indiceLargestSize = k;
//			}
//		}
//
//		std::array<Coord3D_t, 4> fishesCoords;
//		for(size_t k = 0; k < 4; ++k) {
//			fishesCoords[k] = tmpFishPosition[indicesCoords[indiceLargestSize][k]];
//		}
//		auto const areaFish = surfaceSphereFish(fishesCoords);
//
//		// Find if the agent 'a' is within current agent FOV
//		if(		(areaFish > _agent->agentMinAnglePerception) &&
//				((agentCenter < _agent->fov * pi() / 180.0) || (agentCenter > (360. - _agent->fov) * pi() / 180.0))) {
//			//fishPDFNormalization += 1.0;
//			_totalAreaFishes += areaFish;
//			// Compute PDF
//			if(distAgent < repulsionFromAgentsAtDist) {
//				// Repulsion from agent if it is too close
//				for(size_t k = 0; k < evalAnglesNb; ++k) {
//					_fishPDF[k] += areaFish * _fishPDFBessel * exp(kappaFishes * cos(normAngle(evaluedProb[k] - headAngle + pi())));
//				}
//			} else {
//				for(size_t k = 0; k < evalAnglesNb; ++k) {
//					_fishPDF[k] += areaFish * _fishPDFBessel * exp(kappaFishes * cos(evaluedProb[k] - headAngle));
//					//_fishPDF[k] += areaFish * _fishPDFBessel * exp(kappaFishes * cos(evaluedProb[k] - diffDirections));
//					//_fishPDF[k] += _fishPDFBessel * exp(kappaFishes * cos(evaluedProb[k] - diffDirections));
//				}
//			}
//		}
//	}
//
//	// Normalize PDF
//	if(_totalAreaFishes > 0.0) {
//		for(size_t k = 0; k < evalAnglesNb; ++k) {
//			_fishPDF[k] /= _totalAreaFishes;
//		}
//	}
//}
//
//
//void Model3DBehavior::reinit() {
//	_fishPDFBessel = 1.0 / (2.0 * pi() * boost::math::cyl_bessel_i(0.0, kappaFishes));
//
//	for(size_t k = 0; k < evalAnglesNb; ++k) {
//		_neutPDFCenter[k] = exp(kappaNeutCenter * cos(evaluedProb[k])) / (2.0 * pi() * boost::math::cyl_bessel_i(0.0, kappaNeutCenter));
//		_neutPDFWall[k] = exp(kappaNeutWall * cos(evaluedProb[k])) / (2.0 * pi() * boost::math::cyl_bessel_i(0.0, kappaNeutWall));
//		//_neutPDF[k] = (1.0 / (2.0 * pi() * boost::math::cyl_bessel_i(0.0, kappaNeut))) * exp(kappaNeut * cos(evaluedProb[k]));
//		//_neutPDF2[k] = (inertia / (2.0 * pi() * boost::math::cyl_bessel_i(0.0, kappaNeut))) * exp(kappaNeut * cos(evaluedProb[k]));
//	}
//
//	//ModelBehavior::reinit();
//	//_wallPDFBessel = 1.0 / (2.0 * pi() * boost::math::cyl_bessel_i(0.0, kappaWalls));
//}
//
//
//// INFO assume square arena
//void Model3DBehavior::_computeDistToWall() {
//	_distToWall =  std::numeric_limits<double>::max();
//	_distToWall = std::min(_distToWall, std::abs(_agent->headPos.first - 0.0));
//	_distToWall = std::min(_distToWall, std::abs(_agent->headPos.first - 1.0));
//	_distToWall = std::min(_distToWall, std::abs(_agent->headPos.second - 0.0));
//	_distToWall = std::min(_distToWall, std::abs(_agent->headPos.second - 1.0));
//
////	if(std::abs(_agent->headPos.first - 0.0) < thresholdWalls || std::abs(_agent->headPos.first - 1.0) < thresholdWalls) {
////		for(size_t k = 0; k < evalAnglesNb; ++k) {
////			_neutPDFWall[k] = (exp(kappaNeutWall * cos(evaluedProb[k] - pi() / 2.0)) + exp(kappaNeutWall * cos(evaluedProb[k] - pi() * 3.0 / 2.0))) / (2.0 * pi() * boost::math::cyl_bessel_i(0.0, kappaNeutWall)) / 2.0;
////		}
////	} else if(std::abs(_agent->headPos.second - 0.0) < thresholdWalls || std::abs(_agent->headPos.second - 1.0) < thresholdWalls) {
////		for(size_t k = 0; k < evalAnglesNb; ++k) {
////			_neutPDFWall[k] = (exp(kappaNeutWall * cos(evaluedProb[k])) + exp(kappaNeutWall * cos(evaluedProb[k] - pi()))) / (2.0 * pi() * boost::math::cyl_bessel_i(0.0, kappaNeutWall)) / 2.0;
////		}
////	}
//
//	//_distToWall =  std::numeric_limits<double>::max();
//	_neutPDFWall = vector<real_t>(evalAnglesNb, 0.0);
//	real_t wallPDFNormalization = 0.0;
//
//	//std::vector<real_t> const impactLim({_agent->headPos.second, _agent->headPos.first, _agent->headPos.second, _agent->headPos.first});
//	std::vector<real_t> const wallXPosTmp({0.0 - _agent->headPos.first, _agent->headPos.first - _agent->headPos.first,
//			_simulation.arena.size.first - _agent->headPos.first, _agent->headPos.first - _agent->headPos.first});
//	std::vector<real_t> const wallYPosTmp({_agent->headPos.second - _agent->headPos.second, _simulation.arena.size.second - _agent->headPos.second,
//			_agent->headPos.second - _agent->headPos.second, 0.0 - _agent->headPos.second});
//
//	std::vector<real_t> wallRetinaX;
//	std::vector<real_t> wallRetinaY;
//	std::vector<real_t> wallCenterAngleTmp;
//	std::vector<real_t> distWallTmp;
//	for(size_t k = 0; k < wallXPosTmp.size(); ++k) {
//		wallRetinaX.push_back(wallXPosTmp[k] * cos(_agent->direction) + wallYPosTmp[k] * sin(_agent->direction));
//		wallRetinaY.push_back(-wallXPosTmp[k] * sin(_agent->direction) + wallYPosTmp[k] * cos(_agent->direction));
//		wallCenterAngleTmp.push_back(atan2(wallRetinaY[k], wallRetinaX[k]));
//		distWallTmp.push_back(sqrt(wallRetinaX[k] * wallRetinaX[k] + wallRetinaY[k] * wallRetinaY[k]));
//	}
//
//	std::vector<real_t> anglesWall;
//	if(distWallTmp[0] < _agent->wallDistancePerception && distWallTmp[1] > _agent->wallDistancePerception && distWallTmp[3] > _agent->wallDistancePerception) {
//		anglesWall = {wallCenterAngleTmp[0] + pi() / 2., wallCenterAngleTmp[0] - pi() / 2.};
//	} else if(distWallTmp[0] < _agent->wallDistancePerception && distWallTmp[1] < _agent->wallDistancePerception) {
//		anglesWall = {wallCenterAngleTmp[0] + pi() / 2., wallCenterAngleTmp[1] - pi() / 2.};
//	} else if(distWallTmp[1] < _agent->wallDistancePerception && distWallTmp[0] > _agent->wallDistancePerception && distWallTmp[2] > _agent->wallDistancePerception) {
//		anglesWall = {wallCenterAngleTmp[1] + pi() / 2., wallCenterAngleTmp[1] - pi() / 2.};
//	} else if(distWallTmp[1] < _agent->wallDistancePerception && distWallTmp[2] < _agent->wallDistancePerception) {
//		anglesWall = {wallCenterAngleTmp[1] + pi() / 2., wallCenterAngleTmp[2] - pi() / 2.};
//	} else if(distWallTmp[2] < _agent->wallDistancePerception && distWallTmp[1] > _agent->wallDistancePerception && distWallTmp[3] > _agent->wallDistancePerception) {
//		anglesWall = {wallCenterAngleTmp[2] + pi() / 2., wallCenterAngleTmp[2] - pi() / 2.};
//	} else if(distWallTmp[2] < _agent->wallDistancePerception && distWallTmp[3] < _agent->wallDistancePerception) {
//		anglesWall = {wallCenterAngleTmp[2] + pi() / 2., wallCenterAngleTmp[3] - pi() / 2.};
//	} else if(distWallTmp[3] < _agent->wallDistancePerception && distWallTmp[2] > _agent->wallDistancePerception && distWallTmp[0] > _agent->wallDistancePerception) {
//		anglesWall = {wallCenterAngleTmp[3] + pi() / 2., wallCenterAngleTmp[3] - pi() / 2.};
//	} else if(distWallTmp[3] < _agent->wallDistancePerception && distWallTmp[0] < _agent->wallDistancePerception) {
//		anglesWall = {wallCenterAngleTmp[3] + pi() / 2., wallCenterAngleTmp[0] - pi() / 2.};
//	} else {
//		anglesWall = {};
//	}
//
//	for(size_t k = 0; k < anglesWall.size(); ++k) {
//		real_t const wallCenterAngle = normAngle(anglesWall[k]);
//		if(wallCenterAngle < _agent->fov * pi() / 180.0 || wallCenterAngle > 360.0 - _agent->fov * pi() / 180.0) {
//			for(size_t j = 0; j < evalAnglesNb; ++j) {
//				_neutPDFWall[j] += 1.0 / 2.0 * pi() * boost::math::cyl_bessel_i(0.0, kappaNeutWall) * exp(kappaNeutWall * cos(evaluedProb[j] - wallCenterAngle));
//			}
//			++wallPDFNormalization;
//		}
//	}
//
////	for(size_t k = 0; k < 4; ++k) {
////		real_t const wallRetinaX = wallXPosTmp[k] * cos(_agent->direction) + wallYPosTmp[k] * sin(_agent->direction);
////		real_t const wallRetinaY = -wallXPosTmp[k] * sin(_agent->direction) + wallYPosTmp[k] * cos(_agent->direction);
////		real_t const wallCenterAngleTmp = atan2(wallRetinaY, wallRetinaX);
////		real_t const distWallTmp = sqrt(wallRetinaX * wallRetinaX + wallRetinaY * wallRetinaY);
////		//_distToWall = std::min(_distToWall, distWallTmp);
////		if(distWallTmp > _agent->wallDistancePerception) {
////			real_t const wallCenterAngle = normAngle(wallCenterAngleTmp);
////			real_t const distWall = distWallTmp;
////			real_t const impactPos = impactLim[k];
////			if((distWall <= _agent->wallDistancePerception && impactPos >= 0 && impactPos <= _simulation.arena.size.first) &&
////					(wallCenterAngle < _agent->fov * pi() / 180.0 || wallCenterAngle > 360.0 - _agent->fov * pi() / 180.0 )) {
////				++wallPDFNormalization;
////				for(size_t j = 0; j < evalAnglesNb; ++j) {
////					_neutPDFWall[j] += 1.0 / 2.0 * pi() * boost::math::cyl_bessel_i(0.0, kappaNeutWall) * exp(kappaNeutWall * cos(evaluedProb[j] - wallCenterAngle));
////				}
////			}
////		} else {
////			{
////				real_t const wallCenterAngle = normAngle(wallCenterAngleTmp + abs(acos(distWallTmp / _agent->wallDistancePerception)));
////				real_t const impactPos = impactLim[k] + pow(-1.0, ceil(k / 2.0)) * abs(sin(acos(distWallTmp / _agent->wallDistancePerception)) * _agent->wallDistancePerception);
////				if((impactPos >= 0 && impactPos <= _simulation.arena.size.first) &&
////						(wallCenterAngle < _agent->fov * pi() / 180.0 || wallCenterAngle > 360.0 - _agent->fov * pi() / 180.0 )) {
////					++wallPDFNormalization;
////					for(size_t j = 0; j < evalAnglesNb; ++j) {
////						_neutPDFWall[j] += 1.0 / 2.0 * pi() * boost::math::cyl_bessel_i(0.0, kappaNeutWall) * exp(kappaNeutWall * cos(evaluedProb[j] - wallCenterAngle));
////					}
////				}
////			}
////			{
////				real_t const wallCenterAngle = normAngle(wallCenterAngleTmp - abs(acos(distWallTmp / _agent->wallDistancePerception)));
////				real_t const impactPos = impactLim[k] + pow(-1.0, ceil(k / 2.0) + 1.0) * abs(sin(acos(distWallTmp / _agent->wallDistancePerception)) * _agent->wallDistancePerception);
////				if((impactPos >= 0 && impactPos <= _simulation.arena.size.first) &&
////						(wallCenterAngle < _agent->fov * pi() / 180.0 || wallCenterAngle > 360.0 - _agent->fov * pi() / 180.0 )) {
////					++wallPDFNormalization;
////					for(size_t j = 0; j < evalAnglesNb; ++j) {
////						_neutPDFWall[j] += 1.0 / 2.0 * pi() * boost::math::cyl_bessel_i(0.0, kappaNeutWall) * exp(kappaNeutWall * cos(evaluedProb[j] - wallCenterAngle));
////					}
////				}
////			}
////		}
////	}
//
//	// Normalize PDF
//	if(wallPDFNormalization > 0.0) {
//		_gammas = wallsInfluence;
//		normalizeVec(_neutPDFWall, wallPDFNormalization);
//	} else {
//		_gammas = 0.0;
//	}
//
//}
//
//void Model3DBehavior::step() {
//	_computeFishPDF();
//	_computeDistToWall();
//
//	// Compute final PDF and CDF
//	real_t cdfsum = 0.0;
////	if(_distToWall < _agent->wallDistancePerception) {
////	//if(_distToWall < thresholdWalls) {
////		for(size_t k = 0; k < evalAnglesNb; ++k) {
////			real_t const val = (_neutPDFWall[k] + alphasWalls * _totalAreaFishes * _fishPDF[k]) / (1. + alphasWalls * _totalAreaFishes);
////			_PDF[k] = val;
////			_CDF[k] = cdfsum + val;
////			cdfsum += val;
////		}
////	} else {
//		for(size_t k = 0; k < evalAnglesNb; ++k) {
//			real_t const val = (_neutPDFCenter[k] + alphasCenter * _totalAreaFishes * _fishPDF[k]) / (1. + alphasCenter * _totalAreaFishes);
//			_PDF[k] = val;
//			_CDF[k] = cdfsum + val;
//			cdfsum += val;
//		}
////	}
//
//	// Normalize CDF
//	normalizeVec(_CDF, cdfsum);
//
//	for(size_t run = 0; run < 50; ++run) {
//	//for(;;) {
//		real_t newSpeed = _computeAgentSpeed();
//		real_t z = unif01(rne());
//		int dir = -1;
//		for(int k = 0; k < _CDF.size(); ++k) {
//			if(z <= _CDF[k]) {
//				dir = k;
//				break;
//			}
//		}
//
//		real_t newDirection = _agent->direction;
//		real_t deltaDirection = 0.0;
//		if(dir >= 0) {
//			//newDirection += evaluedProb[dir];
//			deltaDirection = evaluedProb[dir];
//		}
//
//		// XXX Put in agent's code
//		// Adjust newDirection according to agent's maxTurningRate
//		if(deltaDirection < pi()) {
//			if(deltaDirection > _agent->maxTurningRate) {
//				deltaDirection = _agent->maxTurningRate;
//			}
//		} else {
//			if(deltaDirection <= 2. * pi() - _agent->maxTurningRate) {
//				deltaDirection = 2. * pi() - _agent->maxTurningRate;
//			}
//		}
//
//		newDirection = normAngle(newDirection + deltaDirection);
//
//		Coord_t newAgentsHeadPos = {
//			_agent->headPos.first + _simulation.dt * newSpeed * cos(newDirection),
//			_agent->headPos.second + _simulation.dt * newSpeed * sin(newDirection)
//		};
//
//		if(_simulation.arena.isInArena(newAgentsHeadPos)) {
//			_agent->direction = newDirection;
//			_agent->speed = newSpeed;
//			_agent->headPos = newAgentsHeadPos;
//			break;
//		} else {
//			//if(dir >= 0) {
//			//	_CDF[dir] = 0.0
//			//}
//		}
//	}
//
//	_agent->updateAgentPosition(_simulation.dt);
//}



/////////////////////////////////////////////////////////
//////////////////////// BM ///////////////////////////// {{{1
/////////////////////////////////////////////////////////

void BM::_computeFishPDF() {
	_fishPDF = vector<real_t>(evalAnglesNb, 0.0);
	//real_t fishPDFNormalization = 0.0;
	_totalAreaFishes = 0.0;

	for(auto& _a: _simulation.agents) {
		Agent* a = _a.first.get();
		if(a == _agent || !a->present)
			continue;

		real_t const relativeAgentsHeadsPosX = a->headPos.first - _agent->headPos.first;
		real_t const relativeAgentsHeadsPosY = a->headPos.second - _agent->headPos.second;
		std::array<Coord3D_t, 6> tmpFishPosition;
		tmpFishPosition[0] = {
			relativeAgentsHeadsPosX,
			relativeAgentsHeadsPosY,
			0.0};

		tmpFishPosition[1] = {
			relativeAgentsHeadsPosX - cos(a->direction) * _agent->length,
			relativeAgentsHeadsPosY - sin(a->direction) * _agent->length,
			0.0};

		tmpFishPosition[2] = {
			relativeAgentsHeadsPosX + cos(a->direction + pi() + _agent->angleFishSide) * _agent->normBody,
			relativeAgentsHeadsPosY + sin(a->direction + pi() + _agent->angleFishSide) * _agent->normBody,
			0.0};

		tmpFishPosition[3] = {
			relativeAgentsHeadsPosX + cos(a->direction + pi() - _agent->angleFishSide) * _agent->normBody,
			relativeAgentsHeadsPosY + sin(a->direction + pi() - _agent->angleFishSide) * _agent->normBody,
			0.0};

		tmpFishPosition[4] = {
			relativeAgentsHeadsPosX - cos(a->direction) * _agent->length / 3.0,
			relativeAgentsHeadsPosY - sin(a->direction) * _agent->length / 3.0,
			_agent->height / 2.0};

		tmpFishPosition[5] = {
			relativeAgentsHeadsPosX - cos(a->direction) * _agent->length / 3.0,
			relativeAgentsHeadsPosY - sin(a->direction) * _agent->length / 3.0,
			-_agent->height / 2.0};

		real_t const relativeDirection = normAngle(a->direction - _agent->direction);

		std::array<real_t, 6> angleRetina;
		for(size_t k = 0; k < 6; ++k) {
			real_t const positionXRetina = tmpFishPosition[k][0] * cos(_agent->direction) + tmpFishPosition[k][1] * sin(_agent->direction);
			real_t const positionYRetina = -tmpFishPosition[k][0] * sin(_agent->direction) + tmpFishPosition[k][1] * cos(_agent->direction);
			angleRetina[k] = atan2(positionYRetina, positionXRetina);
		}

		real_t const headRetinaX = relativeAgentsHeadsPosX * cos(_agent->direction) + relativeAgentsHeadsPosY * sin(_agent->direction);
		real_t const headRetinaY = -relativeAgentsHeadsPosX * sin(_agent->direction) + relativeAgentsHeadsPosY * cos(_agent->direction);
		real_t const headAngle = normAngle(atan2(headRetinaY, headRetinaX));
		real_t const distAgent = sqrt(headRetinaX * headRetinaX + headRetinaY * headRetinaY);
		real_t const agentCenter = normAngle(atan2(cos(angleRetina[0]) + cos(angleRetina[1]), sin(angleRetina[0]) + sin(angleRetina[1])) );
		real_t const diffDirections = normAngle(a->direction - _agent->direction);

		std::array<real_t, 6> agentSize = {
			std::abs(angleRetina[0] - angleRetina[1]),
			std::abs(angleRetina[0] - angleRetina[2]),
			std::abs(angleRetina[0] - angleRetina[3]),
			std::abs(angleRetina[3] - angleRetina[2]),
			std::abs(angleRetina[2] - angleRetina[1]),
			std::abs(angleRetina[3] - angleRetina[1]) };
		size_t indiceLargestSize = 0;
		real_t largestSize = 0.0;
		for(size_t k = 0; k < 6; ++k) {
			if(agentSize[k] > pi())
				agentSize[k] = std::abs(2. * pi() - agentSize[k]);
			if(agentSize[k] > largestSize) {
				largestSize = agentSize[k];
				indiceLargestSize = k;
			}
		}

		std::array<Coord3D_t, 4> fishesCoords;
		for(size_t k = 0; k < 4; ++k) {
			fishesCoords[k] = tmpFishPosition[indicesCoords[indiceLargestSize][k]];
		}
		auto const areaFish = surfaceSphereFish(fishesCoords); //! A_{f_i}

		// Find if the agent 'a' is within current agent FOV
		if(		(areaFish > _agent->agentMinAnglePerception) &&
				((agentCenter < _agent->fov * pi() / 180.0) || (agentCenter > (360. - _agent->fov) * pi() / 180.0))) {
			//fishPDFNormalization += 1.0;
			_totalAreaFishes += areaFish;
			// Compute PDF
			if(distAgent < repulsionFromAgentsAtDist) {
				// Repulsion from agent if it is too close
				for(size_t k = 0; k < evalAnglesNb; ++k) {
					_fishPDF[k] += areaFish * _fishPDFBessel * exp(kappaFishes * cos(normAngle(evaluedProb[k] - headAngle + pi())));
				}
			} else {
				for(size_t k = 0; k < evalAnglesNb; ++k) {
					_fishPDF[k] += areaFish * _fishPDFBessel * exp(kappaFishes * cos(evaluedProb[k] - headAngle));
					//_fishPDF[k] += areaFish * _fishPDFBessel * exp(kappaFishes * cos(evaluedProb[k] - diffDirections));
					//_fishPDF[k] += _fishPDFBessel * exp(kappaFishes * cos(evaluedProb[k] - diffDirections));
				}
			}
		}
	}

	// Normalize PDF
	if(_totalAreaFishes > 0.0) {
		for(size_t k = 0; k < evalAnglesNb; ++k) {
			_fishPDF[k] /= _totalAreaFishes;
		}
	}
}


void BM::reinit() {
	_fishPDFBessel = 1.0 / (2.0 * pi() * boost::math::cyl_bessel_i(0.0, kappaFishes));

	for(size_t k = 0; k < evalAnglesNb; ++k) {
		_neutPDFCenter[k] = exp(kappaNeutCenter * cos(evaluedProb[k])) / (2.0 * pi() * boost::math::cyl_bessel_i(0.0, kappaNeutCenter));
	}
}



real_t BM::_computeAgentSpeed() {
	//static lognormal_distribution<real_t> speedDistrib(_meanSpeed, _varSpeed);
	lognormal_distribution<real_t> speedDistrib(_agent->meanSpeed, _agent->varSpeed);
	return speedDistrib(rne());
}

void BM::step() {
	_computeFishPDF();
	//_computeDistToWall();

	// Compute final PDF and CDF
	real_t cdfsum = 0.0;
		for(size_t k = 0; k < evalAnglesNb; ++k) {
			real_t const val = (_neutPDFCenter[k] + alphasCenter * _totalAreaFishes * _fishPDF[k]) / (1. + alphasCenter * _totalAreaFishes);
			_PDF[k] = val;
			_CDF[k] = cdfsum + val;
			cdfsum += val;
		}

	// Normalize CDF
	normalizeVec(_CDF, cdfsum);

	for(size_t run = 0; run < 50; ++run) {
	//for(;;) {
		real_t newSpeed = _computeAgentSpeed();
		real_t z = unif01(rne());
		int dir = -1;
		for(int k = 0; k < _CDF.size(); ++k) {
			if(z <= _CDF[k]) {
				dir = k;
				break;
			}
		}

		real_t newDirection = _agent->direction;
		real_t deltaDirection = 0.0;
		if(dir >= 0) {
			//newDirection += evaluedProb[dir];
			deltaDirection = evaluedProb[dir];
		}

		// XXX Put in agent's code
		// Adjust newDirection according to agent's maxTurningRate
		if(deltaDirection < pi()) {
			if(deltaDirection > _agent->maxTurningRate) {
				deltaDirection = _agent->maxTurningRate;
			}
		} else {
			if(deltaDirection <= 2. * pi() - _agent->maxTurningRate) {
				deltaDirection = 2. * pi() - _agent->maxTurningRate;
			}
		}

		newDirection = normAngle(newDirection + deltaDirection);

		Coord_t newAgentsHeadPos = {
			_agent->headPos.first + _simulation.dt * newSpeed * cos(newDirection),
			_agent->headPos.second + _simulation.dt * newSpeed * sin(newDirection)
		};

		if(_simulation.arena.isInArena(newAgentsHeadPos)) {
			_agent->direction = newDirection;
			_agent->speed = newSpeed;
			_agent->headPos = newAgentsHeadPos;
			break;
		} else {
			//if(dir >= 0) {
			//	_CDF[dir] = 0.0;
			//}
		}
	}

	_agent->updateAgentPosition(_simulation.dt);
}






/////////////////////////////////////////////////////////
////////////////////// ZonedBM ////////////////////////// {{{1
/////////////////////////////////////////////////////////


void ZonedBM::reinit() {
	try {
		//std::cout << "DEBUG2: " << kappaFishes << " " << kappaNeutCenter << std::endl;
		_fishPDFBessel = 1.0 / (2.0 * pi() * boost::math::cyl_bessel_i(0.0, kappaFishes));
		_wallPDFBessel = 1.0 / (2.0 * pi() * boost::math::cyl_bessel_i(0.0, kappaWalls));

		for(size_t k = 0; k < evalAnglesNb; ++k) {
			_neutPDFCenter[k] = exp(kappaNeutCenter * cos(evaluedProb[k])) / (2.0 * pi() * boost::math::cyl_bessel_i(0.0, kappaNeutCenter));
		}

		std::vector<real_t> speedHistPos = {};
		for(size_t i = 0; i <= speedHistogram.size(); ++i) {
			real_t const val = minSpeed + (maxSpeed - minSpeed) / static_cast<real_t>(speedHistogram.size()) * static_cast<real_t>(i);
			speedHistPos.push_back(val);
		}
		_speedDistribution = std::piecewise_constant_distribution<>(speedHistPos.begin(), speedHistPos.end(), speedHistogram.begin());
	} catch(std::exception& e) {
		std::cout << "EXCEPTION: " << e.what() << " " << kappaNeutCenter << " " << kappaFishes << std::endl;
		throw e;
	}
}



real_t ZonedBM::_computeAgentSpeed() {
	//lognormal_distribution<real_t> speedDistrib(_agent->meanSpeed, _agent->varSpeed);
	//return speedDistrib(rne());
	return _speedDistribution(rne());
}


void ZonedBM::_detectZonesAroundAgent(real_t r) {
	_zonesSensors.clear();
	for(size_t k = 0; k < evalAnglesNb; ++k) {
		real_t newDirection = _agent->direction;
		real_t deltaDirection = evaluedProb[k];

		// Adjust newDirection according to agent's maxTurningRate
		if(deltaDirection < pi()) {
			if(deltaDirection > _agent->maxTurningRate) {
				deltaDirection = _agent->maxTurningRate;
			}
		} else {
			if(deltaDirection <= 2. * pi() - _agent->maxTurningRate) {
				deltaDirection = 2. * pi() - _agent->maxTurningRate;
			}
		}

		newDirection = normAngle(newDirection + deltaDirection);
		Coord_t pos = {
			_agent->headPos.first + r * cos(newDirection),
			_agent->headPos.second + r * sin(newDirection)
		};

		_zonesSensors.push_back(_zdb->zoneId(pos));
	}
}


void ZonedBM::_computeZonesPDF() {
	_zonesPDF.clear();
	real_t normalization = 0.0;

	for(size_t k = 0; k < evalAnglesNb; ++k) {
		if(_zonesSensors[k] == 0) {
			// Not in arena
			_zonesPDF.push_back(0.);
			continue;
		}

		real_t const val = _zonesAffinity[_zonesSensors[k]];
		_zonesPDF.push_back(val);
		normalization += val;
	}

	// Normalize PDF
	for(size_t k = 0; k < evalAnglesNb; ++k) {
		_zonesPDF[k] /= normalization;
	}
}



void ZonedBM::_computeWallsPDF() {
	_neutPDFCenter = vector<real_t>(evalAnglesNb, 0.0);

	// Find closest wall
	real_t minDist = 10000.0;
	size_t minIndex = 0;
	for(size_t i = 0; i < wallsCoord.size(); ++i) {
		real_t const fstX = wallsCoord[i].first.first;
		real_t const fstY = wallsCoord[i].first.second;
		real_t const sndX = wallsCoord[i].second.first;
		real_t const sndY = wallsCoord[i].second.second;
		real_t const dist = std::abs( (sndY - fstY) * _agent->headPos.first - (sndX - fstX) * _agent->headPos.second + sndX * fstY - sndY * fstX) / std::sqrt((sndY - fstY) * (sndY - fstY) + (sndX - fstX) * (sndX - fstX));
		if(dist < minDist) {
			minDist = dist;
			minIndex = i;
		}
	}

	// XXX
	// Find closest corner
	size_t minCornerDist = 0;
	size_t minCornerIndex = 0;
	for(size_t i = 0; i < wallsCoord.size(); ++i) {
		real_t const fstX = wallsCoord[i].first.first;
		real_t const fstY = wallsCoord[i].first.second;
		real_t dist = std::sqrt( (fstX - _agent->headPos.first) * (fstX - _agent->headPos.first) + (fstY - _agent->headPos.second) * (fstY - _agent->headPos.second));
		if(dist < minCornerDist) {
			minCornerDist = dist;
			minCornerIndex = i;
		}
	}
	size_t minCornerPrevIndex = 0;
	if(minCornerIndex == 0) {
		minCornerPrevIndex = wallsCoord.size() - 1;
	} else {
		minCornerPrevIndex = minCornerPrevIndex - 1;
	}


	// Find 2 possible directions
	real_t fstTheta = 0.;
	real_t sndTheta = 0.;
	if(minDist <= minCornerDist) {
		real_t const fstX = wallsDirectionCoord[minIndex].first.first;
		real_t const fstY = wallsDirectionCoord[minIndex].first.second;
		real_t const sndX = wallsDirectionCoord[minIndex].second.first;
		real_t const sndY = wallsDirectionCoord[minIndex].second.second;
		fstTheta = normAngle(atan2(sndY, sndX) - atan2(fstY, fstX));
		sndTheta = normAngle(atan2(fstY, fstX) - atan2(sndY, sndX));
	} else {
		real_t const fstX = wallsDirectionCoord[minCornerIndex].first.first;
		real_t const fstY = wallsDirectionCoord[minCornerIndex].first.second;
		real_t const sndX = wallsDirectionCoord[minCornerIndex].second.first;
		real_t const sndY = wallsDirectionCoord[minCornerIndex].second.second;
		real_t const prevfstX = wallsDirectionCoord[minCornerPrevIndex].first.first;
		real_t const prevfstY = wallsDirectionCoord[minCornerPrevIndex].first.second;
		real_t const prevsndX = wallsDirectionCoord[minCornerPrevIndex].second.first;
		real_t const prevsndY = wallsDirectionCoord[minCornerPrevIndex].second.second;

		fstTheta = normAngle(atan2(sndY, sndX) - atan2(fstY, fstX));
		sndTheta = normAngle(atan2(prevfstY, prevfstX) - atan2(prevsndY, prevsndX));
	}


	// Update PDF
	for(size_t k = 0; k < evalAnglesNb; ++k) {
		_neutPDFCenter[k] = (
			exp(kappaNeutCenter * cos(evaluedProb[k] - fstTheta)) / (2.0 * pi() * boost::math::cyl_bessel_i(0.0, kappaNeutCenter)) +
			exp(kappaNeutCenter * cos(evaluedProb[k] - sndTheta)) / (2.0 * pi() * boost::math::cyl_bessel_i(0.0, kappaNeutCenter))
			) / 2.;
	}

//	_wallsPDF = vector<real_t>(evalAnglesNb, 0.0);
//	_totalAreaWalls = 0.0;
//
//	for(auto& _w: wallsCoord) {
//		real_t const fstX = _w.first.first;
//		real_t const fstY = _w.first.second;
//		real_t const sndX = _w.second.first;
//		real_t const sndY = _w.second.second;
//		real_t const width = 0.01;
//		real_t const height = 0.10;
//		real_t const length = sqrt((fstX - sndX) * (fstX - sndX) + (fstY - sndY) * (fstY - sndY));
//		real_t const direction = normAngle(atan2(sndY, sndX) - atan2(fstY, fstX));
//		real_t const angleSide = acos((length/3.0) / sqrt((length/3.0) * (length/3.0) + (width/2.0) * (width/2.0)));
//		real_t const norm = sqrt((length / 3.0) * (length / 3.0) + (width / 2.0) * (width / 2.0));
//		real_t const relativePosX = fstX - _agent->headPos.first;
//		real_t const relativePosY = fstY - _agent->headPos.second;
//
//		std::array<Coord3D_t, 6> tmpWallPosition;
//		tmpWallPosition[0] = {
//			relativePosX,
//			relativePosY,
//			0.0};
//
//		tmpWallPosition[1] = {
//			relativePosX - cos(direction) * length,
//			relativePosY - sin(direction) * length,
//			0.0};
//
//		tmpWallPosition[2] = {
//			relativePosX + cos(direction + pi() + angleSide) * norm,
//			relativePosY + sin(direction + pi() + angleSide) * norm,
//			0.0};
//
//		tmpWallPosition[3] = {
//			relativePosX + cos(direction + pi() - angleSide) * norm,
//			relativePosY + sin(direction + pi() - angleSide) * norm,
//			0.0};
//
//		tmpWallPosition[4] = {
//			relativePosX - cos(direction) * length / 3.0,
//			relativePosY - sin(direction) * length / 3.0,
//			height / 2.0};
//
//		tmpWallPosition[5] = {
//			relativePosX - cos(direction) * length / 3.0,
//			relativePosY - sin(direction) * length / 3.0,
//			-height / 2.0};
//
//		std::array<real_t, 6> angleRetina;
//		for(size_t k = 0; k < 6; ++k) {
//			real_t const positionXRetina = tmpWallPosition[k][0] * cos(_agent->direction) + tmpWallPosition[k][1] * sin(_agent->direction);
//			real_t const positionYRetina = -tmpWallPosition[k][0] * sin(_agent->direction) + tmpWallPosition[k][1] * cos(_agent->direction);
//			angleRetina[k] = atan2(positionYRetina, positionXRetina);
//		}
//
//		real_t const headRetinaX = relativePosX * cos(_agent->direction) + relativePosY * sin(_agent->direction);
//		real_t const headRetinaY = -relativePosX * sin(_agent->direction) + relativePosY * cos(_agent->direction);
//		real_t const headAngle = normAngle(atan2(headRetinaY, headRetinaX));
//		real_t const distAgent = sqrt(headRetinaX * headRetinaX + headRetinaY * headRetinaY);
//		real_t const agentCenter = normAngle(atan2(cos(angleRetina[0]) + cos(angleRetina[1]), sin(angleRetina[0]) + sin(angleRetina[1])) );
//
//		std::array<real_t, 6> agentSize = {
//			std::abs(angleRetina[0] - angleRetina[1]),
//			std::abs(angleRetina[0] - angleRetina[2]),
//			std::abs(angleRetina[0] - angleRetina[3]),
//			std::abs(angleRetina[3] - angleRetina[2]),
//			std::abs(angleRetina[2] - angleRetina[1]),
//			std::abs(angleRetina[3] - angleRetina[1]) };
//		size_t indiceLargestSize = 0;
//		real_t largestSize = 0.0;
//		for(size_t k = 0; k < 6; ++k) {
//			if(agentSize[k] > pi())
//				agentSize[k] = std::abs(2. * pi() - agentSize[k]);
//			if(agentSize[k] > largestSize) {
//				largestSize = agentSize[k];
//				indiceLargestSize = k;
//			}
//		}
//
//		std::array<Coord3D_t, 4> wallsCoords;
//		for(size_t k = 0; k < 4; ++k) {
//			wallsCoords[k] = tmpWallPosition[indicesCoords[indiceLargestSize][k]];
//		}
//		auto const areaWalls = surfaceSphereFish(wallsCoords);
//
//		// Find if the wall is within current agent FOV
//		if(		(areaWalls > _agent->agentMinAnglePerception) &&
//				((agentCenter < _agent->fov * pi() / 180.0) || (agentCenter > (360. - _agent->fov) * pi() / 180.0))) {
//			_totalAreaWalls += areaWalls;
//			// Compute PDF
//			for(size_t k = 0; k < evalAnglesNb; ++k) {
//				_wallsPDF[k] += areaWalls * _wallPDFBessel * exp(kappaWalls * cos(evaluedProb[k] - headAngle));
//			}
//		}
//	}
//
//	// Normalize PDF
//	if(_totalAreaWalls > 0.0) {
//		for(size_t k = 0; k < evalAnglesNb; ++k) {
//			_wallsPDF[k] /= _totalAreaWalls;
//		}
//	}

}


void ZonedBM::step() {
	_computeFishPDF();

	for(size_t run = 0; run < 50; ++run) {
		real_t newSpeed = _computeAgentSpeed();
		_detectZonesAroundAgent(newSpeed * _simulation.dt);
		_computeZonesPDF();
		if(followWalls)
			_computeWallsPDF();

		// Compute final PDF and CDF
		real_t cdfsum = 0.0;
		for(size_t k = 0; k < evalAnglesNb; ++k) {
			if(_zonesSensors[k] == 0) {
				// Not in arena
				_PDF[k] = 0.;
				_CDF[k] = cdfsum;
			} else {
				//real_t const val = (_neutPDFCenter[k] + alphasCenter * _totalAreaFishes * _fishPDF[k]) / (1. + alphasCenter * _totalAreaFishes);
				// XXX : maybe biaised toward zones (no agent cohesion .. ?)
				//real_t const val = (_neutPDFCenter[k] + alphasCenter * _totalAreaFishes * _fishPDF[k] + gammaZone * _zonesPDF[k]) / (1. + alphasCenter * _totalAreaFishes + gammaZone);

				real_t const val = (_neutPDFCenter[k] + alphasCenter * _totalAreaFishes * _fishPDF[k]) * (_zonesPDF[k]) / (1. + alphasCenter * _totalAreaFishes);
				//real_t const val = (_neutPDFCenter[k] + alphasCenter * _totalAreaFishes * _fishPDF[k]) * (gammaZone * _zonesPDF[k]) / (1. + alphasCenter * _totalAreaFishes);
				//std::cout << "DEB: val=" << val << " "
				//	<< " neut=" << (_neutPDFCenter[k] / (1. + alphasCenter * _totalAreaFishes))
				//	<< " agents=" << (alphasCenter * _totalAreaFishes * _fishPDF[k] / (1. + alphasCenter * _totalAreaFishes))
				//	<< " zones=" << (gammaZone * _zonesPDF[k])
				//	<< " _totalAreaFishes=" << _totalAreaFishes
				//	<< std::endl;



				//real_t const val = (_neutPDFCenter[k] + alphasCenter * _totalAreaFishes * _fishPDF[k] + beta * _totalAreaWalls * _wallsPDF[k]) * (gammaZone * _zonesPDF[k]) / (1. + alphasCenter * _totalAreaFishes + beta * _totalAreaFishes);
				//real_t const val = (_neutPDFCenter[k] + alphasCenter * 1. * _totalAreaFishes * _fishPDF[k] + beta * 1. * _totalAreaWalls * _wallsPDF[k] + gammaZone * 1. * _zonesPDF[k]) / (1. + alphasCenter * _totalAreaFishes + beta * _totalAreaFishes + gammaZone);
				//std::cout << "DEB: val=" << val << " "
				//	<< " neut=" << (_neutPDFCenter[k] / (1. + alphasCenter * _totalAreaFishes + beta * _totalAreaWalls + gammaZone))
				//	<< " agents=" << (alphasCenter * 1. * _totalAreaFishes * _fishPDF[k] / (1. + alphasCenter * _totalAreaFishes + beta * _totalAreaWalls + gammaZone))
				//	<< " walls=" << (beta * 1. * _totalAreaWalls * _wallsPDF[k] / (1. + alphasCenter * _totalAreaFishes + beta * _totalAreaWalls + gammaZone))
				//	<< " zones=" << (gammaZone * 1. * _zonesPDF[k] / (1. + alphasCenter * _totalAreaFishes + beta * _totalAreaWalls + gammaZone))
				//	<< " _totalAreaFishes=" << _totalAreaFishes
				//	<< " _totalAreaWalls=" << _totalAreaWalls
				//	<< std::endl;
				_PDF[k] = val;
				_CDF[k] = cdfsum + val;
				cdfsum += val;
			}
		}

		// Normalize CDF
		normalizeVec(_CDF, cdfsum);

		real_t z = unif01(rne());
		int dir = -1;
		for(int k = 0; k < _CDF.size(); ++k) {
			if(z <= _CDF[k]) {
				dir = k;
				break;
			}
		}

		real_t newDirection = _agent->direction;
		real_t deltaDirection = 0.0;
		if(dir >= 0) {
			//newDirection += evaluedProb[dir];
			deltaDirection = evaluedProb[dir];
		}

		// XXX Put in agent's code
		// Adjust newDirection according to agent's maxTurningRate
		if(deltaDirection < pi()) {
			if(deltaDirection > _agent->maxTurningRate) {
				deltaDirection = _agent->maxTurningRate;
			}
		} else {
			if(deltaDirection <= 2. * pi() - _agent->maxTurningRate) {
				deltaDirection = 2. * pi() - _agent->maxTurningRate;
			}
		}

		newDirection = normAngle(newDirection + deltaDirection);

		Coord_t newAgentsHeadPos = {
			_agent->headPos.first + _simulation.dt * newSpeed * cos(newDirection),
			_agent->headPos.second + _simulation.dt * newSpeed * sin(newDirection)
		};

		if(_simulation.arena.isInArena(newAgentsHeadPos)) {
			_agent->direction = newDirection;
			_agent->speed = newSpeed;
			_agent->headPos = newAgentsHeadPos;
			break;
		}
	}

	_agent->updateAgentPosition(_simulation.dt);



//	// Compute final PDF and CDF
//	real_t cdfsum = 0.0;
//		for(size_t k = 0; k < evalAnglesNb; ++k) {
//			real_t const val = (_neutPDFCenter[k] + alphasCenter * _totalAreaFishes * _fishPDF[k]) / (1. + alphasCenter * _totalAreaFishes);
//			_PDF[k] = val;
//			_CDF[k] = cdfsum + val;
//			cdfsum += val;
//		}
//
//	// Normalize CDF
//	normalizeVec(_CDF, cdfsum);
//
//	for(size_t run = 0; run < 50; ++run) {
//	//for(;;) {
//		real_t newSpeed = _computeAgentSpeed();
//		real_t z = unif01(rne());
//		int dir = -1;
//		for(int k = 0; k < _CDF.size(); ++k) {
//			if(z <= _CDF[k]) {
//				dir = k;
//				break;
//			}
//		}
//
//		real_t newDirection = _agent->direction;
//		real_t deltaDirection = 0.0;
//		if(dir >= 0) {
//			//newDirection += evaluedProb[dir];
//			deltaDirection = evaluedProb[dir];
//		}
//
//		// XXX Put in agent's code
//		// Adjust newDirection according to agent's maxTurningRate
//		if(deltaDirection < pi()) {
//			if(deltaDirection > _agent->maxTurningRate) {
//				deltaDirection = _agent->maxTurningRate;
//			}
//		} else {
//			if(deltaDirection <= 2. * pi() - _agent->maxTurningRate) {
//				deltaDirection = 2. * pi() - _agent->maxTurningRate;
//			}
//		}
//
//		newDirection = normAngle(newDirection + deltaDirection);
//
//		Coord_t newAgentsHeadPos = {
//			_agent->headPos.first + _simulation.dt * newSpeed * cos(newDirection),
//			_agent->headPos.second + _simulation.dt * newSpeed * sin(newDirection)
//		};
//
//		if(_simulation.arena.isInArena(newAgentsHeadPos)) {
//			_agent->direction = newDirection;
//			_agent->speed = newSpeed;
//			_agent->headPos = newAgentsHeadPos;
//			break;
//		} else {
//			//if(dir >= 0) {
//			//	_CDF[dir] = 0.0;
//			//}
//		}
//	}
}




// MODELINE	"{{{1
// vim:noexpandtab:softtabstop=4:shiftwidth=4:fileencoding=utf-8
// vim:foldmethod=marker
