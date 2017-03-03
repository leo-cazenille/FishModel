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

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;


namespace Fishmodel {

py::array_t<real_t> runZonedModel(std::string mapFilename, std::uint32_t nbSteps = 10000, real_t dt = 1./3., std::uint32_t nbFishes = 5, std::uint32_t nbZones = 9, std::uint32_t seed = 0,
		py::array_t<real_t> kappa0 = py::array_t<real_t>(),
		py::array_t<real_t> kappaf = py::array_t<real_t>(),
		py::array_t<real_t> kappaw = py::array_t<real_t>(),
		py::array_t<real_t> alpha = py::array_t<real_t>(),
		py::array_t<real_t> beta = py::array_t<real_t>(),
		py::array_t<real_t> gamma = py::array_t<real_t>(),
		py::array_t<real_t> gammaz = py::array_t<real_t>(),
		py::array_t<real_t> speedHistogram = py::array_t<real_t>(),
		real_t minSpeed = 0.0, real_t maxSpeed = 0.03,
		py::array_t<real_t> wallsCoord = py::array_t<real_t>(),
		py::array_t<real_t> wallsDirectionCoord = py::array_t<real_t>(),
		py::array_t<bool> followWalls = py::array_t<bool>()
		) {
	size_t const nbAgents = nbFishes;

	py::buffer_info infoKappa0 = kappa0.request();
	py::buffer_info infoKappaf = kappaf.request();
	py::buffer_info infoKappaw = kappaw.request();
	py::buffer_info infoAlpha = alpha.request();
	py::buffer_info infoBeta = beta.request();
	py::buffer_info infoGamma = gamma.request();
	py::buffer_info infoGammaz = gammaz.request();
	py::buffer_info infoSpeedHistogram = speedHistogram.request();
	py::buffer_info infoWallsCoord = wallsCoord.request();
	py::buffer_info infoWallsDirectionCoord = wallsDirectionCoord.request();
	py::buffer_info infoFollowWalls = followWalls.request();

	// TODO verify shape of input array_t

	std::vector<real_t> ret(nbSteps * (1 + nbAgents * 3));
	std::vector<size_t> strides = {sizeof(real_t) * (1+nbAgents*3), sizeof(real_t)};
	std::vector<size_t> shape = {nbSteps, 1 + nbAgents * 3};
	size_t ndim = 2;

	// Init Random Number Generator
	rne().seed(seed);

	// Launch model simulation
	if(mapFilename == "")
		throw std::runtime_error("Please specify a mapFilename");
	Arena arena(mapFilename);

	SimulationFactory factory(arena);
	factory.nbFishes = nbFishes;
	factory.nbRobots = 0;
	factory.nbVirtuals = 0;
	factory.nbZones = nbZones;

	factory.behaviorFishes = "ZoneDependantBM";
	factory.behaviorRobots = "ZoneDependantBM";
	factory.behaviorVirtuals = "ZoneDependantBM";

	auto sim = factory.create();
	sim->dt = dt;


	std::vector<std::pair<Coord_t,Coord_t>> wallsCoord_;
	if(infoWallsCoord.shape[0] != 0) {
		for(size_t i = 0; i < infoWallsCoord.shape[0]; ++i) {
			real_t const coord0X = reinterpret_cast<real_t*>(infoWallsCoord.ptr)[i*infoWallsCoord.shape[1]*infoWallsCoord.shape[2] + 0*infoWallsCoord.shape[2] + 0];
			real_t const coord0Y = reinterpret_cast<real_t*>(infoWallsCoord.ptr)[i*infoWallsCoord.shape[1]*infoWallsCoord.shape[2] + 0*infoWallsCoord.shape[2] + 1];
			real_t const coord1X = reinterpret_cast<real_t*>(infoWallsCoord.ptr)[i*infoWallsCoord.shape[1]*infoWallsCoord.shape[2] + 1*infoWallsCoord.shape[2] + 0];
			real_t const coord1Y = reinterpret_cast<real_t*>(infoWallsCoord.ptr)[i*infoWallsCoord.shape[1]*infoWallsCoord.shape[2] + 1*infoWallsCoord.shape[2] + 1];
			std::pair<Coord_t, Coord_t> wall = {{coord0X, coord0Y}, {coord1X, coord1Y}};
			wallsCoord_.push_back(wall);
		}
	}

	std::vector<std::pair<Coord_t,Coord_t>> wallsDirectionCoord_;
	if(infoWallsDirectionCoord.shape[0] != 0) {
		for(size_t i = 0; i < infoWallsDirectionCoord.shape[0]; ++i) {
			real_t const coord0X = reinterpret_cast<real_t*>(infoWallsCoord.ptr)[i*infoWallsCoord.shape[1]*infoWallsCoord.shape[2] + 0*infoWallsCoord.shape[2] + 0];
			real_t const coord0Y = reinterpret_cast<real_t*>(infoWallsCoord.ptr)[i*infoWallsCoord.shape[1]*infoWallsCoord.shape[2] + 0*infoWallsCoord.shape[2] + 1];
			real_t const coord1X = reinterpret_cast<real_t*>(infoWallsCoord.ptr)[i*infoWallsCoord.shape[1]*infoWallsCoord.shape[2] + 1*infoWallsCoord.shape[2] + 0];
			real_t const coord1Y = reinterpret_cast<real_t*>(infoWallsCoord.ptr)[i*infoWallsCoord.shape[1]*infoWallsCoord.shape[2] + 1*infoWallsCoord.shape[2] + 1];
			std::pair<Coord_t, Coord_t> wall = {{coord0X, coord0Y}, {coord1X, coord1Y}};
			wallsDirectionCoord_.push_back(wall);
		}
	}


	// XXX Verify that shape[0] == zdb->nbZones()
	// Set zones parameters
	for(size_t j = 0; j < nbFishes; ++j) {
		ZoneDependantBehavior* zdb = reinterpret_cast<ZoneDependantBehavior*>(sim->fishes[j].second);
		for(size_t i = 1; i < zdb->nbZones() - 1; ++i) {
			ZonedBM* bm = reinterpret_cast<ZonedBM*>(zdb->behavior(i));
			if(infoKappa0.shape[0] != 0)
				bm->kappaNeutCenter = reinterpret_cast<real_t*>(infoKappa0.ptr)[i];
			if(infoKappaf.shape[0] != 0)
				bm->kappaFishes = reinterpret_cast<real_t*>(infoKappaf.ptr)[i];
			if(infoKappaw.shape[0] != 0)
				bm->kappaWalls = reinterpret_cast<real_t*>(infoKappaw.ptr)[i];
			if(infoAlpha.shape[0] != 0)
				bm->alphasCenter = reinterpret_cast<real_t*>(infoAlpha.ptr)[i];
			if(infoBeta.shape[0] != 0)
				bm->beta = reinterpret_cast<real_t*>(infoBeta.ptr)[i];
			if(infoGamma.shape[0] != 0)
				bm->gammaZone = reinterpret_cast<real_t*>(infoGamma.ptr)[i];
			//std::cout << "DEBUGzp1: " << bm->kappaNeutCenter << " " << bm->kappaFishes << " " << bm->alphasCenter << " " << bm->gammaZone << std::endl;
			if(infoGammaz.shape[0] != 0) {
				std::vector<real_t> affinity(infoGammaz.shape[1]);
				for(size_t k = 0; k < infoGammaz.shape[1]; ++k) {
					affinity[k] = reinterpret_cast<real_t*>(infoGammaz.ptr)[i*infoGammaz.shape[1]+k];
				}
				bm->zonesAffinity(affinity);
			}
			if(infoSpeedHistogram.shape[0] != 0) {
				bm->minSpeed = minSpeed;
				bm->maxSpeed = maxSpeed;
				bm->speedHistogram.clear();
				for(size_t k = 0; k < infoSpeedHistogram.shape[1]; ++k) {
					bm->speedHistogram.push_back(reinterpret_cast<real_t*>(infoSpeedHistogram.ptr)[i*infoSpeedHistogram.shape[1]+k]);
				}
			}
			if(infoWallsCoord.shape[0] != 0) {
				bm->wallsCoord = wallsCoord_;
			}
			if(infoWallsDirectionCoord.shape[0] != 0) {
				bm->wallsDirectionCoord = wallsDirectionCoord_;
			}
			bm->followWalls = reinterpret_cast<bool*>(infoFollowWalls.ptr)[i];
			bm->reinit();
			//std::cout << "DEBUGaa: " << j << "," << i << ": " << bm->kappaFishes << " " << bm->kappaNeutCenter << std::endl;
		}
	}


	//std::cout << "DEBUGzp2" << std::endl;
	//std::cout << "DEBUGzp3: " << reinterpret_cast<real_t*>(infoKappa0.ptr)[1] << " " << reinterpret_cast<real_t*>(infoKappaf.ptr)[1] << " " << reinterpret_cast<real_t*>(infoAlpha.ptr)[1] << " " << reinterpret_cast<real_t*>(infoGamma.ptr)[1];
	//for(size_t k = 0; k < infoGammaz.shape[1]; ++k) {
	//	std::cout << " " << reinterpret_cast<real_t*>(infoGammaz.ptr)[1*infoGammaz.shape[1]+k];
	//}
	//for(size_t k = 0; k < infoSpeedHistogram.shape[1]; ++k) {
	//	std::cout << " " << reinterpret_cast<real_t*>(infoSpeedHistogram.ptr)[1*infoSpeedHistogram.shape[1]+k];
	//}
	//std::cout << std::endl;


	std::cout << std::fixed; // << std::setprecision(3);
	for(size_t i = 0; i < nbSteps; ++i) {
		real_t t = static_cast<real_t>(i) * sim->dt; // / 15.0; // XXX
		//std::cout << t << "\t";
		//std::cout << sim->printCurrentPositions() << std::endl;

		ret[0+i*(1+nbAgents*3)] = t;
		for(size_t j = 0; j < nbFishes; ++j) {
			ret[1+i*(1+nbAgents*3)+j*3+0] = sim->fishes[j].first->headPos.first;
			ret[1+i*(1+nbAgents*3)+j*3+1] = sim->fishes[j].first->headPos.second;
			ret[1+i*(1+nbAgents*3)+j*3+2] = sim->fishes[j].first->direction;
		}

		sim->step();
	}

	//std::cout << "DEBUG1: " << kappa0 << std::endl;

	return py::array(py::buffer_info(ret.data(), sizeof(real_t),
				py::format_descriptor<real_t>::value,
				ndim, shape, strides));
}

}


PYBIND11_PLUGIN(model) {
	using namespace pybind11::literals;
	using namespace Fishmodel;
	py::module m("model", "Fish model");

	m.def("runZonedModel", &Fishmodel::runZonedModel, "run model with zones (bm)",
			"mapFilename"_a,
			"nbSteps"_a = 10000,
			"dt"_a = 1./3.,
			"nbFishes"_a = 5,
			"nbZones"_a = 9,
			"seed"_a = 0,
			"kappa0"_a = py::array_t<real_t>(),
			"kappaf"_a = py::array_t<real_t>(),
			"kappaw"_a = py::array_t<real_t>(),
			"alpha"_a = py::array_t<real_t>(),
			"beta"_a = py::array_t<real_t>(),
			"gamma"_a = py::array_t<real_t>(),
			"gammaz"_a = py::array_t<real_t>(),
			"speedHistogram"_a = py::array_t<real_t>(),
			"minSpeed"_a = 0.0,
			"maxSpeed"_a = 0.03,
			"wallsCoord"_a = py::array_t<real_t>(),
			"wallsDirectionCoord"_a = py::array_t<real_t>(),
			"followWalls"_a = py::array_t<bool>()
	);

	return m.ptr();
}


// MODELINE	"{{{1
// vim:noexpandtab:softtabstop=4:shiftwidth=4:fileencoding=utf-8
// vim:foldmethod=marker
