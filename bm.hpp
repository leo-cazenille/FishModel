/**
 * @file
 * @author Leo Cazenille <leo.cazenille@gmail.com>
 *
 *
 */

#ifndef BM_H
#define BM_H

#include <utility>
#include <vector>
#include <string>
#include <memory>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>

#include <iostream>

#include "model.hpp"
#include "zones.hpp"


namespace Fishmodel {


//struct ModelBehavior : public Behavior {
//public:
//	real_t neighborsInfluence = 1.0;
//	real_t wallsInfluence = 1.0;
//
//	real_t inertia = 1.0;
//	real_t perceptionInfluence = 7.0; // XXX DEPRECATED
//
//	real_t kappaNeut = 6.3;
//	//real_t kappaNeut = 3.0;
//	real_t kappaFishes = 20.0;
//	real_t kappaWalls = 20.0;
//
//protected:
//	real_t _fishPDFBessel;
//
//protected:
//	std::vector<real_t> _neutPDF;
//	std::vector<real_t> _neutPDF2;
//	std::vector<real_t> _fishPDF;
//	std::vector<real_t> _wallPDF;
//	std::vector<real_t> _PDF;
//	std::vector<real_t> _CDF;
//	real_t _alphas;
//	real_t _betas;
//	real_t _gammas;
//
//protected:
//	std::vector<size_t> _retinaID; // TODO
//	std::vector<double> _retinaSpeed; // TODO
//
//protected:
//	virtual void _computeFishPDF();
//	virtual void _computeWallPDF();
//	virtual real_t _computeAgentSpeed();
//
//public:
//	virtual void reinit();
//	virtual void step();
//
//	ModelBehavior(Simulation& simulation, Agent* agent = nullptr) :
//		Behavior(simulation, agent),
//		_neutPDF(evalAnglesNb),
//		_neutPDF2(evalAnglesNb),
//		_fishPDF(evalAnglesNb),
//		_wallPDF(evalAnglesNb),
//		_PDF(evalAnglesNb),
//		_CDF(evalAnglesNb) {
//		reinit();
//	}
//};
//
//
//struct Model3DBehavior : public ModelBehavior {
//public:
//	real_t repulsionFromAgentsAtDist = 0.02; //0.05;
//	//real_t thresholdWalls = 0.05; //0.05; //0.20;
//
//	real_t alphasCenter = 55.0;
//	real_t alphasWalls = 20.0;
//	real_t kappaNeutWall = 20.0;
//	real_t kappaNeutCenter = 6.3; //4.0; //6.3;
//
//protected:
//	//real_t _wallPDFBessel;
//	real_t _totalAreaFishes = 0.0;
//
//	real_t _distToWall = 0.0;
//
//
//	std::vector<real_t> _neutPDFCenter;
//
//protected:
//	virtual void _computeFishPDF();
//	virtual void _computeDistToWall();
//
//public:
//	Model3DBehavior(Simulation& simulation, Agent* agent = nullptr) :
//		ModelBehavior(simulation, agent),
//		_neutPDFWall(evalAnglesNb),
//		_neutPDFCenter(evalAnglesNb) {
//		reinit();
//	}
//	virtual void reinit();
//	virtual void step();
//};



struct BM : public Behavior {
public:
	real_t kappaFishes = 20.0;  //! \kappa_f
	real_t alphasCenter = 55.0;   //! \alpha_0
	real_t kappaNeutCenter = 6.3; //! \kappa_0
	real_t repulsionFromAgentsAtDist = 0.02; //! Repulsion from agent if other is too close

protected:
	std::vector<real_t> _neutPDFCenter;  //! f_o(\theta)
	std::vector<real_t> _fishPDF;        //! f_F(\theta)
	std::vector<real_t> _PDF;            //! f(\theta)
	std::vector<real_t> _CDF;            //! CDF of f(\theta)

protected:
	real_t _fishPDFBessel;
	real_t _totalAreaFishes = 0.0;  //! \A_{T_f}  (automatically set)


protected:
	virtual void _computeFishPDF();
	virtual real_t _computeAgentSpeed();


public:
	BM(Simulation& simulation, Agent* agent = nullptr) :
			Behavior(simulation, agent),
			_neutPDFCenter(evalAnglesNb),
			_fishPDF(evalAnglesNb),
			_PDF(evalAnglesNb),
			_CDF(evalAnglesNb) {
		reinit();
	}

public:
	virtual void reinit();
	virtual void step();
};



struct ZonedBM : public BM {
public:
	real_t gammaZone = 55.0;

protected:
	ZoneDependantBehavior* _zdb;

	std::vector<real_t> _zonesCaptors;
	std::vector<real_t> _zonesAffinity;
	std::vector<real_t> _zonesPDF;

protected:
	virtual void _detectZonesAroundAgent(real_t r);
	virtual void _computeZonesPDF();
	virtual real_t _computeAgentSpeed();

public:
	ZonedBM(Simulation& simulation, Agent* agent = nullptr) :
			BM(simulation, agent) {
		reinit();
	}

	ZonedBM(Simulation& simulation, Agent* agent, ZoneDependantBehavior* zdb, std::vector<real_t> zonesAffinity) :
			BM(simulation, agent),
			_zdb(zdb),
			_zonesAffinity(zonesAffinity) {
		reinit();
	}

public:
	virtual void reinit();
	virtual void step();

	inline void zdb(decltype(_zdb) zdb) { _zdb = zdb; }
	inline void zonesAffinity(decltype(_zonesAffinity)& zonesAffinity) { _zonesAffinity = zonesAffinity; }
};



}
#endif

// MODELINE	"{{{1
// vim:noexpandtab:softtabstop=4:shiftwidth=4:fileencoding=utf-8
// vim:foldmethod=marker
