/**
 * @file
 * @author Leo Cazenille <leo.cazenille@gmail.com>
 *
 *
 */

#ifndef BMWITHWALLS_H
#define BMWITHWALLS_H

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


//! BM implementation for Summer School 2017.
//!  There are two zones: center and close to the walls, defined by the parameter $d$ (wallDistanceThreshold).
//!  To reduce the number of parameters, $\kappa_f$, $\alpha_0$ and $\kappa_0$ are the same for both zones.
//!  For simplicity sake, agents speeds are drawn from a log-normal distribution.
struct BMWithWalls : public Behavior {
public:
	real_t kappaFishes = 10.0; //20.0        //! \kappa_f
	real_t alpha = 25.0; //55.0;             //! \alpha_0
	real_t kappaNeutCenter = 6.3;            //! \kappa_0
	real_t wallDistanceThreshold = 0.02;     //! d
	real_t repulsionFromAgentsAtDist = 0.02; //! Repulsion from agent if other is too close

	real_t currentDistanceToNearestWall = 1.0;

	std::vector<std::pair<Coord_t,Coord_t>> wallsCoord;

protected:
	std::vector<real_t> _neutPDFCenter;  //! f_o(\theta)
	std::vector<real_t> _neutPDFWalls;
	std::vector<real_t> _fishPDF;        //! f_F(\theta)
	std::vector<real_t> _PDF;            //! f(\theta)
	std::vector<real_t> _CDF;            //! CDF of f(\theta)

protected:
	real_t _fishPDFBessel;
	real_t _totalAreaFishes = 0.0;  //! \A_{T_f}  (automatically set)


protected:
	virtual void _computeFishPDF();
	virtual void _computeWallsPDF();
	virtual real_t _computeAgentSpeed();


public:
	BMWithWalls(Simulation& simulation, Agent* agent = nullptr) :
			Behavior(simulation, agent),
			_neutPDFCenter(evalAnglesNb),
			_neutPDFWalls(evalAnglesNb),
			_fishPDF(evalAnglesNb),
			_PDF(evalAnglesNb),
			_CDF(evalAnglesNb) {
		reinit();
	}


public:
	virtual void reinit();
	virtual void step();
};



}
#endif

// MODELINE	"{{{1
// vim:noexpandtab:softtabstop=4:shiftwidth=4:fileencoding=utf-8
// vim:foldmethod=marker
