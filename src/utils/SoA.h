/*
 * SoA.h
 *
 *  Created on: 17 Jan 2018
 *      Author: tchipevn
 */

#ifndef DEPENDENCIES_EXTERNAL_AUTOPAS_SRC_SOA_H_
#define DEPENDENCIES_EXTERNAL_AUTOPAS_SRC_SOA_H_

#include <vector>
#include <array>

namespace autopas {

class SoA {
public:
	void pushPos(const std::array<double, 3>& p) {
		_posX.push_back(p[0]);
		_posY.push_back(p[1]);
		_posZ.push_back(p[2]);
	}
	void pushFor(const std::array<double, 3>& f) {
		_forX.push_back(f[0]);
		_forY.push_back(f[1]);
		_forZ.push_back(f[2]);
	}
	std::array<double, 3> readPos(int i) const {
		return std::array<double, 3>({_posX[i], _posY[i], _posZ[i]});
	}
	std::array<double, 3> readFor(int i) const {
		return std::array<double, 3>({_forX[i], _forY[i], _forZ[i]});
	}
	int getNumParticles() const {
		return _posX.size();
	}
private:
	std::vector<double> _posX, _posY, _posZ;
	std::vector<double> _forX, _forY, _forZ;
};

} /* namespace autopas */

#endif /* DEPENDENCIES_EXTERNAL_AUTOPAS_SRC_SOA_H_ */
