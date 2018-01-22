/*
 * CellBlock3D.h
 *
 *  Created on: 19 Jan 2018
 *      Author: tchipevn
 */

#ifndef SRC_CONTAINERS_CELLBLOCK3D_H_
#define SRC_CONTAINERS_CELLBLOCK3D_H_

#include "utils/ThreeDimensionalMapping.h"
#include <array>
#include <vector>
#include <cmath>

namespace autopas {

template<class ParticleCell>
class CellBlock3D {
public:
	typedef std::size_t index_t;

	CellBlock3D(std::vector<ParticleCell>& vec, const std::array<double, 3> bMin,
			const std::array<double, 3> bMax, double interactionLength) {
		rebuild(vec, bMin, bMax, interactionLength);
	}

	ParticleCell & getCell(index_t index1d) const ;
	ParticleCell & getCell(const std::array<index_t, 3>& index3d) const ;

	void rebuild(std::vector<ParticleCell>& vec, const std::array<double, 3> bMin,
			const std::array<double, 3> bMax, double interactionLength);

	// this class doesn't actually know about particles
	ParticleCell & getContainingCell(const std::array<double, 3>& pos) const;

	std::array<index_t, 3> get3DIndexOfPosition(const std::array<double, 3>& pos) const;
	index_t get1DIndexOfPosition(const std::array<double, 3>& pos) const;

	const std::array<index_t, 3>& getCellsPerDimensionWithHalo() const {
		return _cellsPerDimensionWithHalo;
	}

private:
	std::array<index_t ,3> index3D(index_t index1d) const;
	index_t index1D(const std::array<index_t ,3>& index3d) const;

	std::array<index_t, 3> _cellsPerDimensionWithHalo;
	index_t _numCells;
	std::vector<ParticleCell>* _vec1D;

	std::array<double, 3> _boxMin, _boxMax;
	std::array<double, 3> _haloBoxMin, _haloBoxMax;

	double _interactionLength;
//	int _cellsPerInteractionLength; = 1 hardcode to 1 for now, because flag manager will also need to be adapted

	std::array<double, 3> _cellLength;

	// 1 over above. Since this value is needed for sorting particles in cells, it is computed quite often
	std::array<double, 3> _cellLengthReciprocal;

//	CellBorderAndFlagManager _cellBorderAndFlagManager;
};

template<class ParticleCell>
inline typename CellBlock3D<ParticleCell>::index_t CellBlock3D<ParticleCell>::get1DIndexOfPosition(
		const std::array<double, 3>& pos) const {
	return index1D(get3DIndexOfPosition(pos));
}

template<class ParticleCell>
inline std::array<typename CellBlock3D<ParticleCell>::index_t, 3> CellBlock3D<ParticleCell>::get3DIndexOfPosition(
		const std::array<double, 3>& pos) const {
	std::array<typename CellBlock3D<ParticleCell>::index_t, 3> cellIndex;

	std::array<double, 3> localPoint = pos;
	for (int dim = 0; dim < 3; dim++) {
		long int value =(static_cast<long int>(floor((localPoint[dim] - _boxMin[dim]) * _cellLengthReciprocal[dim]))) + 1l;
		long int nonnegativeValue = std::max(value, 0l);
		index_t nonLargerValue = std::min(static_cast<index_t>(nonnegativeValue), _cellsPerDimensionWithHalo[dim] - 1);
		cellIndex[dim] = nonLargerValue;
	}

	return cellIndex;
	// in very rare cases rounding is stupid, thus we need a check...
	// todo when the border and flag manager is there
}

template<class ParticleCell>
inline ParticleCell& CellBlock3D<ParticleCell>::getContainingCell(
		const std::array<double, 3>& pos) const {
	auto ind = get1DIndexOfPosition(pos);
	return getCell(ind);
}

template<class ParticleCell>
inline void CellBlock3D<ParticleCell>::rebuild(std::vector<ParticleCell>& vec, const std::array<double, 3> bMin,
		const std::array<double, 3> bMax, double interactionLength) {
	_vec1D = &vec;
	_boxMin = bMin;
	_boxMax = bMax;
	_interactionLength = interactionLength;

	// compute cell length
	_numCells = 1;
	for (int d = 0; d < 3; ++d) {
		double diff = _boxMax[d] - _boxMin[d];
		index_t cellsPerDim = std::floor((diff) / _interactionLength);
		// at least one central cell
		cellsPerDim = std::max(cellsPerDim, 1ul);

		_cellsPerDimensionWithHalo[d] = cellsPerDim + 2;

		_cellLength[d] = diff / cellsPerDim;
		_cellLengthReciprocal[d] = cellsPerDim / diff; // compute with least rounding possible

		_haloBoxMin[d] = _boxMin[d] - _cellLength[d];
		_haloBoxMax[d] = _boxMax[d] + _cellLength[d];

		_numCells *= _cellsPerDimensionWithHalo[d];
	}



	_vec1D->resize(_numCells);

//	ParticleCell::_cellBorderAndFlagManager.init(
//		_cellsPerDimensionWithHalo,
//		_haloBoxMin, _haloBoxMax,
//		_boxMin, _boxMax,
//		_cellLength);
}

template<class ParticleCell>
inline ParticleCell& CellBlock3D<ParticleCell>::getCell(index_t index1d) const {
	return _vec1D->at(index1d);
}

template<class ParticleCell>
inline ParticleCell& CellBlock3D<ParticleCell>::getCell(const std::array<index_t, 3>& index3d) const {
	return _vec1D->at(index1D(index3d));
}

template<class ParticleCell>
inline std::array<typename CellBlock3D<ParticleCell>::index_t, 3> CellBlock3D<ParticleCell>::index3D(
		index_t index1d) const {
	return ThreeDimensionalMapping::oneToThreeD(index1d, _cellsPerDimensionWithHalo);
}

template<class ParticleCell>
inline typename CellBlock3D<ParticleCell>::index_t CellBlock3D<ParticleCell>::index1D(const std::array<index_t ,3>& index3d) const {
	return ThreeDimensionalMapping::threeToOneD(index3d, _cellsPerDimensionWithHalo);
}

} /* namespace autopas */


#endif /* SRC_CONTAINERS_CELLBLOCK3D_H_ */
