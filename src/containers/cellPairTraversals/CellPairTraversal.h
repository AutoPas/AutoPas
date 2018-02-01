/*
 * CellPairTraversal.h
 *
 *  Created on: 22 Jan 2018
 *      Author: tchipevn
 */

#ifndef SRC_CONTAINERS_CELLPAIRTRAVERSALS_CELLPAIRTRAVERSAL_H_
#define SRC_CONTAINERS_CELLPAIRTRAVERSALS_CELLPAIRTRAVERSAL_H_

#include <vector>
#include <array>

namespace autopas {

template <class ParticleCell, class CellFunctor>
class CellPairTraversals {
public:
	CellPairTraversals(
		std::vector<ParticleCell>& cells,
		const std::array<unsigned long, 3>& dims,
		CellFunctor * cellFunctor): _cells(&cells), _dims(dims), _cellFunctor(cellFunctor) {}

	virtual ~CellPairTraversals() = default;

	virtual void rebuild(std::vector<ParticleCell> &cells,
						 const std::array<unsigned long, 3> &dims) {
		_cells = &cells;
		_dims = dims;
	};

	virtual void traverseCellPairs() = 0;

protected:
	std::vector<ParticleCell> * _cells;
	std::array<unsigned long, 3> _dims;
	CellFunctor * _cellFunctor;
};

} /* namespace autopas */

#endif /* SRC_CONTAINERS_CELLPAIRTRAVERSALS_CELLPAIRTRAVERSAL_H_ */
