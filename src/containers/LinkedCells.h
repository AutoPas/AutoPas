/*
 * LinkedCells.h
 *
 *  Created on: 17 Jan 2018
 *      Author: tchipevn
 */

#ifndef SRC_CONTAINERS_LINKEDCELLS_H_
#define SRC_CONTAINERS_LINKEDCELLS_H_

#include "ParticleContainer.h"
#include "CellBlock3D.h"
#include "utils/inBox.h"
#include "containers/cellPairTraversals/SlicedTraversal.h"

namespace autopas {

template<class Particle, class ParticleCell>
class LinkedCells : public ParticleContainer<Particle, ParticleCell> {
public:
	LinkedCells(const std::array<double, 3> boxMin, const std::array<double, 3> boxMax, double cutoff) :
		ParticleContainer<Particle, ParticleCell>(boxMin, boxMax, cutoff),
		_cellBlock(this->_data, boxMin, boxMax, cutoff) {
	}

	void addParticle(Particle& p) override {
		bool inBox = autopas::inBox(p.getR(), this->getBoxMin(), this->getBoxMax());
		if(inBox) {
			ParticleCell& cell =_cellBlock.getContainingCell(p.getR());
			cell.addParticle(p);
		} else {
			// todo
		}
	}

	void iteratePairwise(Functor<Particle>* f) override {

	}

	template<class ParticleFunctor>
	void iteratePairwise2(ParticleFunctor * f) {
		CellFunctor<Particle, ParticleCell, ParticleFunctor> cellFunctor(f);
//		cellFunctor.processCellAoSN3(this->_data[13]);
		SlicedTraversal<ParticleCell, CellFunctor<Particle, ParticleCell, ParticleFunctor>> traversal(this->_data, _cellBlock.getCellsPerDimensionWithHalo(), &cellFunctor);
		traversal.traverseCellPairs();

	}

private:
	CellBlock3D<ParticleCell> _cellBlock;
	// ThreeDimensionalCellHandler
};

} /* namespace autopas */

#endif /* SRC_CONTAINERS_LINKEDCELLS_H_ */
