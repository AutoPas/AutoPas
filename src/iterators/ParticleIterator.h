/*
 * ParticleIterator.h
 *
 *  Created on: 17 Jan 2018
 *      Author: tchipevn
 */

#ifndef DEPENDENCIES_EXTERNAL_AUTOPAS_SRC_PARTICLEITERATOR_H_
#define DEPENDENCIES_EXTERNAL_AUTOPAS_SRC_PARTICLEITERATOR_H_

#include "SingleCellIterator.h"
#include <vector>

namespace autopas {

template<class Particle, class ParticleCell>
class ParticleIterator {
public:
	ParticleIterator(std::vector<ParticleCell> * cont) : _vectorOfCells(cont), _iteratorAcrossCells(cont->begin()), _iteratorWithinOneCell() {
		if (_iteratorAcrossCells < cont->end()) {
			_iteratorWithinOneCell = SingleCellIterator<Particle, ParticleCell>(&(*_iteratorAcrossCells));
			if (not _iteratorWithinOneCell.isValid()) {
				next_non_empty_cell();
			}
		}
	}

	void operator ++ () {
		if (_iteratorWithinOneCell.isValid()) {
			++_iteratorWithinOneCell;
		}

		// don't merge into if-else, _cell_iterator may becoeme invalid after ++

		if (not _iteratorWithinOneCell.isValid()) {
			next_non_empty_cell();
		}
	}
	void next_non_empty_cell() {
		// find the next non-empty cell
		const int stride = 1; // num threads
		for (_iteratorAcrossCells += stride; _iteratorAcrossCells < _vectorOfCells->end(); _iteratorAcrossCells += stride) {

			const ParticleCell & c = *_iteratorAcrossCells;

			if(c.isNotEmpty()) {
				_iteratorWithinOneCell = SingleCellIterator<Particle, ParticleCell>(&(*_iteratorAcrossCells));
				break;
			}
		}
	}
	Particle& operator *  () const {
		return _iteratorWithinOneCell.operator*();
	}
	Particle* operator -> () const {
		return &(this->operator*());
	}
	void deleteCurrentParticle() {
//		cout << "deleteCurrentParticle is still ToDo" << endl;
	}

	bool isValid() {
		return _vectorOfCells != nullptr and _iteratorAcrossCells < _vectorOfCells->end() and _iteratorWithinOneCell.isValid();
	}

private:
	std::vector<ParticleCell>* _vectorOfCells;
	typename std::vector<ParticleCell>::iterator _iteratorAcrossCells;
	SingleCellIterator<Particle, ParticleCell> _iteratorWithinOneCell;
};

} /* namespace autopas */

#endif /* DEPENDENCIES_EXTERNAL_AUTOPAS_SRC_PARTICLEITERATOR_H_ */
