/*
 * SlicedTraversal.h
 *
 *  Created on: 22 Jan 2018
 *      Author: tchipevn
 */

#ifndef SRC_CONTAINERS_CELLPAIRTRAVERSALS_SLICEDTRAVERSAL_H_
#define SRC_CONTAINERS_CELLPAIRTRAVERSALS_SLICEDTRAVERSAL_H_

#include "CellPairTraversal.h"
#include "utils/ThreeDimensionalMapping.h"

namespace autopas {

template <class ParticleCell, class CellFunctor>
class SlicedTraversal: public CellPairTraversals<ParticleCell, CellFunctor> {
public:
	explicit SlicedTraversal(
		std::vector<ParticleCell>& cells,
		const std::array<unsigned long, 3>& dims,
		CellFunctor * cellfunctor): CellPairTraversals<ParticleCell, CellFunctor>(cells, dims, cellfunctor) {
		computeOffsets();
	}

	void traverseCellPairs() override;

private:
	void processBaseCell(unsigned long baseIndex) const;
	void computeOffsets();

	void processCell(unsigned long cellIndex) const;
	void processCellPair(unsigned long cellIndex1, unsigned long cellIndex2) const;

	std::array<std::pair<unsigned long, unsigned long>, 14> _cellPairOffsets;
	std::array<unsigned long, 8> _cellOffsets;
};

template<class ParticleCell, class CellFunctor>
inline void SlicedTraversal<ParticleCell, CellFunctor>::processBaseCell(
		unsigned long baseIndex) const {

	using std::pair;

	const int num_pairs = _cellPairOffsets.size();
	for(int j = 0; j < num_pairs; ++j) {
		pair<long, long> current_pair = _cellPairOffsets[j];

		unsigned offset1 = current_pair.first;
		unsigned cellIndex1 = baseIndex + offset1;

		unsigned offset2 = current_pair.second;
		unsigned cellIndex2 = baseIndex + offset2;

		ParticleCell& cell1 = this->_cells->at(cellIndex1);
		ParticleCell& cell2 = this->_cells->at(cellIndex2);

		if(cellIndex1 == cellIndex2) {
			this->_cellFunctor->processCell(cell1);
		}
		else {
			this->_cellFunctor->processCellPair(cell1, cell2);
		}
	}
}

template<class ParticleCell, class CellFunctor>
inline void SlicedTraversal<ParticleCell, CellFunctor>::computeOffsets() {
	using ThreeDimensionalMapping::threeToOneD;
	using std::make_pair;

	unsigned long int o   = threeToOneD(0ul, 0ul, 0ul, this->_dims); // origin
	unsigned long int x   = threeToOneD(1ul, 0ul, 0ul, this->_dims); // displacement to the right
	unsigned long int y   = threeToOneD(0ul, 1ul, 0ul, this->_dims); // displacement ...
	unsigned long int z   = threeToOneD(0ul, 0ul, 1ul, this->_dims);
	unsigned long int xy  = threeToOneD(1ul, 1ul, 0ul, this->_dims);
	unsigned long int yz  = threeToOneD(0ul, 1ul, 1ul, this->_dims);
	unsigned long int xz  = threeToOneD(1ul, 0ul, 1ul, this->_dims);
	unsigned long int xyz = threeToOneD(1ul, 1ul, 1ul, this->_dims);

	int i = 0;
	// if incrementing along X, the following order will be more cache-efficient:
	_cellPairOffsets[i++] = make_pair(o, o  );
	_cellPairOffsets[i++] = make_pair(o, y  );
	_cellPairOffsets[i++] = make_pair(y, z  );
	_cellPairOffsets[i++] = make_pair(o, z  );
	_cellPairOffsets[i++] = make_pair(o, yz );

	_cellPairOffsets[i++] = make_pair(x, yz );
	_cellPairOffsets[i++] = make_pair(x, y  );
	_cellPairOffsets[i++] = make_pair(x, z  );
	_cellPairOffsets[i++] = make_pair(o, x  );
	_cellPairOffsets[i++] = make_pair(o, xy );
	_cellPairOffsets[i++] = make_pair(xy, z );
	_cellPairOffsets[i++] = make_pair(y, xz );
	_cellPairOffsets[i++] = make_pair(o, xz );
	_cellPairOffsets[i++] = make_pair(o, xyz);

	i = 0;
	_cellOffsets[i++] =   o;
	_cellOffsets[i++] =   y;
	_cellOffsets[i++] =   z;
	_cellOffsets[i++] =  yz;

	_cellOffsets[i++] =   x;
	_cellOffsets[i++] =  xy;
	_cellOffsets[i++] =  xz;
	_cellOffsets[i++] = xyz;
}


template<class ParticleCell, class CellFunctor>
inline void SlicedTraversal<ParticleCell, CellFunctor>::traverseCellPairs() {
	using std::array;
	array<unsigned long, 3> end;
	for (int d = 0; d < 3; ++d) {
		end[d] = this->_dims[d] - 1;
	}

	for (unsigned long z = 0; z < end[2]; ++z) {
		for (unsigned long y = 0; y < end[1]; ++y) {
			for (unsigned long x = 0; x < end[0]; ++x) {
				unsigned long ind = ThreeDimensionalMapping::threeToOneD(x, y, z, this->_dims);
				processBaseCell(ind);
			}
		}
	}
}

} /* namespace autopas */


#endif /* SRC_CONTAINERS_CELLPAIRTRAVERSALS_SLICEDTRAVERSAL_H_ */
