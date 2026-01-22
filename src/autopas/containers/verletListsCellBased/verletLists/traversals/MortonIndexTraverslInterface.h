#pragma once

namespace autopas {

class MortonIndexTraversalInterface {
 public:
  /**
   * Destructor
   */
  virtual ~MortonIndexTraversalInterface() = default;

  void setCellsByMortonIndex(const std::vector<size_t> &cellsByMortonIndex) {
    _cellsByMortonIndex = &cellsByMortonIndex;
  }

 protected:
  const std::vector<size_t> *_cellsByMortonIndex = nullptr;
};

}  // namespace autopas