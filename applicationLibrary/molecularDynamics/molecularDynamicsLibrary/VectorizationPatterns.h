#pragma once

  /**
   * Choice of the vectorization pattern
  */
  enum class VectorizationPattern { p1xVec, p2xVecDiv2, pVecDiv2x2, pVecx1, pVecxVec };