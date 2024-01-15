/**
 * @file HierarchicalGridsHelpers.h
 * @author F. Hoppe
 * @date 2023-12-07
 * 
 */

#pragma once


namespace autopas {

/**
 * @brief Class of helpers for Hierarchical grids
 * 
 */
class HierarchicalGridsHelpers {
  public:

    static void setNumberOfLevels(unsigned int levels) {

      _numberOfHGLevels = levels;

    }

    static unsigned int getNumberOfLevels() {

      return _numberOfHGLevels;

    }


  //private:
    /**
     * Contains the number of levels to be set up for the HierarchicalGrids
     */
    static unsigned int _numberOfHGLevels;

}; //HierarchicalGridsHelpers
 
} // Namespace autopas