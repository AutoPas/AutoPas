
/**
 * @file
 * @brief VTK File Writer
 */
#include "VTKWriter.h"

#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

using namespace std;

namespace outputWriter {

VTKWriter::VTKWriter() {
  // TODO Auto-generated constructor stub
}

VTKWriter::~VTKWriter() {
  // TODO Auto-generated destructor stub
}

void VTKWriter::initializeOutput(int numParticles) {
  vtkFile = new VTKFile_t("UnstructuredGrid");

  // per point, we add type, position, velocity and force
  PointData pointData;
  DataArray_t id(type::Float32, "id", 1);
  DataArray_t velocity(type::Float32, "velocity", 3);
  DataArray_t forces(type::Float32, "force", 3);
  DataArray_t typeId(type::Int32, "typeId", 1);
  pointData.DataArray().push_back(id);
  pointData.DataArray().push_back(velocity);
  pointData.DataArray().push_back(forces);
  pointData.DataArray().push_back(typeId);

  CellData cellData;  // we don't have cell data => leave it empty

  // 3 coordinates
  Points points;
  DataArray_t pointCoordinates(type::Float32, "points", 3);
  points.DataArray().push_back(pointCoordinates);

  Cells cells;  // we don't have cells, => leave it empty
  // for some reasons, we have to add a dummy entry for paraview
  DataArray_t cells_data(type::Float32, "types", 0);
  cells.DataArray().push_back(cells_data);

  PieceUnstructuredGrid_t piece(pointData, cellData, points, cells, numParticles, 0);
  UnstructuredGrid_t unstructuredGrid(piece);
  vtkFile->UnstructuredGrid(unstructuredGrid);
}

void VTKWriter::writeFile(const std::string &filename, int iteration) {
  stringstream strstr;
  strstr << filename << setfill('0') << setw(4) << iteration << ".vtu";

  std::ofstream file(strstr.str().c_str());
  VTKFile(file, *vtkFile);
  delete vtkFile;
}

void VTKWriter::plotParticle(autopas::MoleculeLJ<> &p) {
  if (vtkFile->UnstructuredGrid().present()) {
    //    cout << "UnstructuredGrid is present" << endl;
  } else {
    //    cout << "ERROR: No UnstructuredGrid present" << endl;
  }

  PointData::DataArray_sequence &pointDataSequence = vtkFile->UnstructuredGrid()->Piece().PointData().DataArray();
  PointData::DataArray_iterator dataIterator = pointDataSequence.begin();

  dataIterator->push_back(p.getID());
  // cout << "Appended id data in: " << dataIterator->Name();

  dataIterator++;
  dataIterator->push_back(p.getV()[0]);
  dataIterator->push_back(p.getV()[1]);
  dataIterator->push_back(p.getV()[2]);
  // cout << "Appended velocity data in: " << dataIterator->Name();

  dataIterator++;
  dataIterator->push_back(p.getF()[0]);
  dataIterator->push_back(p.getF()[1]);
  dataIterator->push_back(p.getF()[2]);
  // cout << "Appended force data in: " << dataIterator->Name();

  dataIterator++;
  dataIterator->push_back(p.getTypeId());

  Points::DataArray_sequence &pointsSequence = vtkFile->UnstructuredGrid()->Piece().Points().DataArray();
  Points::DataArray_iterator pointsIterator = pointsSequence.begin();
  pointsIterator->push_back(p.getR()[0]);
  pointsIterator->push_back(p.getR()[1]);
  pointsIterator->push_back(p.getR()[2]);
}

}  // namespace outputWriter
