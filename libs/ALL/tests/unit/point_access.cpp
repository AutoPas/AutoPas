/*
Copyright 2018-2020 Rene Halver, Forschungszentrum Juelich GmbH, Germany
Copyright 2018-2020 Godehard Sutmann, Forschungszentrum Juelich GmbH, Germany

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or
   other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software without
   specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#define BOOST_TEST_MAIN

#define BOOST_TEST_MODULE point_access_and_creation
#include "ALL_Point.hpp"
#include <boost/test/unit_test.hpp>
#include <list>
#include <vector>

BOOST_AUTO_TEST_SUITE(point_access)

BOOST_AUTO_TEST_CASE(array) {
  int dimension = 3;
  double x = 2.5;
  double y = -5.7;
  double z = 8.3;

  double data[3];
  data[0] = x;
  data[1] = y;
  data[2] = z;

  ALL::Point<double> result(dimension, data);

  BOOST_CHECK_EQUAL(result.getDimension(), dimension);
  BOOST_CHECK_CLOSE(result[0], x, 1e-9);
  BOOST_CHECK_CLOSE(result[1], y, 1e-9);
  BOOST_CHECK_CLOSE(result[2], z, 1e-9);
}

BOOST_AUTO_TEST_CASE(array_weight) {
  int dimension = 3;
  double x = 2.5;
  double y = -5.7;
  double z = 8.3;

  double data[3];
  data[0] = x;
  data[1] = y;
  data[2] = z;
  double weight = 2.5;

  ALL::Point<double> result(dimension, data, weight);

  BOOST_CHECK_EQUAL(result.getDimension(), dimension);
  BOOST_CHECK_CLOSE(result[0], x, 1e-9);
  BOOST_CHECK_CLOSE(result[1], y, 1e-9);
  BOOST_CHECK_CLOSE(result[2], z, 1e-9);
  BOOST_CHECK_CLOSE(result.getWeight(), weight, 1e-9);
}

BOOST_AUTO_TEST_CASE(array_weight_index) {
  int dimension = 3;
  double x = 2.5;
  double y = -5.7;
  double z = 8.3;

  double data[3];
  data[0] = x;
  data[1] = y;
  data[2] = z;

  double weight = 1.4;
  long index = 1l;
  ALL::Point<double> result(dimension, data, weight, index);

  BOOST_CHECK_EQUAL(result.getDimension(), dimension);
  BOOST_CHECK_CLOSE(result[0], x, 1e-9);
  BOOST_CHECK_CLOSE(result[1], y, 1e-9);
  BOOST_CHECK_CLOSE(result[2], z, 1e-9);
  BOOST_CHECK_CLOSE(result.getWeight(), weight, 1e-9);
  BOOST_CHECK_EQUAL(result.get_id(), index);
}

BOOST_AUTO_TEST_CASE(std_vector) {
  int dimension = 3;
  double x = 2.5;
  double y = -5.7;
  double z = 8.3;
  std::vector<double> data(dimension);
  data.at(0) = x;
  data.at(1) = y;
  data.at(2) = z;
  ALL::Point<double> result(data);

  BOOST_CHECK_EQUAL(result.getDimension(), dimension);
  BOOST_CHECK_CLOSE(result[0], x, 1e-9);
  BOOST_CHECK_CLOSE(result[1], y, 1e-9);
  BOOST_CHECK_CLOSE(result[2], z, 1e-9);
}

BOOST_AUTO_TEST_CASE(vector_weight) {
  int dimension = 3;
  double x = 2.5;
  double y = -5.7;
  double z = 8.3;

  std::vector<double> data(dimension);
  data.at(0) = x;
  data.at(1) = y;
  data.at(2) = z;

  double weight = 1.4;
  ALL::Point<double> result(data, weight);

  BOOST_CHECK_EQUAL(result.getDimension(), dimension);
  BOOST_CHECK_CLOSE(result[0], x, 1e-9);
  BOOST_CHECK_CLOSE(result[1], y, 1e-9);
  BOOST_CHECK_CLOSE(result[2], z, 1e-9);
  BOOST_CHECK_CLOSE(result.getWeight(), weight, 1e-9);
}

BOOST_AUTO_TEST_CASE(vector_weight_index) {
  int dimension = 3;
  double x = 2.5;
  double y = -5.7;
  double z = 8.3;

  std::vector<double> data(dimension);
  data.at(0) = x;
  data.at(1) = y;
  data.at(2) = z;

  double weight = 1.4;
  long index = 1l;
  ALL::Point<double> result(data, weight, index);

  BOOST_CHECK_EQUAL(result.getDimension(), dimension);
  BOOST_CHECK_CLOSE(result[0], x, 1e-9);
  BOOST_CHECK_CLOSE(result[1], y, 1e-9);
  BOOST_CHECK_CLOSE(result[2], z, 1e-9);
  BOOST_CHECK_CLOSE(result.getWeight(), weight, 1e-9);
  BOOST_CHECK_EQUAL(result.get_id(), index);
}

BOOST_AUTO_TEST_SUITE_END()
