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

#define BOOST_TEST_MODULE point_tetrahedron
#include "ALL_Point.hpp"
#include <boost/test/unit_test.hpp>
#include <vector>

BOOST_AUTO_TEST_SUITE(point_tetrahedron)

BOOST_AUTO_TEST_CASE(tetrahedron_double) {
  int dimension = 3;
  double x = 2.5;
  double y = -5.7;
  double z = 8.3;

  std::vector<double> A_data({0.0, 0.0, 0.0});
  std::vector<double> B_data({1.0, 0.0, 0.0});
  std::vector<double> C_data({0.0, 1.0, 0.0});
  std::vector<double> D_data({0.0, 0.0, 1.0});
  std::vector<double> P_data({0.25, 0.25, 0.25});

  ALL::Point<double> A(A_data);
  ALL::Point<double> B(B_data);
  ALL::Point<double> C(C_data);
  ALL::Point<double> D(D_data);
  ALL::Point<double> P(P_data);

  BOOST_CHECK(P.inTetrahedron(A, B, C, D));
}

BOOST_AUTO_TEST_CASE(tetrahedron_float) {
  int dimension = 3;
  float x = 2.5f;
  float y = -5.7f;
  double z = 8.3f;

  std::vector<float> A_data({0.0f, 0.0f, 0.0f});
  std::vector<float> B_data({1.0f, 0.0f, 0.0f});
  std::vector<float> C_data({0.0f, 1.0f, 0.0f});
  std::vector<float> D_data({0.0f, 0.0f, 1.0f});
  std::vector<float> P_data({0.25f, 0.25f, 0.25f});

  ALL::Point<float> A(A_data);
  ALL::Point<float> B(B_data);
  ALL::Point<float> C(C_data);
  ALL::Point<float> D(D_data);
  ALL::Point<float> P(P_data);

  BOOST_CHECK(P.inTetrahedron(A, B, C, D));
}

BOOST_AUTO_TEST_CASE(tetrahedron_long_double) {
  int dimension = 3;
  long double x = 2.5l;
  long double y = -5.7l;
  long double z = 8.3l;

  std::vector<long double> A_data({0.0l, 0.0l, 0.0l});
  std::vector<long double> B_data({1.0l, 0.0l, 0.0l});
  std::vector<long double> C_data({0.0l, 1.0l, 0.0l});
  std::vector<long double> D_data({0.0l, 0.0l, 1.0l});
  std::vector<long double> P_data({0.25l, 0.25l, 0.25l});

  ALL::Point<long double> A(A_data);
  ALL::Point<long double> B(B_data);
  ALL::Point<long double> C(C_data);
  ALL::Point<long double> D(D_data);
  ALL::Point<long double> P(P_data);

  BOOST_CHECK(P.inTetrahedron(A, B, C, D));
}

BOOST_AUTO_TEST_CASE(tetrahedron_double_on_edge) {
  int dimension = 3;
  double x = 2.5;
  double y = -5.7;
  double z = 8.3;

  std::vector<double> A_data({0.0, 0.0, 0.0});
  std::vector<double> B_data({1.0, 0.0, 0.0});
  std::vector<double> C_data({0.0, 1.0, 0.0});
  std::vector<double> D_data({0.0, 0.0, 1.0});
  std::vector<double> P_data({0.0, 0.5, 0.0});

  ALL::Point<double> A(A_data);
  ALL::Point<double> B(B_data);
  ALL::Point<double> C(C_data);
  ALL::Point<double> D(D_data);
  ALL::Point<double> P(P_data);

  BOOST_CHECK(!P.inTetrahedron(A, B, C, D));
}

BOOST_AUTO_TEST_CASE(tetrahedron_double_outside) {
  int dimension = 3;
  double x = 2.5;
  double y = -5.7;
  double z = 8.3;

  std::vector<double> A_data({0.0, 0.0, 0.0});
  std::vector<double> B_data({1.0, 0.0, 0.0});
  std::vector<double> C_data({0.0, 1.0, 0.0});
  std::vector<double> D_data({0.0, 0.0, 1.0});
  std::vector<double> P_data({0.25, 0.25, 0.51});

  ALL::Point<double> A(A_data);
  ALL::Point<double> B(B_data);
  ALL::Point<double> C(C_data);
  ALL::Point<double> D(D_data);
  ALL::Point<double> P(P_data);

  BOOST_CHECK(!P.inTetrahedron(A, B, C, D));
}

BOOST_AUTO_TEST_SUITE_END()
