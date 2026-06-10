/*
Copyright 2018-2020 Rene Halver, Forschungszentrum Juelich GmbH, Germany
Copyright 2018-2020 Godehard Sutmann, Forschungszentrum Juelich GmbH, Germany

Redistribution and use in source and binary forms, with or without modification, 
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this 
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this 
   list of conditions and the following disclaimer in the documentation and/or 
   other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may 
   be used to endorse or promote products derived from this software without 
   specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND 
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. 
IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, 
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, 
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, 
OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
POSSIBILITY OF SUCH DAMAGE.
*/

#define BOOST_TEST_MAIN

#define BOOST_TEST_MODULE point_tetrahedron
#include "ALL_Point.hpp"
#include <vector>
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(point_cross)

BOOST_AUTO_TEST_CASE(point_cross_double) {
    std::vector<double> A_data = {2.8, -5.3, 8.8};
    std::vector<double> B_data = {-5.2, 1.1, 0.4};

    ALL::Point<double> A(A_data);
    ALL::Point<double> B(B_data);

    ALL::Point<double> C = A.cross(B);

    BOOST_CHECK_CLOSE(C[0],-11.8,1e-9);
    BOOST_CHECK_CLOSE(C[1],-46.88,1e-9);
    BOOST_CHECK_CLOSE(C[2],-24.48,1e-9);
}    

BOOST_AUTO_TEST_CASE(point_cross_float) {
    std::vector<float> A_data = {2.8f, -5.3f, 8.8f};
    std::vector<float> B_data = {-5.2f, 1.1f, 0.4f};

    ALL::Point<float> A(A_data);
    ALL::Point<float> B(B_data);

    ALL::Point<float> C = A.cross(B);

    BOOST_CHECK_CLOSE(C[0],-11.8f,1e-4);
    BOOST_CHECK_CLOSE(C[1],-46.88f,1e-4);
    BOOST_CHECK_CLOSE(C[2],-24.48f,1e-4);
}    

BOOST_AUTO_TEST_CASE(point_cross_long_double) {
    std::vector<long double> A_data = {2.8l, -5.3l, 8.8l};
    std::vector<long double> B_data = {-5.2l, 1.1l, 0.4l};

    ALL::Point<long double> A(A_data);
    ALL::Point<long double> B(B_data);

    ALL::Point<long double> C = A.cross(B);

    BOOST_CHECK_CLOSE(C[0],-11.8l,1e-14);
    BOOST_CHECK_CLOSE(C[1],-46.88l,1e-14);
    BOOST_CHECK_CLOSE(C[2],-24.48l,1e-14);
}    

BOOST_AUTO_TEST_SUITE_END()

