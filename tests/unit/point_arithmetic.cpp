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

#define BOOST_TEST_MODULE point_arithmetic
#include "ALL_Point.hpp"
#include <vector>
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(point_arithmetic)

BOOST_AUTO_TEST_CASE(point_add_double) {
    std::vector<double> A_data = {2.8, -5.3, 8.8};
    std::vector<double> B_data = {-5.2, 1.1, 0.4};

    ALL::Point<double> A(A_data);
    ALL::Point<double> B(B_data);

    ALL::Point<double> C = A+B;

    BOOST_CHECK_CLOSE(C[0],-2.4,1e-9);
    BOOST_CHECK_CLOSE(C[1],-4.2,1e-9);
    BOOST_CHECK_CLOSE(C[2],9.2,1e-9);
}    

BOOST_AUTO_TEST_CASE(point_add_float) {
    std::vector<float> A_data = {2.8f, -5.3f, 8.8f};
    std::vector<float> B_data = {-5.2f, 1.1f, 0.4f};

    ALL::Point<float> A(A_data);
    ALL::Point<float> B(B_data);

    ALL::Point<float> C = A+B;

    BOOST_CHECK_CLOSE(C[0],-2.4f,1e-4);
    BOOST_CHECK_CLOSE(C[1],-4.2f,1e-4);
    BOOST_CHECK_CLOSE(C[2],9.2f,1e-4);
}    

BOOST_AUTO_TEST_CASE(point_add_long_double) {
    std::vector<long double> A_data = {2.8l, -5.3l, 8.8l};
    std::vector<long double> B_data = {-5.2l, 1.1l, 0.4l};

    ALL::Point<long double> A(A_data);
    ALL::Point<long double> B(B_data);

    ALL::Point<long double> C = A+B;

    BOOST_CHECK_CLOSE(C[0],-2.4l,1e-14);
    BOOST_CHECK_CLOSE(C[1],-4.2l,1e-14);
    BOOST_CHECK_CLOSE(C[2],9.2l,1e-14);
}    

BOOST_AUTO_TEST_CASE(point_sub_double) {
    std::vector<double> A_data = {2.8, -5.3, 8.8};
    std::vector<double> B_data = {-5.2, 1.1, 0.4};

    ALL::Point<double> A(A_data);
    ALL::Point<double> B(B_data);

    ALL::Point<double> C = A-B;

    BOOST_CHECK_CLOSE(C[0],8.0,1e-9);
    BOOST_CHECK_CLOSE(C[1],-6.4,1e-9);
    BOOST_CHECK_CLOSE(C[2],8.4,1e-9);
}    

BOOST_AUTO_TEST_CASE(point_cross_float) {
    std::vector<float> A_data = {2.8f, -5.3f, 8.8f};
    std::vector<float> B_data = {-5.2f, 1.1f, 0.4f};

    ALL::Point<float> A(A_data);
    ALL::Point<float> B(B_data);

    ALL::Point<float> C = A-B;

    BOOST_CHECK_CLOSE(C[0],8.0f,1e-4);
    BOOST_CHECK_CLOSE(C[1],-6.4f,1e-4);
    BOOST_CHECK_CLOSE(C[2],8.4f,1e-4);
}    

BOOST_AUTO_TEST_CASE(point_cross_long_double) {
    std::vector<long double> A_data = {2.8l, -5.3l, 8.8l};
    std::vector<long double> B_data = {-5.2l, 1.1l, 0.4l};

    ALL::Point<long double> A(A_data);
    ALL::Point<long double> B(B_data);

    ALL::Point<long double> C = A-B;

    BOOST_CHECK_CLOSE(C[0],8.0l,1e-14);
    BOOST_CHECK_CLOSE(C[1],-6.4l,1e-14);
    BOOST_CHECK_CLOSE(C[2],8.4l,1e-14);
}    

BOOST_AUTO_TEST_CASE(point_scale_double) {
    std::vector<double> A_data = {2.8, -5.3, 8.8};

    ALL::Point<double> A(A_data);
    double B = 0.5;

    ALL::Point<double> C1 = A*B;
    ALL::Point<double> C2 = B*A;

    BOOST_CHECK_CLOSE(C1[0],1.4,1e-9);
    BOOST_CHECK_CLOSE(C1[1],-2.65,1e-9);
    BOOST_CHECK_CLOSE(C1[2],4.4,1e-9);
    BOOST_CHECK_CLOSE(C2[0],1.4,1e-9);
    BOOST_CHECK_CLOSE(C2[1],-2.65,1e-9);
    BOOST_CHECK_CLOSE(C2[2],4.4,1e-9);
}    

BOOST_AUTO_TEST_CASE(point_scale_float) {
    std::vector<float> A_data = {2.8f, -5.3f, 8.8f};

    ALL::Point<float> A(A_data);
    float B = 0.5f;

    ALL::Point<float> C1 = A*B;
    ALL::Point<float> C2 = B*A;

    BOOST_CHECK_CLOSE(C1[0],1.4f,1e-4);
    BOOST_CHECK_CLOSE(C1[1],-2.65f,1e-4);
    BOOST_CHECK_CLOSE(C1[2],4.4f,1e-4);
    BOOST_CHECK_CLOSE(C2[0],1.4f,1e-4);
    BOOST_CHECK_CLOSE(C2[1],-2.65f,1e-4);
    BOOST_CHECK_CLOSE(C2[2],4.4f,1e-4);
}    

BOOST_AUTO_TEST_CASE(point_scale_long_double) {
    std::vector<long double> A_data = {2.8, -5.3, 8.8};

    ALL::Point<long double> A(A_data);
    long double B = 0.5l;

    ALL::Point<long double> C1 = A*B;
    ALL::Point<long double> C2 = B*A;

    BOOST_CHECK_CLOSE(C1[0],1.4l,1e-14);
    BOOST_CHECK_CLOSE(C1[1],-2.65l,1e-14);
    BOOST_CHECK_CLOSE(C1[2],4.4l,1e-14);
    BOOST_CHECK_CLOSE(C2[0],1.4l,1e-14);
    BOOST_CHECK_CLOSE(C2[1],-2.65l,1e-14);
    BOOST_CHECK_CLOSE(C2[2],4.4l,1e-14);
}    

BOOST_AUTO_TEST_CASE(point_dist_plane) {
    std::vector<double> A_data = {1.,0.,0.};
    std::vector<double> B_data = {0.,1.,0.};
    std::vector<double> C_data = {0.,0.,1.};
    std::vector<double> p_data = {0.5,0.5,0.5};
    ALL::Point<double> A(A_data);
    ALL::Point<double> B(B_data);
    ALL::Point<double> C(C_data);
    ALL::Point<double> p(p_data);

    double dist = p.dist_plane(A,B,C);

    BOOST_CHECK_CLOSE(dist, 0.28867513459481292, 1e-14);
}

BOOST_AUTO_TEST_SUITE_END()

