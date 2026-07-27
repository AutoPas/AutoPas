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

#define BOOST_TEST_MODULE all_class_creation
#include "ALL.hpp"
#include "ALL_Point.hpp"
#include "ALL_CustomExceptions.hpp"
#include <boost/test/unit_test.hpp>
#include <list>
#include <mpi.h>
#include <vector>

BOOST_AUTO_TEST_SUITE(all_class_creation)

BOOST_AUTO_TEST_CASE(empty) {
  int already_init;
  MPI_Initialized(&already_init);
  if (!already_init)
    MPI_Init(NULL, NULL);

  ALL::ALL<double, double> test();

  int already_final;
  MPI_Finalized(&already_final);
  if (!already_final)
    MPI_Finalize();
}

BOOST_AUTO_TEST_CASE(simple) {
  int already_init;
  MPI_Initialized(&already_init);
  if (!already_init)
    MPI_Init(NULL, NULL);

  ALL::LB_t method = ALL::TENSOR;
  int dimension = 3;
  double gamma = 4.0;

  ALL::ALL<double, double> test(method, dimension, gamma);

  // comment(s.schulz): What is supposed to be checked here?
  //                    Just the calling? Then we should make sure to also access the result
  double check_dim = test.getDimension();
  double check_gamma = test.getGamma();

  // comment(s.schulz): Not sure if there is a better way than assert in the boost framework
  assert(check_dim == dimension);
  assert(check_gamma == gamma);

  int already_final;
  MPI_Finalized(&already_final);
  if (!already_final)
    MPI_Finalize();
}

// comment (r.halver): check if the constructor fails when using a dimension not equal to 3
// !to be changed, when the library can deal with dimensions not equal to 3!
BOOST_AUTO_TEST_CASE(dimension_2)    
{
    int already_init;
    MPI_Initialized(&already_init);
    if (!already_init)
        MPI_Init(NULL,NULL);

    ALL::LB_t method = ALL::TENSOR;
    int dimension = 2;
    double gamma = 4.0;

    bool success = false;
    try
    {
        ALL::ALL<double, double> test(method, dimension, gamma);
    }
    catch(ALL::InvalidArgumentException)
    {
        success = true;
    }
    BOOST_TEST(success); 
}

BOOST_AUTO_TEST_SUITE_END()
