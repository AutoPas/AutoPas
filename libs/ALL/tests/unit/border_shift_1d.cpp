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

#define BOOST_TEST_MODULE borderShift1d
#include "ALL_Functions.hpp"
#include <boost/test/unit_test.hpp>
#include <iostream>
#include <list>
#include <mpi.h>
#include <vector>

BOOST_AUTO_TEST_SUITE(borderShift1d)

using ALL::Functions::borderShift1d;

BOOST_AUTO_TEST_CASE(double_normal_case) {

  int remote_rank = 384;
  int local_coord = 3;
  int global_dim = 8;
  double local_work = 0.6;
  double remote_work = 0.4;
  double local_size = 1.0;
  double remote_size = 1.0;
  double gamma = 4.0;
  double min_size = 0.2;

  double shift =
      borderShift1d(remote_rank, local_coord, global_dim, local_work,
                    remote_work, local_size, remote_size, gamma, min_size);

  BOOST_CHECK_CLOSE(shift, -0.05, 1e-9);
}

BOOST_AUTO_TEST_CASE(double_neighbor_null) {

  int remote_rank = MPI_PROC_NULL;
  int local_coord = 3;
  int global_dim = 8;
  double local_work = 0.6;
  double remote_work = 0.4;
  double local_size = 1.0;
  double remote_size = 1.0;
  double gamma = 4.0;
  double min_size = 0.2;

  double shift =
      borderShift1d(remote_rank, local_coord, global_dim, local_work,
                    remote_work, local_size, remote_size, gamma, min_size);

  BOOST_CHECK_CLOSE(shift, 0.0, 1e-9);
}

BOOST_AUTO_TEST_CASE(double_neighbor_larger_than_min_size) {

  int remote_rank = 384;
  int local_coord = 3;
  int global_dim = 8;
  double local_work = 0.6;
  double remote_work = 0.4;
  double local_size = 1.0;
  double remote_size = 1.0;
  double gamma = 4.0;
  double min_size = 0.98;

  double shift =
      borderShift1d(remote_rank, local_coord, global_dim, local_work,
                    remote_work, local_size, remote_size, gamma, min_size);

  BOOST_CHECK_CLOSE(shift, -0.02 * 0.49, 1e-9);
}

BOOST_AUTO_TEST_CASE(double_local_process_on_edge) {

  int remote_rank = 384;
  int local_coord = 3;
  int global_dim = 4;
  double local_work = 0.6;
  double remote_work = 0.4;
  double local_size = 1.0;
  double remote_size = 1.0;
  double gamma = 4.0;
  double min_size = 0.98;

  double shift =
      borderShift1d(remote_rank, local_coord, global_dim, local_work,
                    remote_work, local_size, remote_size, gamma, min_size);

  BOOST_CHECK_CLOSE(shift, 0, 1e-9);
}

BOOST_AUTO_TEST_SUITE_END()
