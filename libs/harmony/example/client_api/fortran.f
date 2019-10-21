c
c Copyright 2003-2016 Jeffrey K. Hollingsworth
c
c This file is part of Active Harmony.
c
c Active Harmony is free software: you can redistribute it and/or modify
c it under the terms of the GNU Lesser General Public License as published
c by the Free Software Foundation, either version 3 of the License, or
c (at your option) any later version.
c
c Active Harmony is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c GNU Lesser General Public License for more details.
c
c You should have received a copy of the GNU Lesser General Public License
c along with Active Harmony.  If not, see <http://www.gnu.org/licenses/>.
c
      program example

      integer i

c     simple performance function array
      integer hd

      integer y(20)

      integer x, int_type, perf

      int_type = 1

c     populate the performance array
      do 40 i=1, 20
         y(i)=(i-9)*(i-9)*(i-9)-75*(i-9)+24
 40   continue

c     initialize the harmony server
      call harmony_startup_f(hd)

c     send in the application description
      call harmony_application_setup_file_f('./simplext.tcl'//CHAR(0),
     &  hd)

c     register the tunable parameter(s)
      call harmony_add_variable_int_f(hd, 'SimplexT'//CHAR(0),
     &  'x'//CHAR(0), 0, x)

c     main loop
      do 30 i = 1, 200
c     update the performance result
         call harmony_performance_update_int_f(hd, y(x))
c     update the variable x
         call harmony_request_variable_f(hd, 'x'//CHAR(0), x)

 30   continue

c     close the session
      call harmony_end_f(hd)
      end program example
