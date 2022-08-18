#
#  Copyright Daniele Melini & Giorgio Spada, 2014-2022
#
#  This file is part of REAR.
#
#  REAR is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  REAR is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with REAR.  If not, see <http://www.gnu.org/licenses/>.
#
#
FC=gfortran
FOPTS=-fopenmp -O
#
GF_OBJS=make_gf.o utils.o plegendre.o icosahedron.o
MAP_OBJS=make_map.o utils.o icosahedron.o
#
all: ../make_gf.exe ../make_map.exe
#
clean:
	$(RM) ../make_gf.exe ../make_map.exe *.o
#
../make_gf.exe: $(GF_OBJS)
	$(FC) $(FOPTS) $(GF_OBJS) -o ../make_gf.exe

../make_map.exe: $(MAP_OBJS)
	$(FC) $(FOPTS) $(MAP_OBJS) -o ../make_map.exe
#
%.o: %.f
	$(FC) $(FOPTS) -w -c $^ -o $@
#
%.o: %.f90
	$(FC) $(FOPTS) -c $^ -o $@
#
%.o: %.f95
	$(FC) $(FOPTS) -c $^ -o $@
#
