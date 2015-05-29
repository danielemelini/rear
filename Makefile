#
#  Copyright Daniele Melini & Giorgio Spada, 2014
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
all: banner
	@echo "  Use:"
	@echo "  	 ""make gf""  to compute Green's Functions"
	@echo "  	 ""make map"" to compute predictions of geodetic observables"
	@echo ""
	@echo "  Fortran compiler: $(FC)"
	@echo "  Compiler options: $(FOPTS)"
	@echo ""

banner:
	@echo ""
	@echo "  <<<<<<<<  REAR: a Regional ElAstic Rebound calculator  >>>>>>>>"
	@echo ""	

gf: banner
	@echo "  ---> Compiling make_gf.f90 ..."
	@echo ""
	@$(FC) $(FOPTS) make_gf.f90 -o make_gf.exe
	@echo ""
	@echo "  ---> Running make_gf.exe ..."
	@echo ""
	@./make_gf.exe
	@$(RM) make_gf.exe
	@echo ""

map: banner
	@echo "  ---> Compiling make_map.f90 ..."
	@echo ""
	@$(FC) $(FOPTS) make_map.f90 -o make_map.exe
	@echo ""
	@echo "  ---> Running make_map.exe ..."
	@echo ""
	@./make_map.exe
	@$(RM) make_map.exe
	@echo ""


