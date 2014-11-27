#    Copyright (C) 2014 University of Southern California,
#                       Andrew D. Smith and Ehsan Behnam
#
#    Authors: Andrew D. Smith
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

PROGS = build_graph_naively naive_batch_query

SMITHLAB_CPP = ../smithlab_cpp
COMMON = ../common

INCLUDEDIRS = $(SMITHLAB_CPP) $(COMMON)
INCLUDEDIRS += $(BOOST)
INCLUDEARGS = $(addprefix -I,$(INCLUDEDIRS))

CXX = g++
CXXFLAGS = -Wall -fmessage-length=50
OPTFLAGS = -O2
DEBUGFLAGS = -g -pg

# Flags for Mavericks
ifeq "$(shell uname)" "Darwin"
CXXFLAGS += -arch x86_64
ifeq "$(shell if [ `sysctl -n kern.osrelease | cut -d . -f 1` -ge 13 ];\
              then echo 'true'; fi)" "true"
CXXFLAGS += -stdlib=libstdc++
else
# Is OMP only implemented for non-Mavericks?
CXXFLAGS += -DHAVE_OMP
LIBS += -fopenmp
endif
endif

ifdef DEBUG
CXXFLAGS += $(DEBUGFLAGS)
endif

ifdef OPT
CXXFLAGS += $(OPTFLAGS)
endif

all: $(PROGS)

build_graph_naively : $(addprefix $(COMMON)/, RegularNearestNeighborGraph.o)

amordad_batch_query : $(addprefix $(COMMON)/, RegularNearestNeighborGraph.o \
	LSHAngleHashTable.o LSHAngleHashFunction.o)

$(PROGS): $(addprefix $(SMITHLAB_CPP)/, smithlab_os.o \
	smithlab_utils.o OptionParser.o) \
	$(addprefix $(COMMON)/, FeatureVector.o)

%.o: %.cpp %.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(INCLUDEARGS)

%: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^ $(INCLUDEARGS) $(LIBS)

install: all
	@mkdir -p ../../bin
	@install -m 755 $(PROGS) ../../bin

clean:
	@-rm -f $(PROGS) *.o *.so *.a

.PHONY: clean
