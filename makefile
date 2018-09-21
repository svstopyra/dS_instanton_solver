CXX = g++
CXXFLAGS = -g -std=c++11
LIB = ar
LIBFLAGS = -rcT
#Set the following to the location of the boost libraries on your system:
BOOST_LOCATION = /export/zupcx28/sstopyra/boost_1_60_0

#Set the following to the location of MPFR on your system:
MPFR_INCLUDE_LOCATION = /export/zupcx28/sstopyra/mpfr-build/include
MPFR_LIB_LOCATION = /export/zupcx28/sstopyra/mpfr-build/lib

#Includes. NB - crucial that odeint_with_event_detection is included BEFORE boost, as it overrides the default ode solver.
INCLUDE_COMMAND = -Iodeint_with_event_detection -I$(BOOST_LOCATION) 

#To compile with double precision, set the following to TRUE (default compilation is 100dp multi-precision).
DOUBLE_PRECISION = FALSE
ifeq ($(DOUBLE_PRECISION),TRUE)
CXXFLAGS += -DUSE_DOUBLE_PRECISION
endif
ifeq ($(DOUBLE_PRECISION),FALSE)
LINK_COMMAND = -L$(MPFR_LIB_LOCATION) -lmpfr
INCLUDE_COMMAND += -I$(MPFR_INCLUDE_LOCATION)
endif
#For more options (e.g. specifying custom float types, different level of precision etc...), edit multi_precision_definitions.h

dSSolverLibrary: dSSolverLibrary_intermediate.a AdS.o dS.o Flat.o multi_precision_definitions.o potentials.o
	$(LIB) $(LIBFLAGS) dSSolverLibrary.a dSSolverLibrary_intermediate.a $(MPFR_LIB_LOCATION)/libmpfr.a
	rm dSSolverLibrary_intermediate.a


dSSolverLibrary_intermediate.a: AdS.o dS.o Flat.o multi_precision_definitions.o potentials.o
	$(LIB) $(LIBFLAGS) dSSolverLibrary_intermediate.a AdS.o dS.o Flat.o multi_precision_definitions.o potentials.o

AdS.o: AdS.cpp project_specific.h multi_precision_definitions.h AdS.h
	$(CXX) $(CXXFLAGS) -c AdS.cpp $(INCLUDE_COMMAND)

dS.o: dS.cpp project_specific.h multi_precision_definitions.h dS.h find_zero.h transition_solver.h numerical_integration.h
	$(CXX) $(CXXFLAGS) -c dS.cpp $(INCLUDE_COMMAND)

Flat.o: Flat.cpp project_specific.h multi_precision_definitions.h Flat.h
	$(CXX) $(CXXFLAGS) -c Flat.cpp $(INCLUDE_COMMAND)

multi_precision_definitions.o: multi_precision_definitions.cpp multi_precision_definitions.h
	$(CXX) $(CXXFLAGS) -c multi_precision_definitions.cpp $(INCLUDE_COMMAND)

#NewtonRaphson.o: NewtonRaphson.cpp project_specific.h Newton_Raphson_Code.h
#	$(CXX) $(CXXFLAGS) -c NewtonRaphson.cpp $(INCLUDE_COMMAND)

potentials.o: potentials.cpp project_specific.h potentials_code.h
	$(CXX) $(CXXFLAGS) -c potentials.cpp $(INCLUDE_COMMAND)

clean:
	rm AdS.o dS.o Flat.o multi_precision_definitions.o potentials.o dSSolverLibrary_intermediate.a dSSolverLibrary.a

