###################################################################
#  Makefile for project 2
#
#  Daniel R. Reynolds
#  SMU Mathematics
#  Math 3316
#  16 September 2015
###################################################################

# compilation definitions
CXX = g++
CXXFLAGS = -O -std=c++11

# makefile targets
all : test_bisect.exe 

test_bisect.exe : test_bisect.cpp bisection.cpp
	${CXX} ${CXXFLAGS} $^ -o $@

clean :
	\rm -f *.o *.txt

realclean : clean
	\rm -f *.exe *~


####### End of Makefile #######
