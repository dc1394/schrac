PROG := schrac
SRCS :=	ci_string.cpp data.cpp diffdata.cpp diffsolver.cpp diracnormalize.cpp \
		eigenvaluesearch.cpp energy.cpp getcomlineoption.cpp goexit.cpp \
		normalization.cpp readinputfile.cpp rho.cpp scfloop.cpp schnormalize.cpp \
		schracmain.cpp simpson.cpp solvelinearequ.cpp vhartree.cpp \
		wavefunctionsave.cpp checkpoint.cpp
OBJS :=	$(SRCS:%.cpp=%.o)
DEPS :=	$(SRCS:%.cpp=%.d)

VPATH  = src src/checkpoint
CXX = clang++
CXXFLAGS = -Wextra -O3 -pipe -std=c++11
LDFLAGS = -L/home/dc1394/oss/boost_1_57_0/stage/clang/lib/ -lboost_program_options \
		  -lgsl -lgslcblas -lm -L/home/dc1394/oss/tbb43_20150316oss/lib/intel64/gcc4.4 -ltbb

all: $(PROG) ; rm -f $(OBJS) $(DEPS)

-include $(DEPS)

$(PROG): $(OBJS)
		$(CXX) $(LDFLAGS) $(CXXFLAGS) -o $@ $^

%.o: %.cpp
		$(CXX) $(CXXFLAGS) -c -MMD -MP $<

clean:
		rm -f $(PROG) $(OBJS) $(DEPS)
