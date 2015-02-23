PROG := schrac
SRCS :=	ci_string.cpp data.cpp diffdata.cpp diffsolver.cpp diracnormalize.cpp \
		eigenvaluesearch.cpp energy.cpp getcomlineoption.cpp normalization.cpp \
		readinputfile.cpp schnormalize.cpp schracmain.cpp simpson.cpp \
		wavefunctionsave.cpp
OBJS :=	$(SRCS:%.cpp=%.o)
DEPS :=	$(SRCS:%.cpp=%.d)

VPATH  = src
CXX = icc
CXXFLAGS = -Wextra -O3 -pipe -std=c++11
LDFLAGS = -L/home/dc1394/oss/boost_1_57_0/stage/icc/lib/ -lboost_program_options \
		  -lgsl -lgslcblas -ltbb

all: $(PROG) ; rm -f $(OBJS) $(DEPS)

-include $(DEPS)

$(PROG): $(OBJS)
		$(CXX) $(LDFLAGS) $(CXXFLAGS) -o $@ $^

%.o: %.cpp
		$(CXX) $(CXXFLAGS) -c -MMD -MP $<

clean:
		rm -f $(PROG) $(OBJS) $(DEPS)
