PROG := schrac
SRCS :=	schracmain.cpp my_getopt.cpp 
OBJS :=	$(SRCS:%.cpp=%.o)
DEPS :=	$(SRCS:%.cpp=%.d)

CXX = icc
CXXFLAGS = -Wextra -O3 -pipe -std=c++11
LDFLAGS = -L/home/dc1394/oss/boost_1_57_0/stage/icc/lib/ -lboost_program_options

all: $(PROG)

-include $(DEPS)

$(PROG): $(OBJS)
		$(CXX) $(LDFLAGS) $(CXXFLAGS) -o $@ $^

%.o: %.cpp
		$(CXX) $(CXXFLAGS) -c -MMD -MP $<

clean:
		rm -f $(PROG) $(OBJS) $(DEPS)
