#
# プログラム名
#
PROG = schrac

#
# ソースコードが存在する相対パス
#
VPATH = src src/checkpoint

#
# コンパイル対象のソースファイル群（カレントディレクトリ以下の*.cppファイル）
#
SRCS = $(shell find * -name "*.cpp")

#
# ターゲットファイルを生成するために利用するオブジェクトファイル
#
OBJDIR = 
ifeq "$(strip $(OBJDIR))" ""
  OBJDIR = .
endif

OBJS = $(addprefix $(OBJDIR)/, $(notdir $(SRCS:.cpp=.o)))

#
# *.cppファイルの依存関係が書かれた*.dファイル
#
DEPS = $(OBJS:.o=.d)

#
# C++コンパイラの指定
#
CXX = g++

#
# C++コンパイラに与える、（最適化等の）オプション
#
CXXFLAGS = -Wall -Wextra -O3 -std=c++17 -mtune=native -march=native

#
# リンク対象に含めるライブラリの指定
#
LDFLAGS = -L/home/dc1394/oss/boost_1_75_0/stage/gcc/lib -lboost_program_options \
		  -L/home/dc1394/oss/tbb/lib/intel64/gcc4.8 -ltbb \
		  -lgsl -lgslcblas -lm

#
# makeの動作
#
all: $(PROG) ; rm -f $(OBJS) $(DEPS)

#
# 依存関係を解決するためのinclude文
#
-include $(DEPS)

#
# プログラムのリンク
#
$(PROG): $(OBJS)
		$(CXX) $^ $(LDFLAGS) $(CXXFLAGS) -o $@

#
# プログラムのコンパイル
#
%.o: %.cpp
		$(CXX) $(CXXFLAGS) -c -MMD -MP $<

#
# make cleanの動作
#
clean:
		rm -f $(PROG) $(OBJS) $(DEPS)
