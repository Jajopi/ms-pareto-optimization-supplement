CXX = g++
BOOST_INCLUDE ?= /opt/homebrew/Cellar/boost/1.90.0_1/include
INCLUDES = -I2i_bench/include -I2i_bench/external $(if $(BOOST_INCLUDE),-I$(BOOST_INCLUDE),)
CXXFLAGS = -std=c++17 -Wall -Wextra -O2 -DBOOST_TIMER_ENABLE_DEPRECATED $(INCLUDES)
TARGET = EliasFano
SOURCES = EliasFano.cpp
OBJECTS = $(SOURCES:.cpp=.o)

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJECTS) $(TARGET)

.PHONY: all clean
