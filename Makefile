src = *.cpp

all:
	g++ $(src) -std=c++11 -g -o rasterizer
