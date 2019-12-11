src = *.cpp

rasterizer_cpp:
	g++ $(src) -std=c++11 -g -o rasterizer
