all:
	g++ src/image.cpp src/main.cpp src/mesh.cpp src/scene_parser.cpp src/texture.cpp src/vecmath.cpp -o main -O3 -Wall -fopenmp -I include/ -static