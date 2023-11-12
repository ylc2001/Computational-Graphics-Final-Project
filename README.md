### Introduction

A simple ray tracing rendering engine was implemented by choosing the **Path Tracing** algorithm. In terms of light propagation, three materials are implemented: **Diffuse Reflection, Specular Reflection, and Refraction**, and the nature of light propagation on the surface of an object can be mixed by these three. Supported geometries are **spheres, planes, rotated surfaces, imported triangular mesh models**, and they can be linearly transformed to change size/direction/position. Implemented analytic intersection of Bezier and B-spline **rotated surfaces** based on Newton's iterative method.

The **depth-of-field effect** is implemented with adjustable focus and aperture size; implemented **antialiasing** by adding random perturbations to the samples; mapping of spheres and rotated surfaces is implemented using UV expansion.

In terms of algorithmic efficiency, **Axis-aligned bounding box** acceleration is implemented for complex meshes and rotated surfaces, and **OpenMP** is used for hardware acceleration to take full advantage of multiple CPU cores programmatically.

### Project Framework

This assignment is implemented on the basis of the code framework of PA1, and the file structure is close to it; compiled using Makefile, run make in the root directory of the project to compile, and it can be successfully compiled and run under Windows system or Linux system. In essence, the compilation commands are as follows:

```shell
g++ src/image.cpp src/main.cpp src/mesh.cpp src/scene_parser.cpp src/texture.cpp src/vecmath.cpp -o main -O3 -Wall -fopenmp -I include/
```

The code framework is described below: 

- include/: All .h or .hpp files that contain definitions of classes for individual objects and objects, as well as the main logic of the algorithm. 
- mesh/: Triangular mesh model file resources in .obj format 
- output/: The directory where the output image is stored 
- src/: All .cpp files, containing the main function code as well as the code related to image storage, scene file reading, triangular mesh model file and texture file reading, and the 3D vector arithmetic library provided by PA1. 
- testcases/: Scene description files needed to generate results
- texture/: Stores mapping resources 

Usage：`./main <input scene file> <output bmp file> <method> <spp>`

where `method` refers to the algorithm used by the program, either rc (ray casting) or pt (path tracing). 

Usage example：` ./main testcases/test01.txt output/test01.bmp pt 15000`

## Algorithm & implementation

### Path Tracing

The idea of the path tracing algorithm is roughly as follows: a ray is fired from the camera via a specified pixel, and if it fails to hit an object, the background color of the scene is used as the color of that pixel. In general, the ray is continuously reflected and refracted by the surface of the object in the scene. If the light source can be hit, the intensity of the light source is multiplied by the reflection color (attenuation rate) of the objects it passes through as the contribution of the light source to the pixel's color; otherwise, when the pre-set iteration depth is exceeded, the iteration stops and all the previously accumulated light intensity is counted as the pixel's color. Each pixel is sampled multiple times (sampling rate) and the result of each calculation is averaged as the final color value for that pixel.

The code for the path tracing algorithm is mainly in render.hpp, and its main logical parts are as follows:

```c++
void render() // called by main() after reading necessary infos
{
    PerspectiveCamera *camera = (PerspectiveCamera *)sceneParser.getCamera();
    int w = camera->getWidth(), h = camera->getHeight();
    Image renderedImage(w, h);
#pragma omp parallel for schedule(dynamic, 1) // OpenMP
    for (int x = 0; x < w; ++x) {
        unsigned short Xi[3] = {0, 0, (unsigned short)(x * x * x)};
        for (int y = 0; y < h; ++y) {
            Vector3f avg_color = Vector3f::ZERO;
            for (int s = 0; s < samps; ++s) { // sample rate
                Ray camRay = camera->generateRay(Vector2f(x + erand48(Xi) - 0.5, y + erand48(Xi) - 0.5)); // generate a ray
                avg_color += Color_func(camRay, sceneParser, Xi); // calculate corresponding color value
            }
            renderedImage.SetPixel(x, y, avg_color / samps); // avg to decide the final color
        }
    }
    renderedImage.SaveImage(fout);
}
```

## Efficiency

### Algorithm acceleration

For triangular meshes and rotated surfaces, both of which are computationally intensive in light intersection, I have implemented axis-aligned box acceleration. A rectangular box aligned with the x, y, and z axes is used to enclose the target object. Intersections are made with the rectangular box, and if there are no intersections, the box will not be intersected by the object it encircles. This reduces a lot of unnecessary calculations and speeds up the rendering of scenes with triangulated meshes and rotated surfaces. It has been tested to reduce rendering time by up to 64% for a given scene with box-around acceleration. The relevant code is located in mesh.cpp and revsurface.hpp.

### Hardware acceleration

Using OpenMP, a lightweight multi-thread library in C++, it is possible to parallelize our program. For image rendering tasks, since each pixel's color values are computed independently of each other, they can be evenly distributed among different threads running. Basically, the computer's CPU has as many cores as it has cores, and the performance is doubled, allowing you to take full advantage of the CPU's performance, with significant acceleration.

