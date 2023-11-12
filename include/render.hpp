#ifndef PATH_TRACER_H
#define PATH_TRACER_H

#include <cmath>
#include <cstring>
#include <iostream>
#include <string>
#include <ctime>

#include "camera.hpp"
#include "constants.h"
#include "group.hpp"
#include "hit.hpp"
#include "image.hpp"
#include "light.hpp"
#include "ray.hpp"
#include "scene_parser.hpp"
#include "utils.hpp"
using namespace std;

std::default_random_engine generator;
std::uniform_real_distribution<double> distr(0.0, 1.0);
double erand48(unsigned short *X)
{
    return distr(generator);
}
double erand48()
{
    return distr(generator);
}

// Ray Casting 光线投射 基础算法
static Vector3f rcColor(Ray ray, const SceneParser &sceneParser, unsigned short *Xi)
{
    Group *baseGroup = sceneParser.getGroup();
    Hit hit;
    // 判断ray是否和场景有交点，并返回最近交点的数据，存储在hit中
    bool isIntersect = baseGroup->intersect(ray, hit, 0);
    if (isIntersect)
    {
        Vector3f color = Vector3f::ZERO;
        // 找到交点之后，累加来自所有光源的光强影响
        for (int li = 0; li < sceneParser.getNumLights(); ++li)
        {
            Light *light = sceneParser.getLight(li);
            Vector3f L, lightColor;
            // 获得光照强度
            light->getIllumination(ray.pointAtParameter(hit.getT()), L, lightColor);
            // 计算局部光强
            color += hit.getMaterial()->Shade(ray, hit, L, lightColor);
        }
        return color;
    }
    else
    {
        // 不存在交点，返回背景色
        return sceneParser.getBackgroundColor();
    }
}

// Path Tracing 路径追踪 from smallpt
static Vector3f ptColor(Ray ray, const SceneParser &scene, unsigned short *Xi)
{
    Group *group = scene.getGroup();
    int depth = 0;
    Vector3f color(0, 0, 0), cf(1, 1, 1);

    while (true)
    {
        if (++depth > TRACE_DEPTH)
        {
            if (erand48(Xi) < cf.max())
                cf = cf * cf.max();
            else
                return color;
        }
        // 判断camRay是否和场景有交点,返回最近交点的数据,存储在hit中.
        Hit hit;
        if (!group->intersect(ray, hit, 0))
        {
            color += scene.getBackgroundColor();
            return color;
        }

        // Path Tracing
        ray.origin += ray.direction * hit.getT();
        Material *material = hit.getMaterial();
        Vector3f refColor(hit.color);
        Vector3f N(hit.getNormal());

        Vector3f nl = Vector3f::dot(N, ray.getDirection()) < 0 ? N : N * -1;

        // Emission
        color += material->emission * cf;
        cf = cf * refColor;
        float random_type = erand48(Xi);
        if (random_type <= material->type.x())
        { // diffuse
            double r1 = 2 * M_PI * erand48(Xi), r2 = erand48(Xi), r2s = sqrt(r2);
            Vector3f w = nl, u = ((fabs(w.x()) > .1 ? Vector3f(0, 1, 0) : Vector3f(1, 0, 0)) % w).normalized(), v = w % u;
            ray.direction = (u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2)).normalized();
        }
        else if (random_type <= material->type.y() + material->type.x())
        { // specular
            float cost = Vector3f::dot(ray.direction, N);
            ray.direction = (ray.direction - N * (cost * 2)).normalized();
        }
        else
        { // refraction
            float n = material->refr;
            float R0 = ((1.0 - n) * (1.0 - n)) / ((1.0 + n) * (1.0 + n));
            if (Vector3f::dot(N, ray.direction) > 0)
            { // inside the medium
                N.negate();
                n = 1 / n;
            }
            n = 1 / n;
            float cost1 = -Vector3f::dot(N, ray.direction);
            float cost2 = 1.0 - n * n * (1.0 - cost1 * cost1);
            float Rprob = R0 + (1.0 - R0) * pow(1.0 - cost1, 5.0);
            if (cost2 > 0 && erand48(Xi) > Rprob)
            { // refraction direction
                ray.direction = ((ray.direction * n) + (N * (n * cost1 - sqrt(cost2)))).normalized();
            }
            else
            { // reflection direction
                ray.direction = (ray.direction + N * (cost1 * 2)).normalized();
            }
        }
    }
}

class PathTracer
{
public:
    const SceneParser &sceneParser;
    int samps;
    const char *fout;
    Vector3f (*Color_func)(Ray ray, const SceneParser &sceneParser, unsigned short *Xi); // 计算一个像素点的颜色
    PathTracer(const SceneParser &sceneParser, int samps, const char *method,
               const char *fout)
        : sceneParser(sceneParser), samps(samps), fout(fout)
    {
        if (!strcmp(method, "rc"))
            Color_func = rcColor;
        else if (!strcmp(method, "pt"))
            Color_func = ptColor;
        else
        {
            cout << "Unknown method: " << method << endl;
            exit(1);
        }
    }

    void render()
    {
        PerspectiveCamera *camera = (PerspectiveCamera *)sceneParser.getCamera();
        int w = camera->getWidth(), h = camera->getHeight();
        cout << "Width: " << w << " Height: " << h << endl;
        Image renderedImage(w, h);
        time_t start = time(NULL);
#pragma omp parallel for schedule(dynamic, 1) // OpenMP
        for (int x = 0; x < w; ++x)
        {
            float elapsed = (time(NULL) - start);
            float progress = (1. + x) / w;
            fprintf(stderr, "\rRendering (%d spp) %5.2f%% Time/Estimate: %.2f/%.2f sec",
                    samps, progress * 100., elapsed, elapsed / progress);
            unsigned short Xi[3] = {0, 0, (unsigned short)(x * x * x)};
            for (int y = 0; y < h; ++y)
            {
                Vector3f avg_color = Vector3f::ZERO;
                for (int s = 0; s < samps; ++s)
                {
                    Ray camRay = camera->generateRay(Vector2f(x + erand48(Xi) - 0.5, y + erand48(Xi) - 0.5));
                    avg_color += Color_func(camRay, sceneParser, Xi);
                }
                renderedImage.SetPixel(x, y, avg_color / samps);
            }
        }
        fprintf(stderr, "\n");
        renderedImage.SaveImage(fout);
    }
};

#endif // !PATH_TRACER_H