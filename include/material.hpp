#ifndef MATERIAL_H
#define MATERIAL_H

#include <cassert>
#include "vecmath.h"
#include "stb_image.h"

#include "ray.hpp"
#include "hit.hpp"
#include "texture.h"
#include <iostream>
#include <math.h>
#include <cstring>

// --V--TODO: Implement Shade function that computes Phong introduced in class.
class Material
{
public:
    explicit Material(const Vector3f &color,
                      const Vector3f &s_color = Vector3f::ZERO, float s = 0,
                      const Vector3f &e_color = Vector3f::ZERO, float r = 1,
                      Vector3f t = Vector3f(1, 0, 0),
                      const char *textureFile = "")
        : color(color),
          specularColor(s_color),
          emission(e_color),
          shininess(s),
          refr(r),
          type(t)
    {
        // texture = nullptr;
        if (strlen(textureFile) > 0)
        {
            texture = stbi_load(textureFile, &w, &h, &c, 0);
            printf("Texture file: %s loaded. Size: %dx%dx%d\n", textureFile, w, h,
                   c);
        }
        else
        {
            texture = nullptr;
        }
    }

    virtual ~Material() = default;

    Vector3f getColor(int u, int v) const
    {
        if (!texture)
            return Vector3f::ZERO;
        u = u < 0 ? 0 : u;
        u = u > w - 1 ? w - 1 : u;
        v = v < 0 ? 0 : v;
        v = v > h - 1 ? h - 1 : v;
        int idx = v * w * c + u * c;
        return Vector3f(texture[idx + 0], texture[idx + 1], texture[idx + 2]) / 255.;
    }

    Vector3f getColor(float u, float v) const
    {
        if (!texture)
            return color;
        else
        {
            u -= int(u);
            v -= int(v);
            u = u < 0 ? 1 + u : u;
            v = v < 0 ? 1 + v : v;
            u = u * w;
            v = h * (1 - v);
            int iu = (int)u, iv = (int)v;
            float alpha = u - iu, beta = v - iv;
            Vector3f ret(0);
            ret += (1 - alpha) * (1 - beta) * getColor(iu, iv);
            ret += alpha * (1 - beta) * getColor(iu + 1, iv);
            ret += (1 - alpha) * beta * getColor(iu, iv + 1);
            ret += alpha * beta * getColor(iu + 1, iv + 1);
            return ret;
        }
    }

    Vector3f Shade(const Ray &ray, const Hit &hit,
                   const Vector3f &dirToLight, const Vector3f &lightColor)
    {
        Vector3f shaded = Vector3f::ZERO;
        Vector3f L = dirToLight;
        Vector3f N = hit.getNormal();
        Vector3f V = -ray.getDirection();
        Vector3f R = 2 * Vector3f::dot(N, L) * N - L;
        float diffuse = Vector3f::dot(L, N);
        if (diffuse < 0)
            diffuse = 0;
        float specular = Vector3f::dot(V, R);
        if (specular < 0)
            specular = 0;
        specular = pow(specular, shininess);
        Vector3f m_result = (color * diffuse) + (specularColor * specular);
        shaded.x() = lightColor.x() * m_result.x();
        shaded.y() = lightColor.y() * m_result.y();
        shaded.z() = lightColor.z() * m_result.z();
        return shaded;
    }

    Vector3f color;         // 颜色
    Vector3f specularColor; // 镜面反射系数
    Vector3f emission;      // 发光系数
    float shininess;        // 高光指数
    float refr;             // 折射率
    Vector3f type;          // 种类
    unsigned char *texture; // 纹理
    int w, h, c;            // 纹理相关变量
};

#endif // MATERIAL_H
