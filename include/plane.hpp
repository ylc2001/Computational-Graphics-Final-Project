#ifndef PLANE_H
#define PLANE_H

#include "object3d.hpp"
#include "vecmath.h"
#include <cmath>

// --V--TODO: Implement Plane representing an infinite plane
// function: ax+by+cz=d
// choose your representation , add more fields and fill in the functions

class Plane : public Object3D
{
public:
    Plane() : n(Vector3f())
    {
        D = 0;
    }

    Plane(const Vector3f &normal, float d, Material *m) : Object3D(m), n(normal)
    {
        D = d;
    }

    ~Plane() override = default;

    bool intersect(const Ray &r, Hit &h, float tmin) override
    {
        Vector3f Ro = r.getOrigin();
        Vector3f Rd = r.getDirection();
        if (Vector3f::dot(n, Rd) > 1e-9 || Vector3f::dot(n, Rd) <-1e-9)
        {
            float t = (D - Vector3f::dot(n, Ro)) / (Vector3f::dot(n, Rd));
            if (t <= tmin || t >= h.getT())
                return false;
            else
            {
                if (Vector3f::dot(n, Rd) > 0)
                    h.set(t, material, -n.normalized(), material->color);
                else
                    h.set(t, material, n.normalized(),  material->color);
                return true;
            }
        }
       return false;
    }

protected:
    Vector3f n;
    float D;
};

#endif // PLANE_H
