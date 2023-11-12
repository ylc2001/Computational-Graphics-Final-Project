#ifndef SPHERE_H
#define SPHERE_H

#include "object3d.hpp"
#include "vecmath.h"
#include <cmath>

// --V--TODO: Implement functions and add more fields as necessary

class Sphere : public Object3D
{
public:
    Sphere()
    {
        // unit ball at the center
        _radius = 1;
        _center = Vector3f(0.0, 0.0, 0.0);
    }

    Sphere(const Vector3f &center, float radius, Material *material)
        : Object3D(material), _center(center), _radius(radius)
    {
        //
    }

    ~Sphere() override = default;

    bool intersect(const Ray &r, Hit &h, float tmin) override
    {
        //
        Vector3f Ro = r.getOrigin();
        Vector3f Rd = r.getDirection();

        float a = Rd.squaredLength();
        float b = 2 * (Rd.x() * (Ro.x() - _center.x()) + Rd.y() * (Ro.y() - _center.y()) + Rd.z() * (Ro.z() - _center.z()));
        float c = (Ro - _center).squaredLength() - _radius * _radius;

        float delta = b * b - 4 * a * c;
        if (delta >= 0)
        {
            float t = (-b - sqrt(delta)) / (2 * a);
            if (t > tmin && t < h.getT())
            {
                Vector3f normal = (Ro + t * Rd) - _center;
                normal.normalize();
                Vector3f OPvec(Ro + Rd * t - _center);
                Vector3f OPnormal = OPvec.normalized();
                float u = 0.5 + atan2(OPnormal.x(), OPnormal.z()) / (2 * M_PI);
                float v = -0.5 + asin(OPnormal.y()) / M_PI;
                h.set(t, material, normal, material->getColor(u, v));
                return true;
            }
        }
        return false;
    }

protected:
    Vector3f _center;
    float _radius;
};

#endif
