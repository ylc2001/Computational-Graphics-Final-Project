#ifndef REVSURFACE_HPP
#define REVSURFACE_HPP

#include "object3d.hpp"
#include "curve.hpp"
#include <tuple>

const int resolution = 10;
const int NEWTON_STEPS = 20;
const float NEWTON_EPS = 1e-4;

class RevSurface : public Object3D
{

    Curve *pCurve;

public:
    RevSurface(Curve *pCurve, Material *material) : Object3D(material), pCurve(pCurve)
    {
        // Check flat.
        for (const auto &cp : pCurve->getControls())
        {
            if (cp.z() != 0.0)
            {
                printf("Profile of revSurface must be flat on xy plane.\n");
                exit(0);
            }
        }
        aabb_bounds[0] = Vector3f(-pCurve->radius, pCurve->ymin - 3, -pCurve->radius);
        aabb_bounds[1] = Vector3f(pCurve->radius, pCurve->ymax + 3, pCurve->radius);
    }

    ~RevSurface() override
    {
        delete pCurve;
    }

    bool intersect(const Ray &r, Hit &h, float tmin) override
    {
        // (PA2 optional TODO): implement this for the ray-tracing routine using G-N iteration.
        float t, theta, mu;
        if (!aabb_intersect(r, t) || t > h.getT())
            return false;
        getUV(r, t, theta, mu);
        Vector3f normal, point;
        if (!newton(r, t, theta, mu, normal, point))
        {
            return false;
        }
        if (!isnormal(mu) || !isnormal(theta) || !isnormal(t))
            return false;
        if (t < 0 || mu < pCurve->range[0] || mu > pCurve->range[1] ||
            t > h.getT())
            return false;
        h.set(t, material, normal.normalized(),
              material->getColor(theta / (2 * M_PI), mu));
        return true;
    }

    bool newton(const Ray &r, float &t, float &theta, float &mu,
                Vector3f &normal, Vector3f &point)
    {
        Vector3f dmu, dtheta;
        for (int i = 0; i < NEWTON_STEPS; ++i)
        {
            if (theta < 0.0)
                theta += 2 * M_PI;
            if (theta >= 2 * M_PI)
                theta = fmod(theta, 2 * M_PI);
            if (mu >= 1)
                mu = 1.0 - FLT_EPSILON;
            if (mu <= 0)
                mu = FLT_EPSILON;
            point = getPoint(theta, mu, dtheta, dmu);
            Vector3f f = r.origin + r.direction * t - point;
            float dist2 = f.squaredLength();
            normal = Vector3f::cross(dmu, dtheta);
            if (dist2 < NEWTON_EPS)
                return true;
            float D = Vector3f::dot(r.direction, normal);
            t -= Vector3f::dot(dmu, Vector3f::cross(dtheta, f)) / D;
            mu -= Vector3f::dot(r.direction, Vector3f::cross(dtheta, f)) / D;
            theta += Vector3f::dot(r.direction, Vector3f::cross(dmu, f)) / D;
        }
        return false;
    }

    Vector3f getPoint(const float &theta, const float &mu, Vector3f &dtheta,
                      Vector3f &dmu)
    {
        Vector3f pt;
        Quat4f rot;
        rot.setAxisAngle(theta, Vector3f::UP);
        Matrix3f rotMat = Matrix3f::rotation(rot);
        CurvePoint cp = pCurve->getPoint(mu);
        pt = rotMat * cp.V;
        dmu = rotMat * cp.T;
        dtheta = Vector3f(-cp.V.x() * sin(theta), 0, -cp.V.x() * cos(theta));
        return pt;
    }

    void getUV(const Ray &r, const float &t, float &theta, float &mu)
    {
        Vector3f pt(r.origin + r.direction * t);
        theta = atan2(-pt.z(), pt.x()) + M_PI;
        mu = (pCurve->ymax - pt.y()) / (pCurve->ymax - pCurve->ymin);
    }

    bool aabb_intersect(const Ray &r, float &t_min)
    {
        Vector3f o(r.getOrigin()), invdir(1 / r.getDirection());
        vector<int> sgn = {invdir.x() < 0, invdir.y() < 0, invdir.z() < 0};
        t_min = INF;
        float tmin, tmax, tymin, tymax, tzmin, tzmax;
        tmin = (aabb_bounds[sgn[0]].x() - o.x()) * invdir.x();
        tmax = (aabb_bounds[1 - sgn[0]].x() - o.x()) * invdir.x();
        tymin = (aabb_bounds[sgn[1]].y() - o.y()) * invdir.y();
        tymax = (aabb_bounds[1 - sgn[1]].y() - o.y()) * invdir.y();
        if ((tmin > tymax) || (tymin > tmax))
            return false;
        if (tymin > tmin)
            tmin = tymin;
        if (tymax < tmax)
            tmax = tymax;
        tzmin = (aabb_bounds[sgn[2]].z() - o.z()) * invdir.z();
        tzmax = (aabb_bounds[1 - sgn[2]].z() - o.z()) * invdir.z();
        if ((tmin > tzmax) || (tzmin > tmax))
            return false;
        if (tzmin > tmin)
            tmin = tzmin;
        if (tzmax < tmax)
            tmax = tzmax;
        t_min = tmin;
        return true;
    }

    Vector3f aabb_bounds[2];
};

#endif // REVSURFACE_HPP
