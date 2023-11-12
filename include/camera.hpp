#ifndef CAMERA_H
#define CAMERA_H

#include "ray.hpp"
#include "vecmath.h"
#include "utils.hpp"
#include <float.h>
#include <cmath>

class Camera
{
public:
    Camera(const Vector3f &center, const Vector3f &direction, const Vector3f &up, int imgW, int imgH)
    {
        this->center = center;
        this->direction = direction.normalized();
        this->horizontal = Vector3f::cross(this->direction, up).normalized();
        this->up = Vector3f::cross(this->horizontal, this->direction);
        this->width = imgW;
        this->height = imgH;
    }

    // Generate rays for each screen-space coordinate
    virtual Ray generateRay(const Vector2f &point) = 0;
    virtual ~Camera() = default;

    int getWidth() const { return width; }
    int getHeight() const { return height; }

    // Extrinsic parameters
    Vector3f center;
    Vector3f direction;
    Vector3f up;
    Vector3f horizontal;

protected:
    // Intrinsic parameters
    int width;
    int height;
};

// --V--TODO: Implement Perspective camera
// You can add new functions or variables whenever needed.
class PerspectiveCamera : public Camera
{

public:
    PerspectiveCamera(const Vector3f &center, const Vector3f &direction,
                      const Vector3f &up, int imgW, int imgH, float angle,
                      float f = 10.0f, float aperture = 1.0f)
        : Camera(center, direction, up, imgW, imgH),
          focal(f),
          aperture(aperture)
    {
        // angle is in radian.
        this->angle = angle;
    }

    Ray generateRay(const Vector2f &point) override
    {
        // Fx = Fy
        int width = this->getWidth();
        int height = this->getHeight();
        float Fx = width / (2 * tan(angle * ((float)width / (float)height) / 2));
        float Fy = height / (2 * tan(angle / 2));
        float u = point.x();
        float v = point.y();
        float cx = width / 2.0;
        float cy = height / 2.0;
        Vector3f Ray_direction = ((u - cx) / Fx) * this->horizontal + ((v - cy) / Fy) * this->up + this->direction;
        Ray_direction.normalize();
        float t = focal / (Vector3f::dot(Ray_direction, this->direction));
        Vector3f Hit_point = center + Ray_direction * t;
        Vector3f random_center = center + ((erand48() - 0.5) * aperture) * up + ((erand48() - 0.5) * aperture) * horizontal;
        Vector3f random_dir = (Hit_point - random_center).normalized();
        return Ray(random_center, random_dir);
    }

    float angle;
    float focal, aperture;
};

#endif // CAMERA_H
