#ifndef MESH_H
#define MESH_H

#include <vector>
#include "object3d.hpp"
#include "triangle.hpp"
#include "vecmath.h"

class Mesh : public Object3D
{

public:
    Mesh(const char *filename, Material *m);

    struct TriangleIndex
    {
        TriangleIndex()
        {
            x[0] = 0;
            x[1] = 0;
            x[2] = 0;
        }
        int &operator[](const int i) { return x[i]; }
        // By Computer Graphics convention, counterclockwise winding is front face
        int x[3]{};
    };

    std::vector<Vector3f> v;        // vertices
    std::vector<Vector3f> v_normal; // 各个顶点的法向量
    std::vector<TriangleIndex> t;   // triangles TriangleIndex
    std::vector<Vector3f> n;        // 面的法向量, 顺序同 t, t_obj
    std::vector<Triangle> t_obj;    // triangle Object3D objects
    bool intersect(const Ray &r, Hit &h, float tmin) override;
    bool aabb_intersect(const Ray &r, float &t_min);
    void updateBound(const Vector3f &vec);
    Vector3f aabb_bounds[2]; // border. [0] is min, [1] is max

private:
    // Normal can be used for light estimation
    void computeNormal();
};

#endif
