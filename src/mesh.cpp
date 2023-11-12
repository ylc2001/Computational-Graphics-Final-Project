#include "mesh.hpp"
#include <fstream>
#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <utility>
#include <sstream>

bool Mesh::intersect(const Ray &r, Hit &h, float tmin)
{
    // Optional: Change this brute force method into a faster one.
    if (aabb_intersect(r, tmin) == false)
        return false;
    bool result = false;
    for (int triId = 0; triId < (int)t.size(); ++triId)
    {
        // TriangleIndex &triIndex = t[triId];
        // Triangle triangle(v[triIndex[0]],
        //                   v[triIndex[1]], v[triIndex[2]], material);
        // triangle.normal = n[triId];
        result |= t_obj[triId].intersect(r, h, tmin);
    }
    return result;
}

bool Mesh::aabb_intersect(const Ray &r, float &t_min)
{
    Vector3f o(r.getOrigin());       // Origin of the ray
    Vector3f dir = r.getDirection(); // ray direction
    float txmin, txmax, tymin, tymax, tzmin, tzmax;
    if (abs(dir.x()) < 0.000001f && (o.x() < aabb_bounds[0].x() || o.x() > aabb_bounds[1].x()))
        return false;
    else
    {
        if (dir.x() >= 0)
        {
            txmin = (aabb_bounds[0].x() - o.x()) / dir.x();
            txmax = (aabb_bounds[1].x() - o.x()) / dir.x();
        }
        else
        {
            txmin = (aabb_bounds[1].x() - o.x()) / dir.x();
            txmax = (aabb_bounds[0].x() - o.x()) / dir.x();
        }
    }
    if (abs(dir.y()) < 0.000001f && (o.y() < aabb_bounds[0].y() || o.y() > aabb_bounds[1].y()))
        return false;
    else
    {
        if (dir.y() >= 0)
        {
            tymin = (aabb_bounds[0].y() - o.y()) / dir.y();
            tymax = (aabb_bounds[1].y() - o.y()) / dir.y();
        }
        else
        {
            tymin = (aabb_bounds[1].y() - o.y()) / dir.y();
            tymax = (aabb_bounds[0].y() - o.y()) / dir.y();
        }
    }
    if (abs(dir.z()) < 0.000001f && (o.z() < aabb_bounds[0].z() || o.z() > aabb_bounds[1].z()))
        return false;
    else
    {
        if (dir.z() >= 0)
        {
            tzmin = (aabb_bounds[0].z() - o.z()) / dir.z();
            tzmax = (aabb_bounds[1].z() - o.z()) / dir.z();
        }
        else
        {
            tzmin = (aabb_bounds[1].z() - o.z()) / dir.z();
            tzmax = (aabb_bounds[0].z() - o.z()) / dir.z();
        }
    }
    float t0, t1;
    //光线进入平面处（最靠近的平面）的最大t值
    t0 = max(tzmin, max(txmin, tymin));
    //光线离开平面处（最远离的平面）的最小t值
    t1 = min(tzmax, min(txmax, tymax));
    if (t0 < t1)
        return true;
    else
        return false;
}

void Mesh::updateBound(const Vector3f &vec)
{
    for (int i = 0; i < 3; ++i)
    {
        aabb_bounds[0][i] = aabb_bounds[0][i] < vec[i] ? aabb_bounds[0][i] : vec[i];
        aabb_bounds[1][i] = aabb_bounds[1][i] < vec[i] ? vec[i] : aabb_bounds[1][i];
    }
}

Mesh::Mesh(const char *filename, Material *material) : Object3D(material)
{

    // Optional: Use tiny obj loader to replace this simple one.
    std::ifstream f;
    f.open(filename);
    if (!f.is_open())
    {
        std::cout << "Cannot open " << filename << "\n";
        return;
    }
    std::string line;
    std::string vTok("v");
    std::string fTok("f");
    std::string texTok("vt");
    char bslash = '/', space = ' ';
    std::string tok;
    int texID;
    while (true)
    {
        std::getline(f, line);
        if (f.eof())
        {
            break;
        }
        if (line.size() < 3)
        {
            continue;
        }
        if (line.at(0) == '#')
        {
            continue;
        }
        std::stringstream ss(line);
        ss >> tok;
        if (tok == vTok)
        {
            Vector3f vec;
            ss >> vec[0] >> vec[1] >> vec[2];
            v.push_back(vec);
            updateBound(vec);
        }
        else if (tok == fTok)
        {
            if (line.find(bslash) != std::string::npos)
            {
                std::replace(line.begin(), line.end(), bslash, space);
                std::stringstream facess(line);
                TriangleIndex trig;
                facess >> tok;
                for (int ii = 0; ii < 3; ii++)
                {
                    facess >> trig[ii] >> texID;
                    trig[ii]--;
                }
                t.push_back(trig);
            }
            else
            {
                TriangleIndex trig;
                for (int ii = 0; ii < 3; ii++)
                {
                    ss >> trig[ii];
                    trig[ii]--;
                }
                t.push_back(trig);
            }
        }
        else if (tok == texTok)
        {
            Vector2f texcoord;
            ss >> texcoord[0];
            ss >> texcoord[1];
        }
    }
    computeNormal();

    f.close();

    Vector3f norm_avg;
    for (int vId = 0; vId < (int)v.size(); vId++)
    {
        norm_avg = Vector3f::ZERO;
        for (int triId = 0; triId < (int)t.size(); ++triId)
        {
            TriangleIndex &triIndex = t[triId];
            if (triIndex[0] == vId || triIndex[1] == vId || triIndex[2] == vId)
            {
                norm_avg += n[triId];
            }
        }
        norm_avg.normalize();
        v_normal.push_back(norm_avg);
    }

    for (int triId = 0; triId < (int)t.size(); ++triId)
    {
        TriangleIndex &triIndex = t[triId];
        Triangle triangle(v[triIndex[0]], v[triIndex[1]], v[triIndex[2]], material);
        triangle.normal = n[triId];
        triangle.set_v_normal(v_normal[triIndex[0]], v_normal[triIndex[1]], v_normal[triIndex[2]]);
        t_obj.push_back(triangle);
    }
}

void Mesh::computeNormal()
{
    n.resize(t.size());
    for (int triId = 0; triId < (int)t.size(); ++triId)
    {
        TriangleIndex &triIndex = t[triId];
        Vector3f a = v[triIndex[1]] - v[triIndex[0]];
        Vector3f b = v[triIndex[2]] - v[triIndex[0]];
        b = Vector3f::cross(a, b);
        n[triId] = b / b.length();
    }
}
