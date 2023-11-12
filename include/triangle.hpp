#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "object3d.hpp"
#include "vecmath.h"
#include <cmath>
#include <iostream>
using namespace std;

// --V--TODO: implement this class and add more fields as necessary,

class Triangle : public Object3D
{

public:
	Triangle() {};

	// a b c are three vertex positions of the triangle
	Triangle(const Vector3f &a, const Vector3f &b, const Vector3f &c, Material *m) : Object3D(m)
	{
		vertices[0] = a;
		vertices[1] = b;
		vertices[2] = c;
		E1 = vertices[1] - vertices[0];
		E2 = vertices[2] - vertices[0];
		normal = Vector3f::cross(E1, E2).normalized();
		normal_set = false;
	}

	Vector3f getNorm(const Vector3f& p)
	{
		if (normal_set == false)
			return normal;
		else
		{
			Vector3f va = (vertices[0] - p), vb = (vertices[1] - p), vc = (vertices[2] - p);
			float ra = Vector3f::cross(vb, vc).length(),
				  rb = Vector3f::cross(vc, va).length(),
				  rc = Vector3f::cross(va, vb).length();
			return (ra * norm_a + rb * norm_b + rc * norm_c).normalized();
		} 
	}

	bool intersect(const Ray &ray, Hit &hit, float tmin) override
	{
		Vector3f Ro = ray.getOrigin();
		Vector3f Rd = ray.getDirection();
		Vector3f norm = normal;
		float d = Vector3f::dot(norm, vertices[0]);
		float RHS = d - Vector3f::dot(Ro, norm);
		if (Vector3f::dot(Rd, norm) > 1e-9 || Vector3f::dot(Rd, norm) < -1e-9)
		{
			float t = RHS / Vector3f::dot(Rd, norm);
			if (t > tmin && t < hit.getT())
			{
				Vector3f hit_point = Ro + t * Rd; // judge if point is in triangle
				Vector3f judge1 = Vector3f::cross(E1, hit_point - vertices[0]);
				Vector3f judge2 = Vector3f::cross(vertices[2] - vertices[1], hit_point - vertices[1]);
				Vector3f judge3 = Vector3f::cross(vertices[0] - vertices[2], hit_point - vertices[2]);
				if (Vector3f::dot(judge1, judge2) >= 0 && Vector3f::dot(judge1, judge3) >= 0)
				{
					if (t <= tmin || t >= hit.getT())
						return false;
					else
					{
						Vector3f p = Ro + t * Rd;
						norm = getNorm(p);
						if (Vector3f::dot(norm, Rd) > 0)
							hit.set(t, material, -norm, material->color);
						else
							hit.set(t, material, norm, material->color);
						return true;
					}
				}
			}
		}
		return false;
	}
	
	void set_v_normal(const Vector3f an, const Vector3f &bn, const Vector3f &cn)
	{
		normal_set = true;
		norm_a = an;
		norm_b = bn;
		norm_c = cn;
	}
	Vector3f normal;				 // 面的法向量
	Vector3f norm_a, norm_b, norm_c; // 各个顶点的法向量(均值法)
	bool normal_set;				 // 当顶点法向有效, 为true
	Vector3f vertices[3];

protected:
	Vector3f E1, E2;
};

#endif // TRIANGLE_H
