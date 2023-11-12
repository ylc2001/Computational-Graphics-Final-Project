#ifndef GROUP_H
#define GROUP_H

#include "object3d.hpp"
#include "ray.hpp"
#include "hit.hpp"
#include <iostream>
#include <vector>

// --V--TODO: Implement Group - add data structure to store a list of Object*
class Group : public Object3D
{

public:
    Group()
    {
        this->num_objects = 0;
    }

    explicit Group(int num_objects)
    {
        this->num_objects = num_objects;
        objects.clear();
    }

    ~Group() override
    {
    }

    bool intersect(const Ray &r, Hit &h, float tmin) override
    {
        bool result = false;
        float current_t;
        float min_t = 10000000000;
        int min_index = -1;
        for (int index = 0; index < (int)objects.size(); index++)
        {
            bool is_intersect = objects[index]->intersect(r, h, tmin);
            if (is_intersect)
            {
                result = true;
                current_t = h.getT();
                if (current_t < min_t)
                {
                    min_t = current_t;
                    min_index = index;
                }
            }
        }
        if (result)
        {
            objects[min_index]->intersect(r, h, tmin);
        }
        return result;
    }

    void addObject(int index, Object3D *obj)
    {
        objects.push_back(obj);
        num_objects++;
    }

    int getGroupSize()
    {
        return num_objects;
    }

private:
    int num_objects;
    std::vector<Object3D *> objects;
};

#endif
