#ifndef SCENE_PARSER_H
#define SCENE_PARSER_H

#include <cassert>
#include "vecmath.h"
#include <vector>

class Camera;
class Light;
class DiskLight;
class Material;
class Object3D;
class Group;
class Sphere;
class Plane;
class Triangle;
class Transform;
class Mesh;
class Curve;
class RevSurface;

#define MAX_PARSER_TOKEN_LENGTH 1024

class SceneParser {
public:

    SceneParser() = delete;
    SceneParser(const char *filename);

    ~SceneParser();

    Camera *getCamera() const {
        return camera;
    }

    Vector3f getBackgroundColor() const {
        return background_color;
    }

    int getNumLights() const {
        return num_lights;
    }

    Light *getLight(int i) const {
        assert(i >= 0 && i < num_lights);
        return lights[i];
    }

    std::vector<DiskLight *> getDiskLights() const {
        std::vector<DiskLight *> dl;
        for (int i = 0; i < num_diskLights; ++i) dl.push_back(diskLights[i]);
        return dl;
    }

    int getNumMaterials() const {
        return num_materials;
    }

    Material *getMaterial(int i) const {
        assert(i >= 0 && i < num_materials);
        return materials[i];
    }

    Group *getGroup() const {
        return group;
    }

private:
    void parseFile();
    void parsePerspectiveCamera();
    void parseBackground();
    void parseLights();
    void parseDiskLights();
    DiskLight *parseDiskLight();
    Light *parsePointLight();
    Light *parseDirectionalLight();
    void parseMaterials();
    Material *parseMaterial();
    Object3D *parseObject(char token[MAX_PARSER_TOKEN_LENGTH]);
    Group *parseGroup();
    Sphere *parseSphere();
    Plane *parsePlane();
    Triangle *parseTriangle();
    Mesh *parseTriangleMesh();
    Transform *parseTransform();
    Curve *parseBezierCurve();
    Curve *parseBsplineCurve();
    RevSurface *parseRevSurface();

    int getToken(char token[MAX_PARSER_TOKEN_LENGTH]);

    Vector3f readVector3f();

    float readFloat();
    int readInt();

    FILE *file;
    Camera *camera;
    Vector3f background_color;
    int num_lights;
    Light **lights;
    int num_diskLights;
    DiskLight **diskLights;
    int num_materials;
    Material **materials;
    Material *current_material;
    Group *group;
};

#endif // SCENE_PARSER_H
