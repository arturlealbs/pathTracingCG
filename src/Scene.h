#ifndef SRC_SCENE_H
#define SRC_SCENE_H
#include <vector>
#include <Triangle.h>

struct Scene {
    OrthoCam camera;
    int width, height;
    Vec background_color;
    double ambient_intensity;
    std::vector<LightObject> lights;
    std::vector<Object> objects;
    int num_samples;
    double tonemapping;
    int seed;

    // Constructor
    Scene(const OrthoCam& cam, int w, int h, const Vec& bg_color, double ambient_int,
          const std::vector<LightObject>& l, const std::vector<Object>& obj,
          int num_samples, double tone_mapping, int s);

};

Scene readSdlFile(const std::string& filename);
void readLights(const std::string& filename, std::vector<LightObject> triangles, Vec color, double intensity);
void readObjects(const std::string& filename, std::vector<Object> triangles, Vec color, 
                 double ambient_coeff_, double diffuse_coeff_, double specular_coeff_,
                 double transparent_coeff_, double specular_exponent_);

#endif // SRC_SCENE_H