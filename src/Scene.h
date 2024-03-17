#ifndef SRC_SCENE_H
#define SRC_SCENE_H


Scene readSdlFile(const std::string& filename);
void readLights(const std::string& filename, std::vector<LightObject> triangles, Vec color, double intensity);
void readObjects(const std::string& filename, std::vector<Object> triangles, Vec color, 
                 double ambient_coeff_, double diffuse_coeff_, double specular_coeff_,
                 double transparent_coeff_, double specular_exponent_);

#endif // SRC_SCENE_H