#include "Scene.h"

#include <OrthoCam.h>
#include <Triangle.h>
#include <Vec.h>
#include <PPM.h>
#include <DisplacementMapping.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

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
    std::string output;

    // Constructor
    Scene(const OrthoCam& cam, int w, int h, const Vec& bg_color, double ambient_int,
          const std::vector<LightObject>& l, const std::vector<Object>& obj,
          int num_samples, double tonemapping, int seed, std::string output)
        : camera(cam), width(w), height(h), background_color(bg_color),
          ambient_intensity(ambient_int), lights(l), objects(obj),
          num_samples(num_samples), tonemapping(tonemapping), seed(seed), output(output) {}


};

struct Displacement{
    PPM texture;
    int face_id;

    Displacement(PPM texture_file, int face_id) : texture(texture_file), face_id(face_id) {};
};

Scene readSdlFile(const std::string& filename) {
    //opening .sdl file
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        exit(EXIT_FAILURE);
    }

    //declaring variables
    Vec eye;
    double left, bottom, right, top;
    int width, height;
    Vec background_color;
    double ambient_intensity;
    int num_samples;
    double tonemapping ;
    int seed;
    std::string output;

    std::vector<LightObject> lights;
    std::vector<Object> objects;

    //reading .sdf file lines
    std::string line;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string token;
        iss >> token;
        if (token == "eye") {
            double x, y, z;
            iss >> x >> y >> z;
            eye = Vec(x, y, z);
        } else if (token == "size") {
            iss >> width >> height;
        } else if (token == "ortho") {
            iss >> left >> bottom >> right >> top;
        } else if (token == "background") {
            double r, g, b;
            iss >> r >> g >> b;
            background_color = Vec(r, g, b);
        } else if (token == "ambient") {
            iss >> ambient_intensity;
        } else if (token == "light") {
            std::string name;
            double r, g, b, intensity;
            iss >> name >> r >> g >> b >> intensity;
            readLights(name, lights, Vec(r, g, b), intensity);
        } else if (token == "npaths") {
            iss >> num_samples;
        } else if (token == "tonemapping") {
            iss >> tonemapping;
        } else if (token == "seed") {
            iss >> seed;
        } else if (token == "object") {
            std::string name;
            double r, g, b, ka, kd, ks, kt, n;
            iss >> name >> r >> g >> b >> ka >> kd >> ks >> kt >> n;
            readObjects(name, objects, Vec(r, g, b), ka, kd, ks, kt, n);
        } else if (token == "output"){
            iss >> output;
        }
    }

    file.close();

    //creating a scene and returning it
    return Scene(OrthoCam(eye, left,right,bottom,top), width, height, background_color, 
                ambient_intensity, lights, objects, num_samples, tonemapping,seed, output);
}

void readLights(const std::string& filename, std::vector<LightObject> triangles, Vec color, double intensity){
    std::string filepath = "../files/" + filename;
    std::ifstream file(filepath);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        exit(1);
    }

    std::vector<Vec> vertices;
    std::string line;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string token;
        iss >> token;
        if (token == "v") {
            Vec vertex;
            iss >> vertex.x >> vertex.y >> vertex.z;
            vertices.push_back(vertex);
        } else if (token == "f") {
            int v0_idx, v1_idx, v2_idx;
            iss >> v0_idx >> v1_idx >> v2_idx;
            LightObject triangle(vertices[v0_idx - 1],vertices[v1_idx - 1],vertices[v2_idx - 1], 
            color, intensity);
     
            triangles.push_back(triangle);
        }
    }

    file.close();
}

void readObjects(const std::string& filename, std::vector<Object> triangles, Vec color, 
                 double ambient_coeff_, double diffuse_coeff_, double specular_coeff_,
                 double transparent_coeff_, double specular_exponent_){
    std::string filepath = "../files/" + filename;
    std::ifstream file(filepath);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        exit(1);
    }

    std::vector<Vec> vertices;
    std::vector<std::pair<Object,Object>> faces;
    std::string line;
    std::vector<Displacement> displacements;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string token;
        iss >> token;

        if (token == "v") {
            Vec vertex;
            iss >> vertex.x >> vertex.y >> vertex.z;

            vertices.push_back(vertex);
        } else if (token == "f") {
            int v0_idx, v1_idx, v2_idx;
            iss >> v0_idx >> v1_idx >> v2_idx;
            Object triangle1(vertices[v0_idx - 1],vertices[v1_idx - 1],vertices[v2_idx - 1], 
            color, ambient_coeff_, diffuse_coeff_,  specular_coeff_, transparent_coeff_, specular_exponent_);

            iss >> v0_idx >> v1_idx >> v2_idx; 
            Object triangle2(vertices[v0_idx - 1], vertices[v1_idx - 1], vertices[v2_idx - 1], 
                             color, ambient_coeff_, diffuse_coeff_, specular_coeff_, transparent_coeff_, specular_exponent_);
            
            faces.push_back({triangle1,triangle2});
        } else if (token == "displacement"){
            std::string texture_file;
            int face_id;
            iss>> texture_file >> face_id;
            Displacement disp = Displacement(PPM(texture_file), face_id);

            displacements.push_back(disp);
        }
    }
    int i = 1;
    for (Displacement disp : displacements){
        std::vector<Object> newTriangles = displacementMapping(faces[disp.face_id - i], disp.texture);
        for (Object triangle: newTriangles){
            triangles.push_back(triangle);
        }
        faces.erase(faces.begin() + (disp.face_id - i));
        ++i;
    }
    
    for (std::pair<Object,Object> face: faces){
        triangles.push_back(face.first);
        triangles.push_back(face.second);
    }
    file.close();
}