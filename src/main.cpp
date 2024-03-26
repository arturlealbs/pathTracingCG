#include <iostream>
#include<vector>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <limits>

class Vec {
public:
    double x, y, z;

    // Construtor
    Vec(double x_ = 0, double y_ = 0, double z_ = 0) : x(x_), y(y_), z(z_) {}

    //Operations
    Vec operator+(const Vec &b) const {
        return Vec(x + b.x, y + b.y, z + b.z);
    }

    Vec operator-(const Vec &b) const {
        return Vec(x - b.x, y - b.y, z - b.z);
    }

    Vec operator*(double b) const {
        return Vec(x * b, y * b, z * b);
    }

    Vec operator/(double b) const {
        return Vec(x / b, y / b, z / b);
    }

    Vec mult(const Vec &b) const {
        return Vec(x * b.x, y * b.y, z * b.z);
    }

    Vec& norm() {
        double length = sqrt(x * x + y * y + z * z);
        return *this = *this * (1 / length);
    }

    double dot(const Vec &b) const {
        return x * b.x + y * b.y + z * b.z;
    }

    Vec cross(const Vec &b) const {
        return Vec(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x);
    }
};

class Ray {
public:
    Vec origin, direction;

    // Construtor
    Ray(Vec origin_, Vec direction_) : origin(origin_), direction(direction_) {}

};

class Triangle {
public:
    Vec v0, v1, v2;
    Vec normal;

    Triangle(Vec v0_, Vec v1_, Vec v2_) : v0(v0_), v1(v1_), v2(v2_) {
        // Calculate normal
        Vec e1 = v1 - v0;
        Vec e2 = v2 - v0;
        normal = e1.cross(e2).norm();
    }

    bool intersect(const Ray &ray) const {
        double t;
        const double EPSILON = 0.0000001;
        Vec edge1, edge2, h, s, q;
        double a, f, u, v;
        edge1 = v1 - v0;
        edge2 = v2 - v0;
        h = ray.direction.cross(edge2);
        a = edge1.dot(h);
        if (a > -EPSILON && a < EPSILON)
            return false; // Ray is parallel to this triangle.
        f = 1.0 / a;
        s = ray.origin - v0;
        u = f * s.dot(h);
        if (u < 0.0 || u > 1.0)
            return false;
        q = s.cross(edge1);
        v = f * ray.direction.dot(q);
        if (v < 0.0 || u + v > 1.0)
            return false;
        // Compute t to find out where the intersection point is on the line.
        t = f * edge2.dot(q);
        if (t > EPSILON) // ray intersection
            return true;
        else // There is a line intersection but not a ray intersection.
            return false;
    }
};

class Object : public Triangle {
public:
    Vec color;
    double ambient_coeff;
    double diffuse_coeff;
    double specular_coeff;
    double transparent_coeff;
    double specular_exponent;

    // Constructor
    Object(Vec v0_, Vec v1_, Vec v2_, Vec color_,
           double ambient_coeff_, double diffuse_coeff_,
           double specular_coeff_, double transparent_coeff_,
           double specular_exponent_)
        : Triangle(v0_, v1_, v2_), color(color_),
          ambient_coeff(ambient_coeff_), diffuse_coeff(diffuse_coeff_),
          specular_coeff(specular_coeff_), transparent_coeff(transparent_coeff_),
          specular_exponent(specular_exponent_) {}

    std::vector<Object> subdivide() const {
        std::vector<Object> subTriangles;

        Vec mid01 = (v0 + v1) * 0.5; // Midpoint of edge v0-v1
        Vec mid12 = (v1 + v2) * 0.5; // Midpoint of edge v1-v2
        Vec mid20 = (v2 + v0) * 0.5; // Midpoint of edge v2-v0

        
        // Create four new triangles by connecting midpoints with original vertices
        subTriangles.push_back(Object(v0, mid01, mid20, color, ambient_coeff, diffuse_coeff,specular_coeff,transparent_coeff,specular_coeff)); // Triangle 1
        subTriangles.push_back(Object(mid01, v1, mid12, color, ambient_coeff, diffuse_coeff,specular_coeff,transparent_coeff,specular_coeff)); // Triangle 2
        subTriangles.push_back(Object(mid12, v2, mid20, color, ambient_coeff, diffuse_coeff,specular_coeff,transparent_coeff,specular_coeff)); // Triangle 3
        subTriangles.push_back(Object(mid01, mid12, mid20, color, ambient_coeff, diffuse_coeff,specular_coeff,transparent_coeff,specular_coeff)); // Triangle 4 (center)

        return subTriangles;
    }
};

class LightObject : public Triangle {
public:
    Vec color;
    double intensity;

    // Constructor 
    LightObject(Vec v0_, Vec v1_, Vec v2_, Vec color_, double intensity_) 
        : Triangle(v0_, v1_, v2_), color(color_), intensity(intensity_) {
        Vec e1 = v1_ - v0_;
        Vec e2 = v2_ - v0_;
        normal = e1.cross(e2).norm();
    }
};

class OrthoCam {
public:
    Vec eye; // Position of the camera
    double left, right, bottom, top; // Orthographic projection window boundaries

    // Constructor
    OrthoCam(const Vec& eye_, double left_, double right_, double bottom_, double top_) :
        eye(eye_), left(left_), right(right_), bottom(bottom_), top(top_) {};

    // Generate a ray from the camera for a given pixel coordinates (x, y)
    Ray generateRay(double x, double y) const {
        double u = left + (right - left) * (x + 0.5); // Calculate u coordinate
        double v = bottom + (top - bottom) * (y + 0.5); // Calculate v coordinate
        return Ray(eye, Vec(u, v, 0) - eye); // Generate ray from eye to (u, v, 0)
    }
};

class PPM {
private:
    int width, height, maxColorValue;
    std::vector<Vec> pixels;

public:
    PPM(const std::string& filename) {
        std::ifstream file(filename, std::ios::binary);
        if (!file.is_open()) {
            std::cerr << "Error opening file: " << filename << std::endl;
            return;
        }

        std::string magicNumber;
        file >> magicNumber >> width >> height >> maxColorValue;

        if (magicNumber != "P6" || maxColorValue != 255) {
            std::cerr << "Unsupported PPM format or color depth." << std::endl;
            return;
        }

        file.ignore();

        pixels.resize(width * height);
        file.read(reinterpret_cast<char*>(pixels.data()), pixels.size() * sizeof(Vec));

        file.close();
    }

    int getWidth() const {
        return width;
    }

    int getHeight() const {
        return height;
    }

    Vec getPixel(int x, int y) const {
        return pixels[y * width + x];
    }
};

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

void writePNM(const std::string& filename, const Vec* image, int width, int height) {
    std::ofstream pnmFile(filename);
    if (!pnmFile.is_open()) {
        std::cerr << "Error: Failed to open file " << filename << " for writing." << std::endl;
        return;
    }

    pnmFile << "P6\n"; // P6 indicates binary encoding, Se não funcionar usa P3
    pnmFile << width << " " << height << "\n";
    pnmFile << "255\n"; // Maximum color value

    // Write pixel values in binary format
    for (int i = 0; i < width * height; i++) {
        Vec pixelColor = image[i] * 255.0;
        for (int i = 0; i < width * height; i++){
        Vec pixelColor = image[i] * 255.0;
        pnmFile << static_cast<int>(pixelColor.x) << " "
                << static_cast<int>(pixelColor.y) << " "
                << static_cast<int>(pixelColor.z) << "\n";
        
    }
    }

    pnmFile.close();
};

Vec getPixelValue(std::pair<double,double> uv, PPM image){
    int pixel_x = int(uv.first * image.getWidth());
    int pixel_y = int(uv.second * image.getWidth());

    return image.getPixel(pixel_x, pixel_y);
}

std::vector<Object> subdivideMesh(std::vector<Object> triangles){
    std::vector<Object> subdividedTriangles;
    for (Object triangle : triangles){
        std::vector<Object> temp = triangle.subdivide();
        for (Object subdividedTriangle : temp){
            subdividedTriangles.push_back(subdividedTriangle);
        }
    }
    return subdividedTriangles;
}

std::pair<double, double> getUVCoords(Vec vertex, double min_x, double max_x, double min_y, double max_y){
    double u_min = 0, u_max = 1;
    double v_min = 0, v_max = 1;

    double u = (vertex.x - min_x) / (max_x - min_x);
    double v = (vertex.y - min_y) / (max_y - min_y);

    u = (u - u_min) / (u_max - u_min);
    v = (v - v_min) / (v_max - v_min);

    return std::make_pair(u,v);
}

void vertexDisplacement(Vec vertex, Vec normal, Vec color){
    float height = (0.3 * color.x) + (0.59 * color.y) + (0.11*color.z);
    vertex = vertex + (normal * height * 0.5);
    
}

std::vector<Object> displacementMapping(std::pair<Object,Object> face,PPM image){
    //Extract from pair
    std::vector<Object> triangles;
    triangles.push_back(face.first);
 ;   triangles.push_back(face.second);

    //Getting the bounds of the UV Map
    double min_x = std::min_element(triangles.begin(), triangles.end(), [](const Triangle& a, const Triangle& b) {
        return std::min({a.v0.x, a.v1.x, a.v2.x}) < std::min({b.v0.x, b.v1.x, b.v2.x});
    })->v0.x;
    double max_x = std::max_element(triangles.begin(), triangles.end(), [](const Triangle& a, const Triangle& b) {
        return std::max({a.v0.x, a.v1.x, a.v2.x}) < std::max({b.v0.x, b.v1.x, b.v2.x});
    })->v0.x;
    double min_y = std::min_element(triangles.begin(), triangles.end(), [](const Triangle& a, const Triangle& b) {
        return std::min({a.v0.y, a.v1.y, a.v2.y}) < std::min({b.v0.y, b.v1.y, b.v2.y});
    })->v0.y;
    double max_y = std::max_element(triangles.begin(), triangles.end(), [](const Triangle& a, const Triangle& b) {
        return std::max({a.v0.y, a.v1.y, a.v2.y}) < std::max({b.v0.y, b.v1.y, b.v2.y});
    })->v0.y;

    //subdividing to match the amount of pixels in the height map
    int i = 1;
    while (i < image.getHeight()){
        triangles = subdivideMesh(triangles);
        i *= 2;
    }

    //Displacing the vertexes of each triangle
    for (Object triangle : triangles){
        vertexDisplacement(triangle.v0, triangle.normal, getPixelValue(getUVCoords(triangle.v0, min_x, max_x, min_y, max_y), image));
        vertexDisplacement(triangle.v1, triangle.normal, getPixelValue(getUVCoords(triangle.v1, min_x, max_x, min_y, max_y), image));
        vertexDisplacement(triangle.v2, triangle.normal, getPixelValue(getUVCoords(triangle.v2, min_x, max_x, min_y, max_y), image));
    }

    return triangles;
}

LightObject* hitLight(Ray& ray, std::vector<LightObject> lights){
    for(LightObject light: lights){
        if (light.intersect(ray)){
            return new LightObject(light);
        }
    }
    return nullptr;
}

Object* hitSomething(Ray& ray, std::vector<Object> objects){
    for(Object object: objects){
        if (object.intersect(ray)){
            return new Object(object);
        }
    }
    return nullptr;
}

Vec getRefraction(double n1, double n2, Vec i, Vec n) {
  float cosI = -i.dot(n);
  float sen2t = std::pow(n1 / n2, 2) * (1 - std::pow(cosI, 2));

  Vec t = (i * ((n1 / n2)) + (((n1 / n2) * cosI - n * std::sqrt(1 - sen2t))));
  return t;
}

void tracePath(Scene s, Ray& ray, Vec pixel_result) {
    // Pseudo-code:
    // if (hit_light) {
    //     ray.result += ray.throughput * light;
    // } else if (hit_something) {
    //     BRDF brdf = SampleBRDF();
    //     ray.throughput *= brdf / brdf_pdf; // brdf_pdf is the probability density function for the sampled BRDF
    //     TracePath(brdf_ray); // brdf_ray is the new ray direction sampled from BRDF
    // } else {
    //     ray.result += ray.throughput * EnvMap();
    // }
    if (hitLight(ray, s.lights) != nullptr){

    }
    else if (hitSomething(ray, s.objects) != nullptr){

    }
    else{

    }
}

Vec* pathTracing(Scene s){
    //The result of the PT should be saved in an array of pixels
    //std::vector<std::vector<Vec>> image(width, std::vector<Vec>(height, Vec()));
    Vec *image = new Vec(s.width * s.height); //usin this to not mix Vec and Vector, and is also more efficient
    for (int y = 0; y < s.height; ++y){
        for (int x = 0; x < s.width; ++x){
            Ray ray = s.camera.generateRay(x, y);
            Vec pixel_result(0.0, 0.0, 0.0);
            for (int i = 0; i < s.num_samples; ++i){
                tracePath(s, ray, pixel_result);
            }
            pixel_result = pixel_result + s.background_color * s.ambient_intensity;
            image[y * s.width + x] = pixel_result / s.num_samples;
        }
     }
    return image;
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

int main(){
    Scene scene = readSdlFile("../scenes/cornellroom.sdl");
    std::cout<<"finished reading file";
    Vec *image = pathTracing(scene);
    std::cout<<"finished path tracing";
    //Write the result to a PPM File
    writePNM(scene.output, image, scene.width, scene.height);
    std::cout<<"finished writing";

    delete[] image;
    return 0;
}