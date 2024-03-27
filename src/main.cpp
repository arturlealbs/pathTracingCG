#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <limits>
#include <random>

const int MAX_DEPTH = 5;
double bias = 1e-3;
//LINEAR ALGEBRA CLASSES
//Classes envolvidas com as operações feitas para cálculos de Path Tracing
class Vec {
public:
    double x, y, z;

    // Construtor
    Vec(double x_ = 0, double y_ = 0, double z_ = 0) : x(x_), y(y_), z(z_) {}

    //Operations
    double length() { return std::sqrt(x * x + y * y + z * z); };

    Vec operator-() const { return {-x, -y, -z}; };
    Vec operator+(const double &val) const {
        return {x + val, y + val, z + val};
    }
    Vec operator-(const double &val) const {
        return {x - val, y - val, z - val};
    }
    Vec operator*(const double &val) const {
        return {x * val, y * val, z * val};
    }

    Vec operator/(const double &val) const {
        return {x / val, y / val, z / val};
    }
    /* Vec operator+(const Vec &b) const {
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
 */
    Vec mult(const Vec &b) const {
        return Vec(x * b.x, y * b.y, z * b.z);
    }

    Vec div(const Vec &b) const{
        return Vec(x / b.x, y / b.y, z / b.z);
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
Vec operator+(const double &s, const Vec &v) { return v + s; };
Vec operator-(const double &s, const Vec &v) { return v - s; };
Vec operator*(const double &s, const Vec &v) { return v * s; };
Vec operator/(const double &s, const Vec &v) { return v / s; };

Vec operator+(const Vec &lhs, const Vec &rhs) {
  return {lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z};
};
Vec operator-(const Vec &lhs, const Vec &rhs) {
  return {lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z};
};
Vec operator*(const Vec &lhs, const Vec &rhs) {
  return {lhs.x * rhs.x, lhs.y * rhs.y, lhs.z * rhs.z};
};
Vec operator/(const Vec &lhs, const Vec &rhs) {
  return {lhs.x / rhs.x, lhs.y / rhs.y, lhs.z / rhs.z};
};

std::ostream &operator<<(std::ostream &os, const Vec &v) {
  os << "(" << v.x << ", " << v.y << ", " << v.z << ")";
  return os;
}



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
    std::pair<double,Vec> intersect(const Ray &ray) const {
        double t;
        const double EPSILON = 0.0000001;
        Vec edge1, edge2, h, s, q;
        double a, f, u, v;
        edge1 = v1 - v0;
        edge2 = v2 - v0;
        h = ray.direction.cross(edge2);
        a = edge1.dot(h);
        if (a > -EPSILON && a < EPSILON)
            return {-1, Vec()}; // Ray is parallel to this triangle.
        f = 1.0 / a;
        s = ray.origin - v0;
        u = f * s.dot(h);
        if (u < 0.0 || u > 1.0)
            return {-1, Vec()};
        q = s.cross(edge1);
        v = f * ray.direction.dot(q);
        if (v < 0.0 || u + v > 1.0)
            return {-1, Vec()};
        // Compute t to find out where the intersection point is on the line.
        t = f * edge2.dot(q);
        auto pontinho = ray.origin + t * ray.direction;
        //std::pair<t, pontinho>;
        //std::pair<-1, null>;
        if (t > EPSILON) // ray intersection
            return {t, pontinho};
        else // There is a line intersection but not a ray intersection.
            return {-1, Vec()};
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
    int hres, wres; 
    double height /* abs(y1 - y0) */, width/* abs(x1 - x0) */; 
    double pixel_h, pixel_w; 
    // Constructor
    OrthoCam(const Vec& eye_, double left_, double right_, double bottom_, double top_, int wres_, int hres_) :
        eye(eye_), left(left_), right(right_), bottom(bottom_), top(top_),wres(wres_), hres(hres_) {
            height = std::abs(top_ - bottom_); 
            width = std::abs( right_ - left_ );
            pixel_h = height/hres; 
            pixel_w = width/wres;
        };
    // (height / hres)
    // Generate a ray from the camera for a given pixel coordinates (x, y)'
    Ray generateRay(double x, double y, double errorW, double errorH) const {
        double u = left + (right - left) * (x + 0.5) / wres; // Calculate u coordinate
        double v = bottom + (top - bottom) * (y + 0.5) / hres; // Calculate v coordinate
        u = u + errorW;
        v = v + errorH;
        return Ray(eye, Vec(u, v, 0) - eye); // Generate ray from eye to (u, v, 0)
    }
};

//DISPLACEMENT CLASSES
//Classes e funções envolvidas na implementação do Displacement Mapping
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

struct Displacement{
    PPM texture;
    int face_id;

    Displacement(PPM texture_file, int face_id) : texture(texture_file), face_id(face_id) {};
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

//SCENE CLASSES
//Classes e funções envolvidas na construção da cena que será renderizada
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

void readLights(const std::string& filename, std::vector<LightObject>& triangles, Vec color, double intensity){
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

void readObjects(const std::string& filename, std::vector<Object>& triangles, Vec color, 
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
        if (token == "#" || token == ""){
            continue;
        } else if (token == "v") {
            Vec vertex;
            iss >> vertex.x >> vertex.y >> vertex.z;

            vertices.push_back(vertex);
        } else if (token == "f") {
            int v0_idx, v1_idx, v2_idx;
            iss >> v0_idx >> v1_idx >> v2_idx;
            Object triangle1(vertices[v0_idx - 1],vertices[v1_idx - 1],vertices[v2_idx - 1], 
            color, ambient_coeff_, diffuse_coeff_,  specular_coeff_, transparent_coeff_, specular_exponent_);

/*             iss >> v0_idx >> v1_idx >> v2_idx; 
            Object triangle2(vertices[v0_idx - 1], vertices[v1_idx - 1], vertices[v2_idx - 1], 
                             color, ambient_coeff_, diffuse_coeff_, specular_coeff_, transparent_coeff_, specular_exponent_);
             */
            triangles.push_back(triangle1);
        } else if (token == "displacement"){
            std::string texture_file;
            int face_id;
            iss>> texture_file >> face_id;
            Displacement disp = Displacement(PPM(texture_file), face_id);

            displacements.push_back(disp);
        }
    }
    /*
    int i = 1;
    for (Displacement disp : displacements){
        std::vector<Object> newTriangles = displacementMapping(faces[disp.face_id - i], disp.texture);
        for (Object triangle: newTriangles){
            triangles.push_back(triangle);
        }
        faces.erase(faces.begin() + (disp.face_id - i));
        ++i;
    }
    
    for (auto&& face: faces){
        triangles.push_back(face.first);
        triangles.push_back(face.second);
    }*/
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
    return Scene(OrthoCam(eye, left,right,bottom,top, width, height), width, height, background_color, 
                ambient_intensity, lights, objects, num_samples, tonemapping,seed, output);
}

//PATH TRACING CLASSES
//Classes para fazer os cálculos ópticos
Vec getRefraction(double n1, double n2, Vec i, Vec n) {
    float cosI = -i.dot(n);
    float sen2t = std::pow(n1 / n2, 2) * (1 - std::pow(cosI, 2));

    Vec t = ((i * (n1 / n2)) - n * std::sqrt(1 - sen2t)) + (((n1 / n2) * cosI));
    return t;
}

Vec difuse(Vec ip, double kd, Vec lightDir, Vec normal, Vec objColor) {
    Vec difuse;
    double prodEscalar = lightDir.dot(normal);
    double aux = kd * prodEscalar;
    difuse.x = ip.x * aux * objColor.x;
    difuse.y = ip.y * aux * objColor.y;
    difuse.z = ip.z * aux * objColor.z;

    return difuse;
};

std::pair <LightObject*, Vec> sampleLightSource(Scene &s, Vec inters){
        
    int randomLightFace = rand() %2;
    LightObject lightFace = s.lights[randomLightFace];
    double alpha = rand() % 100;
    double beta = rand() % 100;
    double gama = rand() % 100;
    double sum = alpha + beta + gama;
    alpha = alpha / sum;
    beta = beta / sum;
    gama = gama / sum;

    Vec v1 = lightFace.v0;
    Vec v2 = lightFace.v1;
    Vec v3 = lightFace.v2;

    Vec lightRand;

    lightRand.x = alpha * v1.x + beta * v2.x + gama * v3.x;
    lightRand.y = v1.y;
    lightRand.z = alpha * v1.z + beta * v2.z + gama * v3.z;

    //Vec toLight = lightRand - inters;
    //toLight = toLight.norm();
    return {&lightFace, lightRand};
};



//Classes e funções para calcular Path Tracing
std::pair<LightObject*, Vec> hitLight(Ray ray, std::vector<LightObject>& lights){
    for(auto&& light: lights){
        auto p = light.intersect(ray);
        if (p.first > 1){
            return {&light, p.second};
        }
    }
    return {nullptr, Vec()};
}

std::pair<Object*,Vec> hitSomething(Ray ray, std::vector<Object>& objects){
    for(auto&& object: objects){
        auto p = object.intersect(ray);
        if (p.first > 1){
            return {&object, p.second};
        }
    }
    return {nullptr, Vec()};
}

bool shadowRay(Scene& scene, Vec point) {
    //bool retorno = false;
    //LightObject* lightTarget = hitLight(ray, scene.lights).first;
    auto lightRand = sampleLightSource(scene, point).second;
    auto rayToLight = (lightRand - point).norm();
    auto p2 = hitSomething(Ray(point, rayToLight), scene.objects);
    Object* objectTarget = p2.first;
    return (objectTarget != nullptr);
}

Vec tracePath(Scene& s, Ray& ray, int depth) {
    // checkar se depth > MAX: return vec();
    if(depth > MAX_DEPTH) return Vec(); // return neutral color


    // Achar primeira colisão
    auto p1 = hitLight(ray, s.lights);
    LightObject* lightTarget = p1.first;
    auto p2 = hitSomething(ray, s.objects);
    Object* objectTarget = p2.first;

    if (lightTarget != nullptr){ 
        // Retorna cor da luz se nenhum objeto faz sombra
        if(objectTarget != nullptr){
            auto light_dist  = (ray.origin - p1.second).length();
            auto object_dist  = (ray.origin - p2.second).length();
            if(light_dist < object_dist) return lightTarget->color; 
        }else return lightTarget->color;
    }  
    if (objectTarget == nullptr)
        // Caso não tenha colisão retorna background
        return s.background_color;

    // propriedades do objeto que intersectou o raio
    auto closest = objectTarget; // objeto
    auto inters = p2.second; // ponto de intersecção
    auto normal = objectTarget -> normal; // normal da face do obj 
    // Fazendo sample da luz e pegando propriedades
    auto p3 = sampleLightSource(s, inters); 
    auto luz = p3.first; // light properties
    auto lightRand = p3.second; // light point 
    auto toLight = (lightRand - inters).norm(); // vector da colisão indo até direção a luz
    
    // propriedades phong
    float ka = closest->ambient_coeff, kd = closest->diffuse_coeff, 
          ks = closest->specular_coeff, kt = closest->transparent_coeff;

    // ---------------------------color calculation----------------------------
    // Rambiente = Ia*kar
    //float iA = s.ambient_intensity;
    Vec ambiente = s.ambient_intensity * ka * closest->color;

    // Rdifuso = Ip*kd(L.N)r
    Vec difusa = difuse(luz->color, kd, toLight, normal, closest->color);

    // Respecular = Ip*ks*(R.V)^n
    Vec rVetor = (2 * toLight.dot(normal) * normal) - toLight; rVetor = rVetor.norm();
    Vec vVetor = -1 * ray.direction; vVetor = vVetor.norm();
    float aux = luz->intensity * ks * pow(rVetor.dot(vVetor), closest->specular_exponent);
    Vec especular;
    especular = luz->color * aux;

    bool shadow = shadowRay(s, inters + bias * normal);

    // A cor direta no ponto 
    Vec corLocal = (shadow)? Vec(): ambiente+difusa+especular;
    return corLocal;
}

std::vector<Vec> pathTracing(Scene& s){
    //The result of the PT should be saved in an array of pixels
    //std::vector<std::vector<Vec>> image(width, std::vector<Vec>(height, Vec()));
    std::vector<Vec>image(s.height * s.width);//new Vec(s.width * s.height); //usin this to not mix Vec and Vector, and is also more efficient 
    std::random_device rd;
    std::mt19937 gen(rd()); // Mersenne Twister pseudo-random generator
    // Define the distribution for doubles
    std::uniform_real_distribution<double> errorH(-s.camera.pixel_h/2, s.camera.pixel_h/2); // Range [0.0, 1.0)
    std::uniform_real_distribution<double> errorW(-s.camera.pixel_w/2, s.camera.pixel_w/2);// pixel.Res/2); // Range [0.0, 1.0)
    for (int y = 0; y < s.height; ++y){
        for (int x = 0; x < s.width; ++x){
            Vec pixel_result(0.0, 0.0, 0.0);
            for (int i = 0; i < s.num_samples; ++i){
                Ray ray = s.camera.generateRay(x, y, errorW(gen), errorH(gen));
                //auto pixelSample = s.camera.[i][j] + {dist(gen), dist(gen), 0.0};
                //auto ray = Ray(s.camera.originGlobal, (pixelSample - s.camera.originGlobal).norm();
                //Ray ray = s.camera.generateRay(x, y, errorW(gen), errorH(gen));
                pixel_result = pixel_result + tracePath(s, ray, 0);
            }
            pixel_result = pixel_result /s.num_samples;
            //pixel_result = pixel_result.div(pixel_result + s.tonemapping); 
            image[(s.height - y - 1) * s.width + x] = pixel_result;
        }
    }
     
    return image;
}

//OUTPUT CLASSES
void writePNM(const std::string& filename, std::vector<Vec> image, int width, int height) {
    std::ofstream ppmFile(filename);
    if (!ppmFile.is_open()) {
        std::cerr << "Error: Failed to open file " << filename << " for writing." << std::endl;
        return;
    }

    ppmFile << "P3\n"; // P6 indicates binary encoding, Se não funcionar usa P3
    ppmFile << width << " " << height << "\n";
    ppmFile << "255\n"; // Maximum color value
    // Write pixel values in binary format
    for (int i = 0; i < width * height; i++) {
        Vec pixelColor = image[i];
        //std::cout << pixelColor.x << " " << pixelColor.y << " " << pixelColor.z << '\n';
        ppmFile << (((int)(pixelColor.x * 255))%256) << " "
                << (((int)(pixelColor.y * 255))%256) << " "
                << (((int)(pixelColor.z * 255))%256) << "\n";

    }
    ppmFile.close();
};

int main(){
    Scene scene = readSdlFile("../scenes/cornellroom.sdl");
    std::cout<<"finished reading file";
    std::vector<Vec> image = pathTracing(scene);
    std::cout<<"finished path tracing";
    //Write the result to a PPM File
    writePNM(scene.output, image, scene.width, scene.height);
    std::cout<<"finished writing";

    return 0;
}