#ifndef TRIANGLE_TRIANGLE_H_
#define TRIANGLE_TRIANGLE_H_

#include <Vec.h>
#include <cmath>

class Triangle {
public:
    Vec v0, v1, v2;
    Vec normal;

    //Constructor
    Triangle(Vec v0_, Vec v1_, Vec v2_);

    bool intersect(const Ray &ray) const;
};

class Object : public Triangle{
public:
    Vec color;
    double ambient_coeff;
    double diffuse_coeff; 
    double specular_coeff;
    double transparent_coeff; 
    double specular_exponent;

    Object(Vec v0_, Vec v1_, Vec v2_, Vec color_, double ambient_coeff_, double diffuse_coeff_, 
           double specular_coeff_, double transparent_coeff_, double specular_exponent_);
    
    std::vector<Object> subdivide() const;
};

class LightObject : public Triangle {
public:
    Vec color;
    double intensity;

    // Constructor with normal calculation
    LightObject(Vec v0_, Vec v1_, Vec v2_, Vec color_, double intensity_);

};
#endif // TRIANGLE_TRIANGLE_H