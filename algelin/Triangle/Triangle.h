#ifndef TRIANGLE_TRIANGLE_H_
#define TRIANGLE_TRIANGLE_H_

#include <Vec.h>
#include <cmath>

class Triangle {
public:
    Vec v0, v1, v2;
    Vec normal;

    //Constructor
    Triangle(Vec v0_, Vec v1_, Vec v2_, Vec color_,
            double ambient_coeff_, double diffuse_coeff_, double specular_coeff_,
            double transparent_coeff_, double specular_exponent_);

    bool intersect(const Ray &ray) const;
};
#endif // TRIANGLE_TRIANGLE_H