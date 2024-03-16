#ifndef VECTOR_VEC_H_
#define VECTOR_VEC_H_

#include <cmath>

class Vec {
public:
    double x, y, z;

    // Construtor
    Vec(double x_ = 0, double y_ = 0, double z_ = 0);

    // Operator
    Vec operator+(const Vec &b) const;
    Vec operator-(const Vec &b) const;
    Vec operator*(double b) const;
    Vec mult(const Vec &b) const;
    Vec& norm();
    double dot(const Vec &b) const;
    Vec cross(const Vec &b) const;
};

#endif // VECTOR_VEC_H