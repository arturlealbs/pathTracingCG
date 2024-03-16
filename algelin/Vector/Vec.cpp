#include "Vec.h"

#include <cmath>

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