#include "Triangle.h"
#include <Vec.h>
#include <Ray.h>
#include <cmath>

class Triangle {
public:
    Vec v0, v1, v2;
    Vec normal;
    Vec color; // RGB values
    double ambient_coeff;
    double diffuse_coeff;
    double specular_coeff;
    double transparent_coeff;
    double specular_exponent;

    Triangle(Vec v0_, Vec v1_, Vec v2_, Vec color_,
            double ambient_coeff_, double diffuse_coeff_, double specular_coeff_,
            double transparent_coeff_, double specular_exponent_)
        : v0(v0_), v1(v1_), v2(v2_), color(color_),
          ambient_coeff(ambient_coeff_), diffuse_coeff(diffuse_coeff_),
          specular_coeff(specular_coeff_), transparent_coeff(transparent_coeff_),
          specular_exponent(specular_exponent_) {
        // Calculate normal
        Vec e1 = v1 - v0;
        Vec e2 = v2 - v0;
        normal = e1.cross(e2).norm();
    }


    // Ray-triangle intersection test
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
            return false; // This ray is parallel to this triangle.
        f = 1.0 / a;
        s = ray.origin - v0;
        u = f * s.dot(h);
        if (u < 0.0 || u > 1.0)
            return false;
        q = s.cross(edge1);
        v = f * ray.direction.dot(q);
        if (v < 0.0 || u + v > 1.0)
            return false;
        //Compute t to find out where the intersection point is on the line.
        t = f * edge2.dot(q);
        if (t > EPSILON) // ray intersection
            return true;
        else // There is a line intersection but not a ray intersection.
            return false;
    }
};