#include "Triangle.h"
#include <Vec.h>
#include <Ray.h>
#include <cmath>
#include <vector>

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
    std::vector<Triangle> subdivide() const {
        std::vector<Triangle> subTriangles;

        Vec mid01 = (v0 + v1) * 0.5; // Midpoint of edge v0-v1
        Vec mid12 = (v1 + v2) * 0.5; // Midpoint of edge v1-v2
        Vec mid20 = (v2 + v0) * 0.5; // Midpoint of edge v2-v0

        // Create four new triangles by connecting midpoints with original vertices
        subTriangles.push_back(Triangle(v0, mid01, mid20)); // Triangle 1
        subTriangles.push_back(Triangle(mid01, v1, mid12)); // Triangle 2
        subTriangles.push_back(Triangle(mid12, v2, mid20)); // Triangle 3
        subTriangles.push_back(Triangle(mid01, mid12, mid20)); // Triangle 4 (center)

        return subTriangles;
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