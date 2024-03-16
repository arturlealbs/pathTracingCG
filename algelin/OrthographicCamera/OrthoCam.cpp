#include <OrthoCam.h>

class OrthoCam {
public:
    Vec eye; // Position of the camera
    double left, right, bottom, top; // Orthographic projection window boundaries

    // Constructor
    OrthoCam(const Vec& eye_, double left_, double right_, double bottom_, double top_) :
        eye(eye_), left(left_), right(right_), bottom(bottom_), top(top_) {}

    // Generate a ray from the camera for a given pixel coordinates (x, y)
    Ray generateRay(double x, double y) const {
        double u = left + (right - left) * (x + 0.5); // Calculate u coordinate
        double v = bottom + (top - bottom) * (y + 0.5); // Calculate v coordinate
        return Ray(eye, Vec(u, v, 0) - eye); // Generate ray from eye to (u, v, 0)
    }
};