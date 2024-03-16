#ifndef RAY_RAY_H_
#define RAY_RAY_H_

#include <Vec.h>

class Ray {
public:
    Vec origin, direction;

    // Construtor
    Ray(Vec origin_, Vec direction_) : origin(origin_), direction(direction_) {}

};

#endif // RAY_RAY_H_