#include <Triangle.h>
#include <Vec.h>
#include <vector>
#include <PPM.h>

Vec displacementMapping(std::pair<Object,Object> face,PPM image){
    //Extract from pair
    std::vector<Object> triangles;
    triangles.push_back(face.first);
    triangles.push_back(face.second);

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

    
}


