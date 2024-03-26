#include "DisplacementMapping.h"

#include <Triangle.h>
#include <Vec.h>
#include <vector>
#include <PPM.h>

std::vector<Object> displacementMapping(std::pair<Object,Object> face,PPM image){
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

std::pair<double, double> getUVCoords(Vec vertex, double min_x, double max_x, double min_y, double max_y){
    double u_min = 0, u_max = 1;
    double v_min = 0, v_max = 1;

    double u = (vertex.x - min_x) / (max_x - min_x);
    double v = (vertex.y - min_y) / (max_y - min_y);

    u = (u - u_min) / (u_max - u_min);
    v = (v - v_min) / (v_max - v_min);

    return std::make_pair(u,v);
}

Vec getPixelValue(std::pair<double,double> uv, PPM image){
    int pixel_x = int(uv.first * image.getWidth());
    int pixel_y = int(uv.second * image.getWidth());

    return image.getPixel(pixel_x, pixel_y);
}

void vertexDisplacement(Vec vertex, Vec normal, Vec color){
    float height = (0.3 * color.x) + (0.59 * color.y) + (0.11*color.z);
    vertex = vertex + (normal * height * 0.5);
    
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