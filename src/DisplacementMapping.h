#ifndef SRC_DISPLACEMENTMAPPING_H
#define SRC_DISPLACEMENTMAPPING_H

#include <Triangle.h>
#include <Vec.h>
#include <vector>
#include <PPM.h>

std::vector<Object> displacementMapping(std::pair<Object, Object> face, PPM image);

std::pair<double, double> getUVCoords(Vec vertex, double min_x, double max_x, double min_y, double max_y);

Vec getPixelValue(std::pair<double, double> uv, PPM image);

void vertexDisplacement(Vec vertex, Vec normal, Vec color);

std::vector<Object> subdivideMesh(std::vector<Object> triangles);

#endif // SRC_DISPLACEMENTMAPPING_H