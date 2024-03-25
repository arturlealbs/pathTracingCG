#ifndef SRC_SCENE_H
#define SRC_SCENE_H

#include <Ray.h>
#include <Vec.h>
#include <OrthoCam.h>

Vec* pathTracing(OrthoCam camera, int width, int height, int numSamples, Vec background_color, double ambient_intensity);
void tracePath(Ray& ray, Vec pixel_result);

#endif // SRC_SCENE_H