#ifndef SRC_PATHTRACING_H
#define SRC_PATHTRACING_H

#include <Ray.h>
#include <Vec.h>
#include <OrthoCam.h>
#include <Scene.h>

Vec* pathTracing(Scene s);
void tracePath(Scene s, Ray& ray, Vec pixel_result);

#endif // SRC_PATHTRACING_H