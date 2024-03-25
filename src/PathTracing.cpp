#include "PathTracing.h"

#include <Ray.h>
#include <Vec.h>
#include <OrthoCam.h>

Vec* pathTracing(OrthoCam camera, int width, int height, int numSamples, Vec background_color, double ambient_intensity){
    //The result of the PT should be saved in an array of pixels
    //std::vector<std::vector<Vec>> image(width, std::vector<Vec>(height, Vec()));
    Vec *image = new Vec(width * height); //usin this to not mix Vec and Vector, and is also more efficient

    for (int y = 0; y < height; ++y){
        for (int x = 0; x < width; ++x){
            Ray ray = camera.generateRay(x, y);
            Vec pixel_result(0.0, 0.0, 0.0);
            for (int i = 0; i < numSamples; ++i){
                tracePath(ray, pixel_result);
            }
            pixel_result = pixel_result + background_color * ambient_intensity;
            image[y * width + x] = pixel_result / numSamples;
        }
     }
    return image;
}

void tracePath(Ray& ray, Vec pixel_result) {
    // Pseudo-code:
    // if (hit_light) {
    //     ray.result += ray.throughput * light;
    // } else if (hit_something) {
    //     BRDF brdf = SampleBRDF();
    //     ray.throughput *= brdf / brdf_pdf; // brdf_pdf is the probability density function for the sampled BRDF
    //     TracePath(brdf_ray); // brdf_ray is the new ray direction sampled from BRDF
    // } else {
    //     ray.result += ray.throughput * EnvMap();
    // }
}
