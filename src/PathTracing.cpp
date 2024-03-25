#include "PathTracing.h"

#include <Triangle.h>
#include <Ray.h>
#include <Vec.h>
#include <OrthoCam.h>
#include <Scene.h>
#include <vector>

Vec* pathTracing(Scene s){
    //The result of the PT should be saved in an array of pixels
    //std::vector<std::vector<Vec>> image(width, std::vector<Vec>(height, Vec()));
    Vec *image = new Vec(s.width * s.height); //usin this to not mix Vec and Vector, and is also more efficient

    for (int y = 0; y < s.height; ++y){
        for (int x = 0; x < s.width; ++x){
            Ray ray = s.camera.generateRay(x, y);
            Vec pixel_result(0.0, 0.0, 0.0);
            for (int i = 0; i < s.num_samples; ++i){
                tracePath(s, ray, pixel_result);
            }
            pixel_result = pixel_result + s.background_color * s.ambient_intensity;
            image[y * s.width + x] = pixel_result / s.num_samples;
        }
     }
    return image;
}

void tracePath(Scene s, Ray& ray, Vec pixel_result) {
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
    
    if (hitLight(ray, s.lights) != nullptr){}
    else if (hitSomething(ray, s.objects) != nullptr){}
    else{}
}

LightObject* hitLight(Ray& ray, std::vector<LightObject> lights){
    for(LightObject light: lights){
        if (light.intersect(ray)){
            return new LightObject(light);
        }
    }
    return nullptr;
}

Object* hitSomething(Ray& ray, std::vector<Object> objects){
    for(Object object: objects){
        if (object.intersect(ray)){
            return new Object(object);
        }
    }
    return nullptr;
}