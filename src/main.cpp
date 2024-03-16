#include <Ray.h>
#include <Vec.h>
#include <OrthoCam.h>

void TracePath(Ray& ray) {
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

int main(){
    //Setting up the Camera
    Vec eye(0.0, 0.0, 5.7); //Eye values of the sdl file
    double left = -1.0, bottom = -1.0, right = 1.0, top = 1.0; //Ortho values of the .sdl file

    OrthoCam camera(eye, left, right, bottom, top);

    //Defining Window Size
    int width = 200, height = 200;  //Size 200 200 in .sdl file
    
    //Defining Number of Samples
    int numSamples = 10; //npath 10 in .sdl size

    //
    Vec background_color(0.0,0.0,0.0);

    //
    double ambient_intensity = 0.5;

    //
    Vec *image = new Vec(width * height);
    //std::vector<std::vector<Vec>> image(width, std::vector<Vec>(height, Vec()));

     //PathTracing Pseudocode
     //TODO Implement all the necessary functions  for the PathTracing
     for (int y = 0; y < height; ++y){
        for (int x = 0; x < width; ++x){
            Ray ray = camera.generateRay(x, y);
            Vec pixel_result(0.0, 0.0, 0.0);
            for (int i = 0; i < numSamples; ++i){
                TracePath(ray);
            }
            pixel_result = pixel_result + background_color * ambient_intensity;
            image[y * width + x] = pixel_result / numSamples;
        }
     }


    return 0;
}