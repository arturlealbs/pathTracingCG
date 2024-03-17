#include <iostream>
#include <fstream>
#include <string>
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

void writePPM(const std::string& filename, const Vec* image, int width, int height) {
    std::ofstream ppmFile(filename);
    if (!ppmFile.is_open()) {
        std::cerr << "Error: Failed to open file " << filename << " for writing." << std::endl;
        return;
    }

    // Write PPM header
    ppmFile << "P3\n";
    ppmFile << width << " " << height << "\n";
    ppmFile << "255\n"; // Maximum color value

    // Write pixel values
    for (int i = 0; i < width * height; i++){
            Vec pixelColor = image[i] * 255.0;
            ppmFile << static_cast<int>(pixelColor.x) << " "
                    << static_cast<int>(pixelColor.y) << " "
                    << static_cast<int>(pixelColor.z) << "\n";
        
    }

    ppmFile.close();
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

    //Defining Background Collor
    Vec background_color(0.0,0.0,0.0); //backgrounud values of the .sdl file

    //Ambient light intensity
    double ambient_intensity = 0.5; //ambient 0.5 of .sdl file

    //Execute path tracing algorithm
    Vec *image = pathTracing(camera, width, height, numSamples, background_color, ambient_intensity);
    
    //Write the result to a PPM File
    writePPM("output.ppm", image, width, height);

    delete[] image;
    return 0;
}