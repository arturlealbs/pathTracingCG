#include <iostream>
#include <fstream>
#include <string>
#include <Ray.h>
#include <Vec.h>
#include <OrthoCam.h>
#include <Scene.h>
#include <PathTracing.h>


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
    Scene scene = readSdlFile("scenes\cornellroom.sdl");
  
    Vec *image = pathTracing(scene);
    
    //Write the result to a PPM File
    writePPM("output.ppm", image, scene.width, scene.height);

    delete[] image;
    return 0;
}