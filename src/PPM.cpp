#include "PPM.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <Vec.h>

class PPM {
private:
    int width, height, maxColorValue;
    std::vector<Vec> pixels;

public:
    PPM(const std::string& filename) {
        std::ifstream file(filename, std::ios::binary);
        if (!file.is_open()) {
            std::cerr << "Error opening file: " << filename << std::endl;
            return;
        }

        std::string magicNumber;
        file >> magicNumber >> width >> height >> maxColorValue;

        if (magicNumber != "P6" || maxColorValue != 255) {
            std::cerr << "Unsupported PPM format or color depth." << std::endl;
            return;
        }

        file.ignore();

        pixels.resize(width * height);
        file.read(reinterpret_cast<char*>(pixels.data()), pixels.size() * sizeof(Vec));

        file.close();
    }

    int getWidth() const {
        return width;
    }

    int getHeight() const {
        return height;
    }

    Vec getPixel(int x, int y) const {
        return pixels[y * width + x];
    }
};

void writePNM(const std::string& filename, const Vec* image, int width, int height) {
    std::ofstream pnmFile(filename);
    if (!pnmFile.is_open()) {
        std::cerr << "Error: Failed to open file " << filename << " for writing." << std::endl;
        return;
    }

    pnmFile << "P6\n"; // P6 indicates binary encoding, Se não funcionar usa P3
    pnmFile << width << " " << height << "\n";
    pnmFile << "255\n"; // Maximum color value

    // Write pixel values in binary format
    for (int i = 0; i < width * height; i++) {
        Vec pixelColor = image[i] * 255.0;
        for (int i = 0; i < width * height; i++){
        Vec pixelColor = image[i] * 255.0;
        pnmFile << static_cast<int>(pixelColor.x) << " "
                << static_cast<int>(pixelColor.y) << " "
                << static_cast<int>(pixelColor.z) << "\n";
        
    }
    }

    pnmFile.close();
};