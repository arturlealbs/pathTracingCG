#ifndef SRC_PPM_H
#define SRC_PPM_H

#include <fstream>
#include <string>

#include <Vec.h>

class PPM {
private:
    int width, height, maxColorValue;
    std::vector<Vec> pixels;

public:
    PPM(const std::string& filename);

    int getWidth() const;

    int getHeight() const;

    Vec getPixel(int x, int y) const;
};
void writePNM(const std::string& filename, const Vec* image, int width, int height);

#endif // SRC_PPM_H