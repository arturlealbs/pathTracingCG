#include <iostream>
#include <Ray.h>
#include <Vec.h>
#include <OrthoCam.h>
#include <Scene.h>
#include <PathTracing.h>
#include <PPM.h>


int main(){
    Scene scene = readSdlFile("../scenes/cornellroom.sdl");
  
    Vec *image = pathTracing(scene);
    
    //Write the result to a PPM File
    writePNM(scene.output, image, scene.width, scene.height);

    delete[] image;
    return 0;
}