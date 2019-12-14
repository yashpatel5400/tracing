#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>

#include "tgaimage.h"
#include "model.h"
#include "our_gl.hpp"

const int width  = 1000;
const int height = 1000;

class GouraudShader : public IShader {
public:
    // written by vertex shader, read by fragment shader
    Model model;
    std::vector<TGAColor> varying_colors;
    Vec3f varying_intensity, light;
    Matrix viewportCoords, perspCoords, camCoords;

    GouraudShader(
        Model& model_, 
        Vec3f& light_,
        Vec3f& camera, 
        Vec3f& center
    ) : model(model_),
        light(light_) {
        camCoords = lookAt(camera, center, Vec3f(0,1,0));;
        perspCoords = perspectiveShift(camera.z);
        viewportCoords = viewport(width / 8, height / 8, width * 3/4, height * 3/4);
    }

    Vec3i vertex(int iface, int nthvert) {
        auto worldCoord = model.vert(iface, nthvert);
        auto normal = model.normal(iface, nthvert);
        auto uv = model.uv(iface, nthvert);
            
        // worldCoords.push_back(worldCoord);
        Vec3i screenCoord = fromHomogenous(viewportCoords * perspCoords * camCoords * toHomogenous(worldCoord));
        varying_colors.push_back(model.diffuse(uv));
        varying_intensity[nthvert] = light * normal;
        return screenCoord;
    }

    bool fragment(Vec3f bar) {
        float intensity = varying_intensity * bar;
        color = TGAColor(
            varying_colors[0][2] * bar[0] + varying_colors[1][2] * bar[1] + varying_colors[2][2] * bar[2],
            varying_colors[0][1] * bar[0] + varying_colors[1][1] * bar[1] + varying_colors[2][1] * bar[2],
            varying_colors[0][0] * bar[0] + varying_colors[1][0] * bar[1] + varying_colors[2][0] * bar[2],
            255
        ) * intensity;
        return false;
    }

    void clear() {
        varying_colors.clear();
    }
};

int main() {
    TGAImage image(width, height, TGAImage::RGB);

    Model model("african_head/african_head.obj");
    auto light = Vec3f(1,-1,1).normalize();
    Vec3f camera(1.0, 2.0, 3.0);
    Vec3f center(0, 0.0, 0);
    GouraudShader shader(model, light, camera, center);

    std::vector<float> zbuffer(width * height, -std::numeric_limits<float>::max());
    for (int f = 0; f < model.nfaces(); f++) {
        std::vector<Vec3i> screenCoords;
        for (int i = 0; i < 3; i++) {
            screenCoords.push_back(shader.vertex(f, i));
        }
        triangle(screenCoords, shader, zbuffer, image);
        shader.clear();
    }

    image.flip_vertically();
    image.write_tga_file("output.tga");
    return 0;
}
