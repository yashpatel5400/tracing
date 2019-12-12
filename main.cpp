#include <iostream>
#include <cmath>
#include <string>
#include <vector>

#include "tgaimage.h"
#include "model.h"

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red   = TGAColor(255, 0,   0,   255);

void line(int x0, int y0, int x1, int y1, TGAImage& image, TGAColor color) {
    float m = float(y1 - y0) / float(x1 - x0);

    // want to sample along the dimension that varies more
    bool isTransposed = false;
    if (std::abs(x1 - x0) < std::abs(y1 - y0)) { // relatively flat line
        std::swap(x0, y0);
        std::swap(x1, y1);
        isTransposed = true;
    }

    for (int x = x0; x <= x1; x++) {
        float alpha = float(x - x0) / (x1 - x0);
        int y = y0 * (1 - alpha) + y1 * alpha;
        if (isTransposed) {
            image.set(y, x, color);
        } else {
            image.set(x, y, color);
        }
    }
}

int main() {
    const int width = 500;
    const int height = 500;
    TGAImage image(width, height, TGAImage::RGB);
    Model model("african_head/african_head.obj");

    for (int f = 0; f < model.nfaces(); f++) {
        for (int v = 0; v < 3; v++) {
            auto v0 = model.vert(f, v);
            auto v1 = model.vert(f, (v + 1) % 3);

            int x0 = (v0.x + 1.0) * float(width)/2;
            int x1 = (v1.x + 1.0) * float(width)/2;
            int y0 = (v0.y + 1.0) * float(height)/2;
            int y1 = (v1.y + 1.0) * float(height)/2;
            line(x0, y0, x1, y1, image, white);
        }
    }

    image.flip_vertically();
    image.write_tga_file("output.tga");
    return 0;
}
