#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>

#include "tgaimage.h"
#include "model.h"

const TGAColor white  = TGAColor(255, 255, 255, 255);
const TGAColor red    = TGAColor(255,   0,   0, 255);
const TGAColor green  = TGAColor(  0, 255,   0, 255);
const TGAColor blue   = TGAColor(  0,   0, 255, 255);

void line(Vec2i v0, Vec2i v1, TGAImage& image, TGAColor color) {    
    int x0 = v0.x;
    int y0 = v0.y;
    int x1 = v1.x;
    int y1 = v1.y;

    // want to sample along the dimension that varies more
    bool isTransposed = false;
    if (std::abs(x1 - x0) < std::abs(y1 - y0)) { // relatively flat line
        std::swap(x0, y0);
        std::swap(x1, y1);
        isTransposed = true;
    }

    if (x0 > x1) {
        std::swap(x0, x1);
        std::swap(y0, y1);
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

void triangle(Vec2i v0, Vec2i v1, Vec2i v2, TGAImage& image, TGAColor color) {
    float m0 = float(v1.y - v0.y) / float(v1.x - v0.x);
    float m1 = float(v2.y - v1.y) / float(v2.x - v1.x);
    float m2 = float(v0.y - v2.y) / float(v0.x - v2.x);

    std::vector<int> xs = {v0.x, v1.x, v2.x};
    std::vector<int> ys = {v0.y, v1.y, v2.y};
    std::vector<float> ms = {m0, m1, m2};

    int minX = std::min(std::min(xs[0], xs[1]), xs[2]);
    int maxX = std::max(std::max(xs[0], xs[1]), xs[2]);

    int minY = std::min(std::min(ys[0], ys[1]), ys[2]);
    int maxY = std::max(std::max(ys[0], ys[1]), ys[2]);

    bool onLine = false;
    Vec2i startPoint, endPoint;

    for (int x = minX; x <= maxX; x++) {
        for (int y = minY; y <= maxY; y++) {
            for (int i = 0; i < 3; i++) {
                int pred = ms[i] * (x - xs[i]) + ys[i];
                if (y == pred) {
                    if (onLine) {
                        endPoint = Vec2i(x, y);
                        onLine = false;
                        line(startPoint, endPoint, image, color);
                    } else {
                        startPoint = Vec2i(x, y);
                        onLine = true;
                    }
                }
            }
        }
    }

    // std::vector<Vec2i> vertices = { v0, v1, v2 };
    // std::sort(vertices.begin(), vertices.end(), [](const Vec2i& a, const Vec2i& b) {
    //     return a.y < b.y;
    // });
 
    // line(v0, v1, image, color);
    // line(v1, v2, image, color);
    // line(v2, v0, image, color);
}

int main() {
    const int width = 200;
    const int height = 200;
    TGAImage image(width, height, TGAImage::RGB);
    Vec2i t0[3] = {Vec2i(10, 70),   Vec2i(50, 160),  Vec2i(70, 80)}; 
    Vec2i t1[3] = {Vec2i(180, 50),  Vec2i(150, 1),   Vec2i(70, 180)}; 
    Vec2i t2[3] = {Vec2i(180, 150), Vec2i(120, 160), Vec2i(130, 180)}; 
    triangle(t0[0], t0[1], t0[2], image, red); 
    triangle(t1[0], t1[1], t1[2], image, white); 
    triangle(t2[0], t2[1], t2[2], image, green);

    // Model model("african_head/african_head.obj");

    // for (int f = 0; f < model.nfaces(); f++) {
    //     // for (int v = 0; v < 3; v++) {
            
    //     //     // int x0 = (v0.x + 1.0) * float(width)/2;
    //     //     // int x1 = (v1.x + 1.0) * float(width)/2;
    //     //     // int y0 = (v0.y + 1.0) * float(height)/2;
    //     //     // int y1 = (v1.y + 1.0) * float(height)/2;
    //     //     // line(x0, y0, x1, y1, image, white);
    //     // }
    //     auto v0 = model.vert(f, 0);
    //     auto v1 = model.vert(f, 1);
    //     auto v2 = model.vert(f, 2);
            
    //     triangle(v0, v1, v2, image, white);
    // }

    image.flip_vertically();
    image.write_tga_file("output.tga");
    return 0;
}
