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

void triangleLineSweep(Vec2i v0, Vec2i v1, Vec2i v2, TGAImage& image, TGAColor color) {
    if (v0.y > v1.y) { std::swap(v0, v1); }
    if (v1.y > v2.y) { std::swap(v1, v2); }
    if (v0.y > v1.y) { std::swap(v0, v1); }

    for (int y = v0.y; y <= v1.y; y++) {
        float topAlpha = float(y - v0.y) / (v2.y - v0.y);
        float bottomAlpha = float(y - v0.y) / (v1.y - v0.y);

        auto pointA = std::round(topAlpha * v2.x + (1 - topAlpha) * v0.x);
        auto pointB = std::round(bottomAlpha * v1.x + (1 - bottomAlpha) * v0.x);

        auto start = std::min(pointA, pointB);
        auto end = std::max(pointA, pointB);

        for (int x = start; x <= end; x++) {
            image.set(x, y, color);
        }
    }

    for (int y = v1.y + 1; y <= v2.y; y++) {
        float topAlpha = float(y - v0.y) / (v2.y - v0.y);
        float bottomAlpha = float(y - v1.y) / (v2.y - v1.y);

        auto pointA = std::round(topAlpha * v2.x + (1 - topAlpha) * v0.x);
        auto pointB = std::round(bottomAlpha * v2.x + (1 - bottomAlpha) * v1.x);

        auto start = std::min(pointA, pointB);
        auto end = std::max(pointA, pointB);

        for (int x = start; x <= end; x++) {
            image.set(x, y, color);
        }
    }
}

void triangle(Vec2i v0, Vec2i v1, Vec2i v2, TGAImage& image, TGAColor color) {
    std::vector<int> xs = {v0.x, v1.x, v2.x};
    std::vector<int> ys = {v0.y, v1.y, v2.y};

    int minX = std::min(std::min(xs[0], xs[1]), xs[2]);
    int maxX = std::max(std::max(xs[0], xs[1]), xs[2]);
    int minY = std::min(std::min(ys[0], ys[1]), ys[2]);
    int maxY = std::max(std::max(ys[0], ys[1]), ys[2]);

    auto AB = v1 - v0;
    auto AC = v2 - v0;

    for (int x = minX; x <= maxX; x++) {
        for (int y = minY; y <= maxY; y++) {
            Vec2i p(x, y);
            auto PA = v0 - p;
            Vec3f first(AB.x, AC.x, PA.x);
            Vec3f second(AB.y, AC.y, PA.y);
            
            Vec3f res = cross(first, second);
            bool foundNeg = false;

            if (res.z != 0) {
                auto u = res.x / res.z;
                auto v = res.y / res.z;
                auto t = 1 - u - v;

                if (u >= 0 && v >= 0 && t >= 0) {
                    image.set(x, y, color);
                }
            }
        }   
    }
}

int main() {
    const int width = 500;
    const int height = 500;
    TGAImage image(width, height, TGAImage::RGB);
    // Vec2i t0[3] = {Vec2i(10, 70),   Vec2i(50, 160),  Vec2i(70, 80)}; 
    // Vec2i t1[3] = {Vec2i(180, 50),  Vec2i(150, 1),   Vec2i(70, 180)}; 
    // Vec2i t2[3] = {Vec2i(180, 150), Vec2i(120, 160), Vec2i(130, 180)}; 
    // triangle(t0[0], t0[1], t0[2], image, red); 
    // triangle(t1[0], t1[1], t1[2], image, white); 
    // triangle(t2[0], t2[1], t2[2], image, green);

    Model model("african_head/african_head.obj");

    for (int f = 0; f < model.nfaces(); f++) {
        auto w0 = model.vert(f, 0);
        auto w1 = model.vert(f, 1);
        auto w2 = model.vert(f, 2);
        
        auto v0 = Vec2i(float(w0.x + 1) * width / 2, float(w0.y + 1) * height / 2);
        auto v1 = Vec2i(float(w1.x + 1) * width / 2, float(w1.y + 1) * height / 2);
        auto v2 = Vec2i(float(w2.x + 1) * width / 2, float(w2.y + 1) * height / 2);

        triangle(v0, v1, v2, image, TGAColor(rand()%255, rand()%255, rand()%255, 255));
    }

    image.flip_vertically();
    image.write_tga_file("output.tga");
    return 0;
}
