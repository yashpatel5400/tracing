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

const int width = 1000;
const int height = 1000;

Vec4f toHomogenous(Vec3f vec) {
    Vec4f m;
    m[0] = vec.x;
    m[1] = vec.y;
    m[2] = vec.z;
    m[3] = 1.0;
    return m;
}

Vec3f fromHomogenous(Vec4f vec) {
    Vec3f res;
    res.x = vec[0] / vec[3];
    res.y = vec[1] / vec[3];
    res.z = vec[2] / vec[3];
    return res;
}

Matrix perspectiveShift(float cameraPos) {
    Matrix m = Matrix::identity();
    m[3][2] = -1.0 / cameraPos;
    return m;
}

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

void triangle(const std::vector<Vec3f>& worldCoords, std::vector<TGAColor>& colors, std::vector<float>& zbuffer, TGAImage& image, float intensity) {
    auto v0 = Vec2i((worldCoords[0].x+1.)*width/2., (worldCoords[0].y+1.)*height/2.);
    auto v1 = Vec2i((worldCoords[1].x+1.)*width/2., (worldCoords[1].y+1.)*height/2.);
    auto v2 = Vec2i((worldCoords[2].x+1.)*width/2., (worldCoords[2].y+1.)*height/2.);

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
                auto u = float(res.x) / res.z;
                auto v = float(res.y) / res.z;
                auto t = 1.0 - u - v;
                
                int curPos = y * width + x;
                float z = worldCoords[0].z * t + worldCoords[1].z * u + worldCoords[2].z * v;
                auto c = TGAColor(
                    colors[0][2] * t + colors[1][2] * u + colors[2][2] * v,
                    colors[0][1] * t + colors[1][1] * u + colors[2][1] * v,
                    colors[0][0] * t + colors[1][0] * u + colors[2][0] * v,
                    255
                ) * intensity;
                if (u >= 0 && v >= 0 && t >= 0 && z > zbuffer[curPos]) {
                    image.set(x, y, c);
                    zbuffer[curPos] = z;
                }
            }
        }   
    }
}

int main() {
    TGAImage image(width, height, TGAImage::RGB);
    Model model("african_head/african_head.obj");

    Vec3f light(0, 0, -1);
    std::vector<float> zbuffer(width * height, -std::numeric_limits<float>::max());
    
    Vec3f cameraPos(0, 0, 3.0);
    auto shiftMat = perspectiveShift(cameraPos.z);
    for (int f = 0; f < model.nfaces(); f++) {
        std::vector<Vec3f> worldCoords = { model.vert(f, 0), model.vert(f, 1), model.vert(f, 2) };
        std::vector<TGAColor> colors = { model.diffuse(model.uv(f, 0)), model.diffuse(model.uv(f, 1)), model.diffuse(model.uv(f, 2)) };

        for (int i = 0; i < 3; i++) {
            worldCoords[i] = fromHomogenous(shiftMat * toHomogenous(worldCoords[i]));
            worldCoords[i].x = std::min(std::max(worldCoords[i].x, -1.f), 1.f);
            worldCoords[i].y = std::min(std::max(worldCoords[i].y, -1.f), 1.f);
            worldCoords[i].z = std::min(std::max(worldCoords[i].z, -1.f), 1.f);
        }
        
        auto norm = cross((worldCoords[2] - worldCoords[0]), (worldCoords[1] - worldCoords[0]));
        norm.normalize();
        auto intensity = light * norm;
        
        if (intensity > 0) {
            triangle(worldCoords, colors, zbuffer, image, intensity);
        }
    }

    image.flip_vertically();
    image.write_tga_file("output.tga");
    return 0;
}
