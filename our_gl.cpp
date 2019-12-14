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
const int depth  = 255;

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

Matrix lookAt(Vec3f eye,  Vec3f center, Vec3f up) {
    Vec3f z = (eye - center).normalize();
    Vec3f x = cross(up, z).normalize();
    Vec3f y = cross(z, x).normalize();

    Matrix m = Matrix::identity();
    for (int i = 0; i < 3; i++) {
        m[0][i] = x[i];
        m[1][i] = y[i];
        m[2][i] = z[i];
        m[i][3] = -center[i];
    }
    return m;
}

Matrix viewport(int x, int y, int w, int h) {
    Matrix m = Matrix::identity();
    m[0][3] = x+w/2.f;
    m[1][3] = y+h/2.f;
    m[2][3] = depth/2.f;

    m[0][0] = w/2.f;
    m[1][1] = h/2.f;
    m[2][2] = depth/2.f;
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

void triangle(const std::vector<Vec3i>& screenCoords, IShader& shader, std::vector<float>& zbuffer, TGAImage& image) {
    auto v0 = screenCoords[0];
    auto v1 = screenCoords[1];
    auto v2 = screenCoords[2];

    std::vector<int> xs = {v0.x, v1.x, v2.x};
    std::vector<int> ys = {v0.y, v1.y, v2.y};

    int minX = std::max(std::min(std::min(xs[0], xs[1]), xs[2]), 0);
    int maxX = std::min(std::max(std::max(xs[0], xs[1]), xs[2]), width - 1);
    int minY = std::max(std::min(std::min(ys[0], ys[1]), ys[2]), 0);
    int maxY = std::min(std::max(std::max(ys[0], ys[1]), ys[2]), height - 1);

    auto AB = v1 - v0;
    auto AC = v2 - v0;

    for (int x = minX; x <= maxX; x++) {
        for (int y = minY; y <= maxY; y++) {
            Vec3i p(x, y, 0);
            auto PA = v0 - p;
            Vec3f first(AB.x, AC.x, PA.x);
            Vec3f second(AB.y, AC.y, PA.y);
            
            Vec3f res = cross(first, second);
            bool foundNeg = false;

            if (res.z != 0) {
                auto u = float(res.x) / res.z;
                auto v = float(res.y) / res.z;
                auto t = 1.0 - u - v;
                Vec3f bar(t, u, v);
                
                int curPos = y * width + x;
                float z = screenCoords[0].z * t + screenCoords[1].z * u + screenCoords[2].z * v;
                if (u >= 0 && v >= 0 && t >= 0 && z > zbuffer[curPos]) {
                    shader.fragment(bar);
                    image.set(x, y, shader.color);
                    zbuffer[curPos] = z;
                }
            }
        }   
    }
}