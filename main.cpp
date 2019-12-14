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

void triangle(const std::vector<Vec3i>& worldCoords, std::vector<TGAColor>& colors, std::vector<float>& zbuffer, TGAImage& image, const std::vector<float>& intensity) {
    auto v0 = worldCoords[0];
    auto v1 = worldCoords[1];
    auto v2 = worldCoords[2];

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
                
                int curPos = y * width + x;
                float z = worldCoords[0].z * t + worldCoords[1].z * u + worldCoords[2].z * v;
                float aggIntensity = t * intensity[0] + u * intensity[1] + v * intensity[2];

                auto c = TGAColor(
                    colors[0][2] * t + colors[1][2] * u + colors[2][2] * v,
                    colors[0][1] * t + colors[1][1] * u + colors[2][1] * v,
                    colors[0][0] * t + colors[1][0] * u + colors[2][0] * v,
                    255
                ) * aggIntensity;
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

    auto light = Vec3f(1,-1,1).normalize();
    Vec3f cameraPos(1.0, 2.0, 3.0);
    Vec3f center(0, 0.0, 0);

    std::vector<float> zbuffer(width * height, -std::numeric_limits<float>::max());

    auto perspCoords = perspectiveShift(cameraPos.z);
    auto camCoords = lookAt(cameraPos, center, Vec3f(0,1,0));;
    auto viewportCoords = viewport(width / 8, height / 8, width * 3/4, height * 3/4);

    for (int f = 0; f < model.nfaces(); f++) {
        std::vector<Vec3f> worldCoords;
        std::vector<Vec3i> screenCoords;
        std::vector<TGAColor> colors;
        std::vector<float> intensity;

        for (int i = 0; i < 3; i++) {
            auto worldCoord = model.vert(f, i);
            auto normal = model.normal(f, i);
            auto uv = model.uv(f, i);
            
            worldCoords.push_back(worldCoord);
            screenCoords.push_back(fromHomogenous(viewportCoords * perspCoords * camCoords * toHomogenous(worldCoord)));
            colors.push_back(model.diffuse(uv));
            intensity.push_back(light * normal);
        }
        triangle(screenCoords, colors, zbuffer, image, intensity);
    }

    image.flip_vertically();
    image.write_tga_file("output.tga");
    return 0;
}
