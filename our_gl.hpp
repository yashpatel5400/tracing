#include <vector>
#include "tgaimage.h"
#include "model.h"

class IShader {
public:
    virtual Vec3i vertex(int iface, int nthvert) = 0;
    virtual bool fragment(Vec3f bar) = 0;

    TGAColor color;
};

Vec4f toHomogenous(Vec3f vec);
Vec3f fromHomogenous(Vec4f vec);
Matrix perspectiveShift(float cameraPos);
Matrix lookAt(Vec3f eye,  Vec3f center, Vec3f up);
Matrix viewport(int x, int y, int w, int h);
void line(Vec2i v0, Vec2i v1, TGAImage& image, TGAColor color);
void triangleLineSweep(Vec2i v0, Vec2i v1, Vec2i v2, TGAImage& image, TGAColor color);
void triangle(const std::vector<Vec3i>& screenCoords, IShader& shader, std::vector<float>& zbuffer, TGAImage& image);
