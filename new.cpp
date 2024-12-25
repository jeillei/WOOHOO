#define OLC_PGE_APPLICATION
#include "olcPixelGameEngine.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <algorithm>
//g++ -o my_program new.cpp -lX11 -lXrandr -lGL -lpng


using namespace std;

struct vec3d
{
    float x, y, z;

};

struct triangle
{
    vec3d p[3];
    olc::Pixel color;
};

struct mesh
{
    vector<triangle> tris;

    bool loadfromfile(string sFilename)
    {
        ifstream f(sFilename);
        if (!f.is_open())
        {
            return false;
        }
        vector<vec3d> verts;
        while(!f.eof())
        {
            char line [128];
            f.getline(line, 128);

            stringstream s;
            s<< line;

            char junk;

            if (line[0] == 'v')
            {
                vec3d v;
                s >> junk >> v.x >> v.y >> v.z;
                verts.push_back(v);
            }

            if (line[0] == 'f')
            {
                int f[3];
                s >> junk >> f[0] >> f[1] >> f[2];
                tris.push_back({ verts[f[0]-1], verts[f[1]-1], verts[f[2]-1] });
            }
        }
        return true;
    }
};

struct Matrix
{
    float m[4][4] = { 0 };
};

class Engine3D : public olc::PixelGameEngine
{
public:
    Engine3D()
    {
        sAppName = "3d Demo";
    }


private:
    mesh meshCube;
    mesh meshpyramid;
    mesh meshhexagonal;

    mesh christmas;
    mesh meshshape;

    Matrix matproj;

    vec3d vCamera = {0,0,0};


    float fTheta = 0;
    void multiplymatvec(vec3d &i, vec3d &o, Matrix &m)
    {
        o.x = i.x * m.m[0][0] + i.y * m.m[1][0] + i.z * m.m[2][0] + m.m[3][0];
        o.y = i.x * m.m[0][1] + i.y * m.m[1][1] + i.z * m.m[2][1] + m.m[3][1];
        o.z = i.x * m.m[0][2] + i.y * m.m[1][2] + i.z * m.m[2][2] + m.m[3][2];
        float w = i.x * m.m[0][3] + i.y * m.m[1][3] + i.z * m.m[2][3] + m.m[3][3];
        
        if (w != 0.0f){ 
            o.x /= w; o.y /= w; o.z /= w;
        }
    }

public:
    bool OnUserCreate() override
    {

        // create a unit cube here
        // imagine a 3d drawing of a cube on paper with origin at the bottom left corner
        meshCube.tris = {

        // SOUTH
		{ 0.0f, 0.0f, 0.0f,    0.0f, 1.0f, 0.0f,    1.0f, 1.0f, 0.0f },
		{ 0.0f, 0.0f, 0.0f,    1.0f, 1.0f, 0.0f,    1.0f, 0.0f, 0.0f },

		// EAST                                                      
		{ 1.0f, 0.0f, 0.0f,    1.0f, 1.0f, 0.0f,    1.0f, 1.0f, 1.0f },
		{ 1.0f, 0.0f, 0.0f,    1.0f, 1.0f, 1.0f,    1.0f, 0.0f, 1.0f },

		// NORTH                                                     
		{ 1.0f, 0.0f, 1.0f,    1.0f, 1.0f, 1.0f,    0.0f, 1.0f, 1.0f },
		{ 1.0f, 0.0f, 1.0f,    0.0f, 1.0f, 1.0f,    0.0f, 0.0f, 1.0f },

		// WEST                                                      
		{ 0.0f, 0.0f, 1.0f,    0.0f, 1.0f, 1.0f,    0.0f, 1.0f, 0.0f },
		{ 0.0f, 0.0f, 1.0f,    0.0f, 1.0f, 0.0f,    0.0f, 0.0f, 0.0f },

		// TOP                                                       
		{ 0.0f, 1.0f, 0.0f,    0.0f, 1.0f, 1.0f,    1.0f, 1.0f, 1.0f },
		{ 0.0f, 1.0f, 0.0f,    1.0f, 1.0f, 1.0f,    1.0f, 1.0f, 0.0f },

		// BOTTOM                                                    
		{ 1.0f, 0.0f, 1.0f,    0.0f, 0.0f, 1.0f,    0.0f, 0.0f, 0.0f },
		{ 1.0f, 0.0f, 1.0f,    0.0f, 0.0f, 0.0f,    1.0f, 0.0f, 0.0f },

        };

        meshpyramid.tris = {
        // Base vertices (4 points for the square base)
        { -1.0f, 0.0f, -1.0f,    1.0f, 0.0f, -1.0f,    1.0f, 0.0f, 1.0f },
        { -1.0f, 0.0f, -1.0f,    1.0f, 0.0f, 1.0f,    -1.0f, 0.0f, 1.0f },

        // Apex of the pyramid (one point at the top)
        { 0.0f, 2.0f, 0.0f },

        // Side faces (4 triangles connecting apex to each edge of the base)
        { 0.0f, 2.0f, 0.0f,    -1.0f, 0.0f, -1.0f,    1.0f, 0.0f, -1.0f }, // front
        { 0.0f, 2.0f, 0.0f,    1.0f, 0.0f, -1.0f,    1.0f, 0.0f, 1.0f },  // right
        { 0.0f, 2.0f, 0.0f,    1.0f, 0.0f, 1.0f,    -1.0f, 0.0f, 1.0f }, // back
        { 0.0f, 2.0f, 0.0f,    -1.0f, 0.0f, 1.0f,    -1.0f, 0.0f, -1.0f },// left
        };


        meshhexagonal.tris = {
            // First hexagonal base (6 triangles forming the first hexagon)
            { 0.0f, 0.0f, 0.0f,   1.0f, 0.0f, 0.0f,   0.5f, 0.866f, 0.0f },  // First triangle
            { 0.0f, 0.0f, 0.0f,   0.5f, 0.866f, 0.0f,   -0.5f, 0.866f, 0.0f }, // Second triangle
            { 0.0f, 0.0f, 0.0f,   -0.5f, 0.866f, 0.0f,   -1.0f, 0.0f, 0.0f }, // Third triangle
            { 0.0f, 0.0f, 0.0f,   -1.0f, 0.0f, 0.0f,    -0.5f, -0.866f, 0.0f }, // Fourth triangle
            { 0.0f, 0.0f, 0.0f,   -0.5f, -0.866f, 0.0f,  0.5f, -0.866f, 0.0f }, // Fifth triangle
            { 0.0f, 0.0f, 0.0f,   0.5f, -0.866f, 0.0f,   1.0f, 0.0f, 0.0f },  // Sixth triangle

            // Second hexagonal base (6 triangles forming the second hexagon at a different Z-level)
            { 0.0f, 0.0f, 2.0f,   1.0f, 0.0f, 2.0f,   0.5f, 0.866f, 2.0f },  // First triangle (top)
            { 0.0f, 0.0f, 2.0f,   0.5f, 0.866f, 2.0f,   -0.5f, 0.866f, 2.0f }, // Second triangle (top)
            { 0.0f, 0.0f, 2.0f,   -0.5f, 0.866f, 2.0f,   -1.0f, 0.0f, 2.0f }, // Third triangle (top)
            { 0.0f, 0.0f, 2.0f,   -1.0f, 0.0f, 2.0f,    -0.5f, -0.866f, 2.0f }, // Fourth triangle (top)
            { 0.0f, 0.0f, 2.0f,   -0.5f, -0.866f, 2.0f,  0.5f, -0.866f, 2.0f }, // Fifth triangle (top)
            { 0.0f, 0.0f, 2.0f,   0.5f, -0.866f, 2.0f,   1.0f, 0.0f, 2.0f },  // Sixth triangle (top)

            // Side faces (6 rectangles connecting top and bottom hexagons)
            { 0.0f, 0.0f, 0.0f,   1.0f, 0.0f, 0.0f,   1.0f, 0.0f, 2.0f },  // Front-right side
            { 0.0f, 0.0f, 0.0f,   1.0f, 0.0f, 0.0f,   0.5f, 0.866f, 2.0f },  // Right side
            { 0.0f, 0.0f, 0.0f,   0.5f, 0.866f, 0.0f,   0.5f, 0.866f, 2.0f },  // Top-right side
            { 0.0f, 0.0f, 0.0f,   0.5f, 0.866f, 0.0f,   -0.5f, 0.866f, 2.0f }, // Back-right side
            { 0.0f, 0.0f, 0.0f,   -0.5f, 0.866f, 0.0f,   -0.5f, 0.866f, 2.0f }, // Top-back side
            { 0.0f, 0.0f, 0.0f,   -0.5f, 0.866f, 0.0f,   -1.0f, 0.0f, 2.0f },  // Back-left side
            { 0.0f, 0.0f, 0.0f,   -1.0f, 0.0f, 0.0f,    -1.0f, 0.0f, 2.0f },  // Top-left side
            { 0.0f, 0.0f, 0.0f,   -1.0f, 0.0f, 0.0f,    -0.5f, -0.866f, 2.0f },// Left side
            { 0.0f, 0.0f, 0.0f,   -0.5f, -0.866f, 0.0f,  -0.5f, -0.866f, 2.0f }, // Bottom-left side
            { 0.0f, 0.0f, 0.0f,   -0.5f, -0.866f, 0.0f,   0.5f, -0.866f, 2.0f }, // Bottom-right side
            { 0.0f, 0.0f, 0.0f,   0.5f, -0.866f, 0.0f,   0.5f, -0.866f, 2.0f }, // Bottom-right side
            { 0.0f, 0.0f, 0.0f,   0.5f, -0.866f, 0.0f,   1.0f, 0.0f, 2.0f }   // Right side
        };

        meshshape.loadfromfile("VideoShip.obj");
        christmas.loadfromfile("tree.obj");


        //projection matrix
        float fnear = 0.1f;
        float ffar = 1000.0f;
        float ffov = 90.0f;
        float fAsRa = (float)ScreenHeight() / (float)ScreenWidth();
        float fFOVRAd = 1.0f/tanf(ffov * 0.5f / 180.0f * 3.14159f);

        matproj.m[0][0] = fAsRa * fFOVRAd;
        matproj.m[1][1] = fAsRa;
        matproj.m[2][2] = ffar/(ffar-fnear);
        matproj.m[3][2] = (-ffar * fnear) / (ffar - fnear);
        matproj.m[2][3] = 1.0f;
        matproj.m[3][3] = 0.0f;

        return true;
    }
    bool OnUserUpdate(float fElapsedTime) override
    {
        // Clear Screen
        Clear(olc::Pixel(0, 0, 0));  // Clear the screen with black color

        Matrix matRotZ, matRotX, matRotY;
        fTheta += 0.0001f + fElapsedTime;

        // Rotation Z 
		matRotZ.m[0][0] = cosf(fTheta);
		matRotZ.m[0][1] = sinf(fTheta);
		matRotZ.m[1][0] = -sinf(fTheta);
		matRotZ.m[1][1] = cosf(fTheta);
		matRotZ.m[2][2] = 1;
		matRotZ.m[3][3] = 1;

		// Rotation X
		matRotX.m[0][0] = 1;
		matRotX.m[1][1] = cosf(fTheta * 0.5f);
		matRotX.m[1][2] = sinf(fTheta * 0.5f);
		matRotX.m[2][1] = -sinf(fTheta * 0.5f);
		matRotX.m[2][2] = cosf(fTheta * 0.5f);
		matRotX.m[3][3] = 1;

        // Rotation Y
        matRotY.m[0][0] = cosf(fTheta);
        matRotY.m[0][2] = sinf(fTheta);
        matRotY.m[1][1] = 1;
        matRotY.m[2][0] = -sinf(fTheta);
        matRotY.m[2][2] = cosf(fTheta);
        matRotY.m[3][3] = 1;


        //store triangle to rasterize later
        vector<triangle> toraster;

        for(auto tri: christmas.tris)
        {
            triangle triProjected, triTranslated, triRotatedZ, triRotatedZX;

            // Rotate in Z-Axis
			multiplymatvec(tri.p[0], triRotatedZ.p[0], matRotY);
			multiplymatvec(tri.p[1], triRotatedZ.p[1], matRotY);
			multiplymatvec(tri.p[2], triRotatedZ.p[2], matRotY);

			// Rotate in X-Axis
			multiplymatvec(triRotatedZ.p[0], triRotatedZX.p[0], matRotZ);
			multiplymatvec(triRotatedZ.p[1], triRotatedZX.p[1], matRotZ);
			multiplymatvec(triRotatedZ.p[2], triRotatedZX.p[2], matRotZ);

            //flip upside down

            triRotatedZ.p[0].y *= -1;
            triRotatedZ.p[1].y *= -1;
            triRotatedZ.p[2].y *= -1;

            //offset into screen
            triTranslated = triRotatedZ;
            triTranslated.p[0].z = triRotatedZ.p[0].z + 8.0f;
            triTranslated.p[1].z = triRotatedZ.p[1].z + 8.0f;
            triTranslated.p[2].z = triRotatedZ.p[2].z + 8.0f;

            // triTranslated.p[0].y = triRotatedZX.p[0].y + 8.0f;
            // triTranslated.p[1].y = triRotatedZX.p[1].y + 8.0f;
            // triTranslated.p[2].y = triRotatedZX.p[2].y + 8.0f;

            //normal lines
            vec3d normal, line1, line2;
            line1.x = triTranslated.p[1].x - triTranslated.p[0].x;
            line1.y = triTranslated.p[1].y - triTranslated.p[0].y;
            line1.z = triTranslated.p[1].z - triTranslated.p[0].z;

            line2.x = triTranslated.p[2].x - triTranslated.p[0].x;
            line2.y = triTranslated.p[2].y - triTranslated.p[0].y;
            line2.z = triTranslated.p[2].z - triTranslated.p[0].z;

            normal.x = line1.y * line2.z - line1.z * line2.y;
            normal.y = line1.z * line2.x - line1.x * line2.z;
            normal.z = line1.x * line2.y - line1.y * line2.x;

            float l = sqrtf(normal.x*normal.x + normal.y*normal.y + normal.z*normal.z);
            normal.x/=l;normal.y/=l;normal.z/=l;


            

            //only trianlges that are away blocked from the eye
            // if(normal.z < 0)
            if ((normal.x * (triTranslated.p[0].x - vCamera.x) + 
            normal.y * (triTranslated.p[0].y - vCamera.y) +
            normal.z * (triTranslated.p[0].z - vCamera.z)) < 0.0f)
            {

                vec3d light_direction = { 0.0f, 0.0f, -1.0f};
                float l = sqrtf(light_direction.x*light_direction.x + light_direction.y*light_direction.y + light_direction.z*light_direction.z);
                light_direction.x/=l; light_direction.y/=l;light_direction.z/=l;

                float dotProduct = max(0.1f, normal.x * light_direction.x + normal.y * light_direction.y + normal.z * light_direction.z);
                olc::Pixel triColor = olc::Pixel(144 * dotProduct, 255 * dotProduct, 144 * dotProduct);

                triTranslated.color = triColor;

                //project 3D -> 2D
                multiplymatvec(triTranslated.p[0],triProjected.p[0],matproj);
                multiplymatvec(triTranslated.p[1],triProjected.p[1],matproj);
                multiplymatvec(triTranslated.p[2],triProjected.p[2],matproj);

                triProjected.color = triTranslated.color;

                triProjected.p[0].x +=1.0f;triProjected.p[0].y +=1.0f;
                triProjected.p[1].x +=1.0f;triProjected.p[1].y +=1.0f;
                triProjected.p[2].x +=1.0f;triProjected.p[2].y +=1.0f;

                triProjected.p[0].x *= 0.5f *(float)ScreenWidth();triProjected.p[0].y *= 0.5f *(float)ScreenHeight();
                triProjected.p[1].x *= 0.5f *(float)ScreenWidth();triProjected.p[1].y *= 0.5f *(float)ScreenHeight();
                triProjected.p[2].x *= 0.5f *(float)ScreenWidth();triProjected.p[2].y *= 0.5f *(float)ScreenHeight();

                toraster.push_back(triProjected);
            //     FillTriangle(
            // (int)triProjected.p[0].x, (int)triProjected.p[0].y,
            // (int)triProjected.p[1].x, (int)triProjected.p[1].y,
            // (int)triProjected.p[2].x, (int)triProjected.p[2].y,
            // triColor);
            }

        }

        sort(toraster.begin(), toraster.end(), [](triangle &t1, triangle &t2)
		{
			float z1 = (t1.p[0].z + t1.p[1].z + t1.p[2].z) / 3.0f;
			float z2 = (t2.p[0].z + t2.p[1].z + t2.p[2].z) / 3.0f;
			return z1 < z2;
		});


        for (auto &tri: toraster)
        {
            FillTriangle(
            (int)tri.p[0].x, (int)tri.p[0].y,
            (int)tri.p[1].x, (int)tri.p[1].y,
            (int)tri.p[2].x, (int)tri.p[2].y,
            tri.color);

            // DrawTriangle((int)tri.p[0].x, (int)tri.p[0].y,
            // (int)tri.p[1].x, (int)tri.p[1].y,
            // (int)tri.p[2].x, (int)tri.p[2].y,
            // olc::Pixel (0,0,0));
        }
        return true;
    } 
};   

int main()
{
    Engine3D demo;
    if(demo.Construct(256, 240, 4, 4)) 
    {
        demo.Start();
    }
    return 0;
}
