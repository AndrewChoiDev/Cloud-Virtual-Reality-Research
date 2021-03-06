#pragma once

#include <FL/Fl.H>
#include <FL/Fl_Gl_Window.H>
#include <GL/glut.h>

#include "ppc.h"
class PPC;

class FrameBuffer : public Fl_Gl_Window {
public:
	unsigned int *pix; // pixel array
	float *zb;
	int w, h;
	V3 *xyz;
	unsigned int *trID;
	int isHW;
	float tstep, rstep;
	FrameBuffer(int u0, int v0, int _w, int _h, unsigned int _id);
	void draw();
	void KeyboardHandle();
	int handle(int guievent);
	void SetBGR(unsigned int bgr);
	void Set(int u, int v, unsigned int color);
	void SetXYZ(int u, int v, V3 p);
	void SetZB(int u, int v, float z);
	unsigned int Get(int u, int v);
	float GetZ(int u, int v);
	void LoadTiff(char* fname);
	void SaveAsTiff(char *fname);
	void SetChecker(unsigned int col0, unsigned int col1, int csize);
	void Draw3DSegment(V3 P0, V3 P1, PPC *ppc, V3 c0, V3 c1);
	void Draw2DSegment(V3 p0, V3 c0, V3 p1, V3 c1);
	void DrawSquarePoint(float uf, float vf, int psize, unsigned int color);
	void DrawDiagonalCross(float uf, float vf, int csize, unsigned int color);
	void DrawDiagonalCrossPadded(float uf, float vf, int csize, unsigned int color);
	void DrawStraightCross(float uf, float vf, int csize, unsigned int color);
	void DrawStraightCrossPadded(float uf, float vf, int csize, unsigned int color);
	void DrawDisk(float uf, float vf, float r, unsigned int color);
	int Farther(int u, int v, float currz);
	void ClearZB();
	void Draw3DPoint(V3 P, PPC *ppc, unsigned int color, int psize);
	// given an image with its camera ppc0, create a new image fb1, with camera ppc1
	// ppc0 and ppc1 have the same eye, that is V3 C
	void GeneralizedRotation(PPC *ppc0, FrameBuffer *fb1, PPC *ppc1);
	unsigned int Lookup(float uf, float vf);
	void SetFromImage(FrameBuffer *fb);
	void ClearItemBuffer();
	void SetItemBuffer(int u, int v, unsigned int tri);
};