#pragma once

#include "V3.h"
#include "M33.h"
#include "ppc.h"
#include "framebuffer.h"
#include "AABB.h"

class TMesh {
public:
	V3* verts;
	int vertsN;
	unsigned int* tris;
	int trisN;
	unsigned char* visTris;


	int onFlag;
	int msiFlag;
	int reflectorFlag;
	V3* colors;
	V3* normals;
	FrameBuffer* texture;
	TMesh() : verts(0), vertsN(0), tris(0), trisN(0), colors(0), normals(0),
		onFlag(1), msiFlag(1), texture(0), reflectorFlag(0), visTris(0) {};
	void SetToCube(V3 cc, float sideLength, unsigned int color0, unsigned int color1);
	void SetQuad(V3* qverts, V3* qcolors);
	void Allocate(int _vertsN, int _trisN);
	void DrawCubeQuadFaces(FrameBuffer* fb, PPC* ppc, unsigned int color);
	void DrawWireFrame(FrameBuffer* fb, PPC* ppc, unsigned int color);
	void LoadBin(char* fname);
	V3 GetCenter();
	void SetCenter(V3 center);
	void Translate(V3 tv);
	void Rotate(V3 aO, V3 aDir, float theta);
	void RenderFilled(FrameBuffer* fb, PPC* ppc, V3 C, V3 L, float ka);
	V3 SetEEQ(V3 v0, V3 v1, V3 v2);
	M33 SetEEQs(V3 pv0, V3 pv1, V3 pv2);
	M33 SetSSIM(V3 pv0, V3 pv1, V3 pv2);
	void Light(V3 C, V3 L, float ka);
	void RayTrace(FrameBuffer* fb, PPC* ppc);
	int IntersectTriangleWithRay(V3 v0, V3 v1, V3 v2, V3 O, V3 dir, float& t, V3& barycc);
	int RayTrace(V3 rO, V3 rdir, V3& rc, float& rz, V3& rrO, V3& rrdir);
	int RayTriangleIntersection(V3 rO, V3 rdir, V3 V0, V3 V1, V3 V2,
		V3& currc, float& currz);
	void SetAABB(AABB& aabb);
	void Scale(float scf);
	void RenderHW();
	void SetRectangle(float rw, float rh);
	void SetTesselatedRectangle(float rw, float rh, int colsN, int rowsN, V3 col);
	void RenderRC(FrameBuffer* fb, PPC* ppc, float vcr);
	void ClearVisibleTriangles();
	void AddVisibleTriangles(FrameBuffer* fb);
	int CountVisibleTriangles();
};