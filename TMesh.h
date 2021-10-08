#pragma once

#include "V3.h"
#include "M33.h"
#include "ppc.h"
#include "framebuffer.h"
#include "AABB.h"

class TMesh {
public:
	const char VISIBLE_TRI = 1;
	const char INVISIBLE_TRI = 0;

	V3* verts;
	int vertsN;
	unsigned int* tris; // size = trisN * 3
	int trisN;
	unsigned char* visTris; // size = trisN


	int onFlag;
	int msiFlag;
	int reflectorFlag;
	int renderOnlyVisTrisFlag;
	V3* colors;
	V3* normals;
	FrameBuffer* texture;
	TMesh() : verts(0), vertsN(0), tris(0), trisN(0), colors(0), normals(0),
		onFlag(1), msiFlag(1), texture(0), reflectorFlag(0), renderOnlyVisTrisFlag(0), visTris(0) {
		ClearVisibleTriangles();
	};
	TMesh& operator= (const TMesh& o) {
		this->verts = o.verts;
		this->tris = o.tris;
		this->visTris = o.visTris;
		this->colors = o.colors;
		this->normals = o.normals;
		this->texture = o.texture;

		this->vertsN = o.vertsN;
		this->trisN = o.trisN;

		this->onFlag = o.onFlag;
		this->msiFlag = o.msiFlag;
		this->reflectorFlag = o.reflectorFlag;
		this->renderOnlyVisTrisFlag = o.renderOnlyVisTrisFlag;

		return *this;
	}
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
	void SetAllTrianglesVisible();
	void AddVisibleTriangles(FrameBuffer* fb);
	int CountVisibleTriangles();
	void explodeMesh();
	bool isMeshExploded();
	void colorWithIndices();
	void setVisibleTrianglesHWFrameBuffer(FrameBuffer* hwfb);
	void addVisibleTrianglesHWFrameBuffer(FrameBuffer* hwfb);
	TMesh constructVisibleMesh();
};