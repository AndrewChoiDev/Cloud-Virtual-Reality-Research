#pragma once

#include "gui.h"
#include "framebuffer.h"
#include "ppc.h"
#include "TMesh.h"
#include "CGInterface.h"
#include "RandomCamera.h"
#include "Path.h"

class Scene {
public:

	CGInterface * cgi;
	ShaderOneInterface *soi;

	GUI *gui;
	FrameBuffer *fb, *fb3, *hwfb, *gpufb;
	PPC *ppc, *ppc3;
	TMesh *tmeshes;
	int tmeshesN;
	int renderHW;
	Scene();
	void DBG();
	void NewButton();
	void Render();
	void Render(FrameBuffer *rfb, PPC *rppc); // projection followed by rasterization
	void RenderHW();
	void RenderGPU();
	void RenderRayTracing(FrameBuffer *rfb, PPC *rppc); // ray tracing
	int RayTrace(V3 rO, V3 rdir, int rayOrder, V3& rc, float &currz);
	float vf; // ppc visualization focal length
	V3 L; // point light source
	float ka; // ambient lighting coefficient
	float p1; // test parameter for fragment shader

	void Render(RandomCamera *rc, FrameBuffer *rcfb);
	void RCSetup();
	float morphFraction; // morphing parameters
	void PerSessionHWInit();
	int hwPerSessionInit;

	GLuint hwtex;
	void SetupSM();
	int smSetup;
	void SetupTeapot();

	void ResampleERI(FrameBuffer *eri, FrameBuffer *ppi, PPC *ppc);
	void ResamplePPC(FrameBuffer *ppi0, PPC *ppc0, FrameBuffer *ppi1, PPC *ppc1);
	int IntersectRayWithSphere(V3 sO, float r, V3 rO, V3 rDir, V3 &P);
	int IntersetRayWithPlane(V3 rayO, V3 rayDir, V3 pPoint, V3 pNormal, V3 &intP);
	unsigned int LookupRayIntoERI(V3 ray, FrameBuffer *eri);
	V3 ReflectRay(V3 ray, V3 n);
	FrameBuffer *eri;
	V3 dsc;
	float dsr;

	void MakeDotPattern();

	// AR image tracking demo
	FrameBuffer *imageToTrack, *currFrame, *visfb;
	V3 ittCorners[4]; // bottom left, bottom right, top right, top left
	PPC *currPPC;
	// find the retangular image itt into frame and set view of current frame
	void TrackImage(FrameBuffer *itt, V3 *corners, PPC *ppc0, FrameBuffer *frame, PPC *ppc);
	void TrackImageSetup();
	void DrawImageRectangle(FrameBuffer *frame, PPC *ppc, V3 *corners, V3 color);
	float ColorDifference(FrameBuffer *itt, V3 *corners, FrameBuffer *frame, PPC *ppc);
	void ImageTracking();
	void TransparentDisplay();
	void TransparentDisplayPlanarProxy();
	float piPhi, piDelta;
	void VisTransp(FrameBuffer *visfb, FrameBuffer *tex, float d, float h,
		PPC *userCamera, PPC* tabletCamera, PPC *ppc, int isConv);
	int IntWall(V3 rayO, V3 ray, float d, float ph, FrameBuffer *tex,
		unsigned int &col);
	void VisDisp(FrameBuffer *disp, FrameBuffer *tex, float d, float ph, 
		PPC *userCam, PPC *tabletCam, int cond);

	// Edge VR Project
	void EdgeVRSetup();
	Path *path;
	void PlaybackPath(float fps);
	void CollectVisibleTrianglesOnPath(float t0, float t1, float fps);
};

extern Scene *scene;