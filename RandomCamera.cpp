#include "stdafx.h"

#include "RandomCamera.h"

RandomCamera::RandomCamera(float hfov, int w, int h, float d0, float d1,
	float _maxD0, float _maxD1, float _dnear, float _dfar) : maxD0(_maxD0),
	maxD1(_maxD1) {

	ppc0 = new PPC(hfov, w, h);
	// scale ppc frustum to put ppc image plane on d0
	float f = ppc0->GetF();
	float scf = d0 / f;
	ppc0->a = ppc0->a * scf;
	ppc0->b = ppc0->b * scf;
	ppc0->c = ppc0->c * scf;

	dnear = d0 / _dnear;
	dfar = d0 / _dfar;

	ppc1 = new PPC(hfov, w, h);
	// scale ppc frustum to put ppc image plane on d0
	scf = d1 / f;
	ppc1->a = ppc1->a * scf;
	ppc1->b = ppc1->b * scf;
	ppc1->c = ppc1->c * scf;

	rayb0 = new V3[w*h];
	rayb1 = new V3[w*h];

}

void RandomCamera::ClearRays() {

	for (int uv = 0; uv < ppc0->w*ppc0->h; uv++)
		rayb0[uv] = rayb1[uv] = V3(0.0f, 0.0f, 0.0f);
}

void RandomCamera::PerturbRay(int u, int v, float pixz, V3 &pP, V3 &pR) {

	// ppc ray origin on near image plane
	V3 pix(.5f + (float)u, .5f + (float)v, 0.0f);

	// random 2D direction in near image plane
	float x = (float)rand() / (float)RAND_MAX;
	float y = (float)rand() / (float)RAND_MAX;
	V3 hdir(x*2.0f-1.0f, y*2.0f-1.0f, 0.0f);
	hdir = hdir.Normalized();

	float distDamp;
	if (pixz > dnear)
		distDamp = 1.0f;
	else if (pixz < dfar)
		distDamp = 0.0f;
	else
		distDamp = (pixz - dfar) / (dnear - dfar);

	V3 ppix0 = pix + hdir*maxD0*distDamp;
	pP = ppc0->GetImagePlanePoint(ppix0[0], ppix0[1]);

	V3 ppix1 = pix + hdir*maxD1*distDamp;
	pR = ppc1->GetImagePlanePoint(ppix1[0], ppix1[1]);

	rayb0[(ppc0->h - 1 - v)*ppc0->w + u] = pP;
	rayb1[(ppc0->h - 1 - v)*ppc0->w + u] = pR;

}

void RandomCamera::PointBasedRender(FrameBuffer *rcfb, int pSize,
	FrameBuffer *fb, PPC *ppc) {

	fb->SetBGR(0xFFFFFFFF);
	fb->ClearZB();

	for (int v = 0; v < ppc0->h; v++) {
		for (int u = 0; u < ppc0->w; u++) {
			int uv = (ppc0->h - 1 - v)*ppc0->w + u;
			float currz = rcfb->GetZ(u, v);
			if (currz == 0.0f)
				continue;
			V3 P = rayb0[uv] + (rayb1[uv] - rayb0[uv]).Normalized() / currz;
			fb->Draw3DPoint(P, ppc, rcfb->Get(u, v), pSize);
		}
	}

}



// new constructor; ppc1 is not used;
RandomCamera::RandomCamera(float hfov, int w, int h, float _vcr, float _df) 
	: vcr(_vcr), df(_df) {

	ppc0 = new PPC(hfov, w, h);
	// scale ppc frustum to put ppc image plane on d0
	float f = ppc0->GetF();
	float scf = df / 2.0f / f;
	ppc0->a = ppc0->a * scf;
	ppc0->b = ppc0->b * scf;
	ppc0->c = ppc0->c * scf;

	rayb0 = new V3[w*h];
	rayb1 = new V3[w*h];

}

void RandomCamera::PerturbRay(int u, int v, V3 &ray0, V3 &ray1) {

	// random 2D direction in near image plane
	float x = (float)rand() / (float)RAND_MAX;
	float y = (float)rand() / (float)RAND_MAX;
	V3 pdir(x*2.0f - 1.0f, y*2.0f - 1.0f, 0.0f);
	float amp = (float)rand() / (float)RAND_MAX;
amp = 1.0f;
	pdir = pdir.Normalized()*amp;
	ray0 = ppc0->C + pdir * vcr;
	ray1 = ppc0->GetImagePlanePoint(.5f + (float)u, .5f + (float)v) + pdir * vcr / 2.0f;

	rayb0[(ppc0->h - 1 - v)*ppc0->w + u] = ray0;
	rayb1[(ppc0->h - 1 - v)*ppc0->w + u] = ray1;

}
