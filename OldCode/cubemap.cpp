#include "stdafx.h"
#include "cubemap.h"

CubeMap::CubeMap(int w) {

	ppcs = new PPC2[6];
	PPC2 tmp(90.0f, w, w);
	ppcs[0] = tmp;
	tmp.Pan(90.0f);
	ppcs[1] = tmp;
	tmp.Pan(90.0f);
	ppcs[2] = tmp;
	tmp.Pan(90.0f);
	ppcs[3] = tmp;
	tmp.Pan(90.0f);
	tmp.Tilt(90.0f);
	ppcs[4] = tmp;
	tmp.Tilt(90.0f);
	tmp.Tilt(90.0f);
	ppcs[5] = tmp;

	faces = new FrameBuffer*[6];
	int i = 0;
	faces[i++] = new FrameBuffer(30, 30, w, w, 0);
	faces[i++] = new FrameBuffer(30, 30, w, w, 0);
	faces[i++] = new FrameBuffer(30, 30, w, w, 0);
	faces[i++] = new FrameBuffer(30, 30, w, w, 0);
	faces[i++] = new FrameBuffer(30, 30, w, w, 0);
	faces[i++] = new FrameBuffer(30, 30, w, w, 0);

	unsigned int cols[2] = {0xFFFFFFFF, 0xF000000};
	i = 0;
	int checkerSize = 16;
	cols[1] = 0xFF0000FF;
	faces[i++]->SetChecker(cols, checkerSize);
	cols[1] = 0xFF00FF00;
	faces[i++]->SetChecker(cols, checkerSize);
	cols[1] = 0xFFFF0000;
	faces[i++]->SetChecker(cols, checkerSize);
	cols[1] = 0xFF00FFFF;
	faces[i++]->SetChecker(cols, checkerSize);
	cols[1] = 0xFFFFFF00;
	faces[i++]->SetChecker(cols, checkerSize);
	cols[1] = 0xFFFF00FF;
	faces[i++]->SetChecker(cols, checkerSize);


}

unsigned int CubeMap::LookUp(V3 dir) {

	for (int fi = 0; fi < 6; fi++) {
		V3 pdir;
		if (!ppcs[fi].Project(dir, pdir))
			continue;
		if (!ppcs[fi].InsideImage(pdir))
			continue;
		return faces[fi]->Get((int)pdir[0], (int)pdir[1]);
	}

	return 0xFF888888;

}

void CubeMap::SetPPCI(PPC2 *ppc, FrameBuffer *fb) {

	for (int v = 0; v < ppc->h; v++) {
		for (int u = 0; u < ppc->w; u++) {
			V3 dir = ppc->GetRay(.5f + (float)u, .5f + (float)v);
			fb->Set(u, v, LookUp(dir));
		}
	}

}
