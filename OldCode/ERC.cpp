#include "stdafx.h"
#include "ERC.h"

ERCI::ERCI(int _w) {

	w = _w;
	h = w / 2;
	fb = new FrameBuffer(10, 10, w, h, 0);

}

V3 ERCI::GetRay(float uf, float vf) {

	float theta = uf / (float)w * 360.0f;
	float phi = 90.0f - vf / (float)h *180.0f;
	V3 R = V3(0.0f, 0.0f, 1.0f);
	R = R.RotatePointAboutAxis(V3(0.0f, 0.0f, 0.0f), V3(0.0f, -1.0f, 0.0f), theta);
	V3 a = V3(0.0f, 1.0f, 0.0f) ^ R;
	R = R.RotatePointAboutAxis(V3(0.0f, 0.0f, 0.0f), a, -phi);
	return R;
}

void ERCI::SetFromCubeMap(CubeMap *cubeMap) {

	for (int v = 0; v < h; v++) {
		for (int u = 0; u < w; u++) {
			V3 currRay = GetRay(.5f + (float)u, .5f + (float)v);
			fb->Set(u, v, cubeMap->LookUp(currRay));
		}
	}

}

void ERCI::Project(V3 dir, V3 &pdir) {

	float phi = asinf(dir[1])*180.0f/3.14159627f;
	float theta = atanf(dir[2] / dir[0])*180.0f / 3.14159627f;
	theta = (dir[0] > 0.0f) ? 270.0f + theta : 90.0f + theta;

	float uf = theta / 360.0f*(float)w;
	float vf = (90.0f - phi) / 180.0f * (float)h;

	pdir[0] = uf;
	pdir[1] = vf;
	pdir[2] = 1.0f;

//	float theta = uf / (float)w * 360.0f;
//	float phi = 90.0f - vf / (float)h *180.0f;


}

unsigned int ERCI::LookUp(V3 dir) {

	V3 pdir;
	Project(dir, pdir);
	return fb->Get((int)pdir[0], (int)pdir[1]);

}

void ERCI::SetPPCI(PPC2 *ppc, FrameBuffer *ppci) {

	for (int v = 0; v < ppc->h; v++) {
		for (int u = 0; u < ppc->w; u++) {
			V3 dir = ppc->GetRay(.5f + (float)u, .5f + (float)v);
			dir = dir.UnitVector();
			ppci->Set(u, v, LookUp(dir));
		}
	}

}

