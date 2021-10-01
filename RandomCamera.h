#pragma once

#include "V3.h"
#include "ppc.h"

class RandomCamera {
public:
// OLD
	PPC *ppc0, *ppc1; // cameras with identical rays that use d0 and d1 as image planes
	float maxD0, maxD1; // distortion amounts in pixels on d0 and d1 planes
	float dnear, dfar; // distances that define the range of depths where distortion
							// tapers off to 0
	RandomCamera(float hfov, int w, int h, float d0, float d1,
		float _maxD0, float _maxD1, float _dnear, float _dfar);
	void PerturbRay(int u, int v, float z, V3 &r0, V3 &r1);
// END OLD

// NEW
	// PPC *ppc0; // uses ppc0 from above; not ppc1
	float vcr; // view circle radius
	float df;// depth where ray perturbation tapers off to 0
	RandomCamera(float hfov, int w, int h, float _vcr, float _df);
	void PerturbRay(int u, int v, V3 &r0, V3 &r1);
// END NEW

	V3 *rayb0, *rayb1; // stores ray tails and tips
	void ClearRays();
	void PointBasedRender(FrameBuffer *rcfb, int pSize, 
		FrameBuffer *fb, PPC *ppc);
};