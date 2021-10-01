#pragma once
#include "v3.h"
#include "framebuffer.h"
#include "cubemap.h"

// declaration of equirectangular camera and image class
// all columns start/end with north/south pole
// middle row is equator
// columns are meridians

class ERCI {
public:
	int w; // number of pixels on row
	int h; // number of pixels on column is w/2
	FrameBuffer *fb; // equirectangular image, i.e. 2D array of pixels
	ERCI(int _w); // constructor, takes resolution, height is half the width
	V3 GetRay(float uf, float vf); // aka unprojection, returns ray through 
			// immage point (uf, vf)
	void SetFromCubeMap(CubeMap *cubeMap); // sets equirectangular image
			// by looking up directions in given cubeMap
	void Project(V3 dir, V3 &pdir); // returns pixel coordinates pdir[0], pdir[1]
			// of intersection between image plane and ray pdir
			// pdir has to be normalized
	unsigned int LookUp(V3 dir); // return color at image location of projection
			// of dir; nearest neighbor lookup
	void SetPPCI(PPC2 *ppc, FrameBuffer *ppci);  // set conventional image ppci
			// with camera ppc from "this" ERCI
		
};