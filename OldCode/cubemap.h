#pragma once

#include "Cam.h"

#include "framebuffer.h"

// declaration of CubeMap class
// six images and six cameras
// the images have square image planes
// the cameras have view directions aimed at center of faces of a cube
// the cameras have 90 degree horizontal field of view, and therefore also
//		90 degree vertical field of view, since the resolution is "square"

class CubeMap {
public:
	PPC2 *ppcs;
	FrameBuffer **faces;
	CubeMap(int w); // build cubemap; set colors matching a world that is a cube
		// with checkered faces; 
	unsigned int LookUp(V3 dir); // return color along direction dir
		// this is done by projecting dir with each camera, and exactly
		//		one of them will have the projected point in front of the head
		//		and inside the image frame
	void SetPPCI(PPC2 *ppc, FrameBuffer *fb);
		// set conventional image with ppc camera and fb pixels from "this"
		//		cubemap
};