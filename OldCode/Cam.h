#pragma once

#include "v3.h"

// this is the declaration of the planar pinhole camera class
// this is the code we wrote in class; there is a similar declaration in ppc.h
// see PHC.pdf

class PPC2 {
public:
	V3 a, b, c, C;
	int w, h;
	PPC2(float hfov, int _w, int _h); // construct a planar pinhole camera with 
							//	given horizontal field of view and resolution
	int Project(V3 P, V3& Pp); // given a 3D point P, return projected point Pp
								// Pp[0] is u or x screen coordinate
								// Pp[1] is v pr y screen coordinate
								// Pp[2] is a quantity proportional with 1/w (1/z)
								// returns 0 if P behind head
	int InsideImage(V3 Pp); // returns 1 if projection is inside image frame
	void Pan(float pan); // pans camera <pan> degrees, i.e. rotates left
	void Tilt(float tilt); // tilts camera <tilt> degrees, i.e. rotates up
	PPC2() {}; // default constructor, needed to be able to allocate array of cameras
	V3 GetRay(float uf, float vf); // unprojection, or conversion from pixel coordinates
				// uf, vf to ray direction (does not normalize)
};