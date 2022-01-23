#pragma once

#include "ppc.h"

struct CubeMapFaceDescriptor {
	V3 direction;
	V3 up;
	std::string name;
	V3 verts[4];
};

class Skybox
{
public: 
	static void renderAndSaveImages(PPC* ppc);
	static void renderFaces(V3 camPos);

private:
	static void renderFace(V3 camPos, CubeMapFaceDescriptor faceDescr);
};

