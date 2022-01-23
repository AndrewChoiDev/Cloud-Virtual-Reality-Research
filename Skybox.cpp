#include "stdafx.h"
#include "Skybox.h"
#include "V3.h"
#include "ppc.h"
#include <sstream>


CubeMapFaceDescriptor faceDescriptors[6] = {
	{ 
		V3(0.f, 0.f, 1.f), V3(0.f, 1.f, 0.f), 
		"back",
		{V3(1, -1, 1), V3(-1, -1, 1), V3(1, 1, 1), V3(-1, 1, 1)}
	},
	{ 
		V3(0.f, 0.f, -1.f), V3(0.f, 1.f, 0.f), 
		"front",
		{V3(-1, -1, -1), V3(1, -1, -1), V3(-1, 1, -1), V3(1, 1, -1)},
	},

	{ 
		V3(-1.f, 0.f, 0.f), V3(0.f, 1.f, 0.f), 
		"left",
		{V3(-1, -1, 1), V3(-1, -1, -1), V3(-1, 1, 1), V3(-1, 1, -1)}
	},
	{ 
		V3(1.f, 0.f, 0.f), V3(0.f, 1.f, 0.f), 
		"right",
		{V3(1, -1, -1), V3(1, -1, 1), V3(1, 1, -1), V3(1, 1, 1)}
	},

	{ 
		V3(0.f, 1.f, 0.f), V3(0.f, 0.f, 1.f), 
		"top",
		{V3(-1, 1, -1), V3(1, 1, -1), V3(-1, 1, 1), V3(1, 1, 1)}
	},
	{ 
		V3(0.f, -1.f, 0.f), V3(0.f, 0.f, -1.f), 
		"bottom",
		{V3(-1, -1, 1), V3(1, -1, 1), V3(-1, -1, -1), V3(1, -1, -1)}
	}
};

void Skybox::renderAndSaveImages(PPC* ppc)
{
	PPC ppc0(*ppc);

	auto center = ppc->C;

	auto skyboxResWidth = 512;

	*ppc = PPC(90.f, skyboxResWidth, skyboxResWidth);

	auto tempRender = FrameBuffer(50, 50, skyboxResWidth, skyboxResWidth, 25);
	tempRender.isHW = 1;
	tempRender.show();


	for (int i = 0; i < 6; i++) {
		auto faceDesc = faceDescriptors[i];

		std::stringstream ss;
		ss << "mydbg/skybox_far_region/sb_" << faceDesc.name << ".tiff";

		ppc->SetPose(center, center + faceDesc.direction, faceDesc.up);

		tempRender.redraw();
		Fl::check();
		tempRender.SaveAsTiff((char*)(ss.str().c_str()));
	}
	tempRender.hide();

	*ppc = ppc0;
}

void Skybox::renderFaces(V3 camPos)
{
	renderFace(camPos, faceDescriptors[0]);
	renderFace(camPos, faceDescriptors[1]);
	renderFace(camPos, faceDescriptors[2]);
	renderFace(camPos, faceDescriptors[3]);
	renderFace(camPos, faceDescriptors[4]);
	renderFace(camPos, faceDescriptors[5]);
}

float texCoords[8] = {
	0, 0,
	1, 0,
	0, 1,
	1, 1
};

int indices[6] = {
	0, 1, 2,
	1, 2, 3
};

unsigned int setupTexture(void* pixels, int width, int height) {
	GLuint texture_id;

	glGenTextures(1, &texture_id);
	glBindTexture(GL_TEXTURE_2D, texture_id);

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);


	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, pixels);

	glBindTexture(GL_TEXTURE_2D, 0);

	return texture_id;
}

unsigned int setupFaceTexture() {
	std::stringstream ss;
	ss << "mydbg/skybox_far_region/sb_" << faceDescriptors[1].name << ".tiff";

	FrameBuffer tempBuffer(1, 1, 512, 512, 40);

	tempBuffer.LoadTiff((char*)(ss.str().c_str()));
	//tempBuffer.SaveAsTiff("mydbg/skybox_far_region/round_trip.tiff");

	return setupTexture(tempBuffer.pix, tempBuffer.w, tempBuffer.h);
}

void Skybox::renderFace(V3 camPos, CubeMapFaceDescriptor faceDescr)
{
	//glActiveTexture(GL_TEXTURE0);
	//glClientActiveTexture(GL_TEXTURE0);
	//glEnable(GL_TEXTURE_2D);

	//auto texID = setupFaceTexture();
	//glBindTexture(GL_TEXTURE_2D, texID);

	glDisable(GL_DEPTH_TEST);

	auto trisN = 2;

	V3 verts[4];
	memcpy(verts, faceDescr.verts, sizeof(V3) * 4);

	for (int i = 0; i < 4; i++) {
		verts[i] = verts[i] + camPos;
	}

	auto tris = indices;

	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_COLOR_ARRAY);
	//glEnableClientState(GL_TEXTURE_COORD_ARRAY);

	glVertexPointer(3, GL_FLOAT, 0, (float*)verts);

	V3 debugColors[4];
	for (int i = 0; i < 4; i++) {
		debugColors[i] = V3(texCoords[i*2], texCoords[i*2+1], 0);
	}
	glColorPointer(3, GL_FLOAT, 0, (float*)debugColors);

	//glTexCoordPointer(2, GL_FLOAT, 0, texCoords);

	glDrawElements(GL_TRIANGLES, 3 * trisN, GL_UNSIGNED_INT, tris);
	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_COLOR_ARRAY);
	//glDisableClientState(GL_TEXTURE_COORD_ARRAY);

	glEnable(GL_DEPTH_TEST);
	
	//glDisable(GL_TEXTURE_2D);
}
