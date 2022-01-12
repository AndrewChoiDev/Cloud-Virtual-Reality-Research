#include "stdafx.h"

#include "framebuffer.h"
#include "math.h"
#include <iostream>
#include "scene.h"

#include <tiffio.h>


using namespace std;

FrameBuffer::FrameBuffer(int u0, int v0,
	int _w, int _h, unsigned int _id) : Fl_Gl_Window(u0, v0, _w, _h, 0) {

	isHW = 0;
	w = _w;
	h = _h;
	pix = new unsigned int[w*h];
	zb = new float[w*h];
	xyz = new V3[w*h];
	trID = new unsigned int[w*h];

	tstep = 1.0f;
	rstep = 1.0f;

}

void FrameBuffer::ClearItemBuffer() {

	for (int uv = 0; uv < w*h; uv++)
		trID[uv] = 0xFFFFFF;

}

void FrameBuffer::ClearZB() {

	for (int uv = 0; uv < w*h; uv++)
		zb[uv] = 0.0f;

}

void FrameBuffer::draw() {

	if (!isHW) {
		glDrawPixels(w, h, GL_RGBA, GL_UNSIGNED_BYTE, pix);
		return;
	}
	if (isHW == 1) {
		scene->RenderHW();
		glReadPixels(0, 0, w, h, GL_RGBA, GL_UNSIGNED_BYTE, pix);
		//cerr << "pixels read!" << endl;
		return;
	}

	scene->PerSessionHWInit();
	
	// for HW and GPU framebuffers
	if (isHW == 1) {
		scene->RenderHW();
		if (scene->smSetup) {
			scene->SetupSM();
		}
	}
	else {
		scene->RenderGPU();
	}

	glReadPixels(0, 0, w, h, GL_RGBA, GL_UNSIGNED_BYTE, pix);
	cerr << "pixels read!" << endl;

}

int FrameBuffer::handle(int event) {

	switch (event)
	{
	case FL_KEYBOARD: {
		KeyboardHandle();
		return 0;
	}
	case FL_LEFT_MOUSE: {
		cerr << "INFO: left mouse click" << endl << "create node" << endl;
		int u = Fl::event_x();
		int v = Fl::event_y();
		if (u < 0 || u > w - 1 || v < 0 || v > h - 1)
			return 0;
		V3 pP(.5f + (float)u, .5f + (float)v, GetZ(u, v));
		V3 P = scene->ppc->UnProject(pP);

		scene->streetLinesGraphSystem.handleMousePressPoint(P);

		//scene->path->AddGroundNode(P);

		scene->Render();
		return 0;
	}
	case FL_MOVE: {
		int u = Fl::event_x();
		int v = Fl::event_y();
		if (u < 0 || u > w - 1 || v < 0 || v > h - 1)
			return 0;
		V3 cv; cv.SetFromColor(Get(u, v));
		V3 rayDir;
		rayDir = scene->ppc->c + scene->ppc->a*(float)u + scene->ppc->b*(float)v;
		rayDir = rayDir.Normalized();
		V3 pP(.5f+(float)u, .5f+(float)v, GetZ(u, v));
		V3 P = scene->ppc->UnProject(pP);
#if 0
		cerr << u << " " << v << ", RGB: " << cv[0] << " " <<
			cv[1] << " " << cv[2] << ", ray: " << rayDir[0] << " "
			<< rayDir[1] << " " << rayDir[2] << "\r";
#endif
		cerr << u << " " << v << ", RGB: " << cv[0] << " " <<
			cv[1] << " " << cv[2] << ", 3D point: " << P[0] << " "
			<< P[1] << " " << P[2] << "\r";
		return 0;
	}
	default:
		return 0;
	}
	return 0;
}

void FrameBuffer::KeyboardHandle() {
	int key = Fl::event_key();
	int needRender = 0;
	switch (key) {
	case FL_Left: {
		scene->ppc->PanLeftRight(rstep);
		cerr << "INFO: pan left" << endl;
		needRender = 1;
		break;
	}
	case FL_Right: {
		scene->ppc->PanLeftRight(-rstep);
		cerr << "INFO: pan right" << endl;
		needRender = 1;
		break;
	}
	case FL_Up: {
		scene->ppc->TiltUpDown(rstep);
		cerr << "INFO: tilt up" << endl;
		needRender = 1;
		break;
	}
	case FL_Down: {
		scene->ppc->TiltUpDown(-rstep);
		cerr << "INFO: tilt down" << endl;
		needRender = 1;
		break;
	}
	case ',': {
		scene->ppc->RollLeftRight(-rstep);
		needRender = 1;
		cerr << "INFO: roll left" << endl;
		break;
	}
	case '.': {
		scene->ppc->RollLeftRight(+rstep);
		cerr << "INFO: roll right" << endl;
		needRender = 1;
		break;
	}
	case 'a': {
		scene->ppc->TranslateRightLeft(-tstep);
		cerr << "INFO: translate left" << endl;
		needRender = 1;
		break;
	}
	case 'd': {
		scene->ppc->TranslateRightLeft(+tstep);
		cerr << "INFO: translate right" << endl;
		needRender = 1;
		break;
	}
	case 'w': {
		scene->ppc->TranslateFrontBack(tstep);
		cerr << "INFO: translate forward" << endl;
		needRender = 1;
		break;
	}
	case 's': {
		scene->ppc->TranslateFrontBack(-tstep);
		cerr << "INFO: translate backward" << endl;
		needRender = 1;
		break;
	}
	case 'r': {
		scene->ppc->TranslateUpDown(tstep);
		cerr << "INFO: translate up" << endl;
		needRender = 1;
		break;
	}
	case 'f': {
		scene->ppc->TranslateUpDown(-tstep);
		cerr << "INFO: translate down" << endl;
		needRender = 1;
		break;
	}
	case '4': {
		scene->streetLinesGraphSystem.path.Save("mydbg/path_vd.txt");
		cerr << "INFO: Path Saved" << endl;
		break;
	}
	case '5': {
		scene->streetLinesGraphSystem.path.Load("mydbg/path_vd.txt");
		cerr << "INFO: Path Loaded" << endl;
		scene->Render();
		break;
	}
	case '6': {
		//scene->PlaybackPathHWSideBySide(30.0f);
		//scene->PlaybackPathHWOffsets(30.0f);
		//scene->tmeshes[12].setVisibleTrianglesHWFrameBuffer(scene->hwfb);
		//scene->tmeshes[12] = scene->tmeshes[12].constructVisibleMesh();
		//cerr << "INFO: Visible Triangles Set From FrameBuffer" << endl;
		//needRender = 1;
		//scene->PlaybackPathChunks(30.0f);
		
		//scene->PlaybackPath(30.0f);
		scene->PlaybackDeltaPathChunks(30.0f);
		break;
	}
	case '7': {
		//scene->streetLinesGraphSystem.
		break;
	}
	case '-':
		tstep /= 1.1f;
		cerr << endl << tstep << endl;
		break;
	case '=':
		tstep *= 1.1f;
		cerr << endl << tstep << endl;
		break;
	case 'p':
		scene->ppc->Save("mydbg/view.txt");
		break;
	//default:
	//	cerr << "INFO: do not understand keypress" << endl;
	//	return;
	}
	scene->streetLinesGraphSystem.handleKey(key, needRender);

//	scene->ResampleERI(scene->eri, scene->fb, scene->ppc);
	if (needRender)
		scene->Render();

}


void FrameBuffer::Set(int u, int v, unsigned int color) {

	if (u < 0 || u > w - 1 || v < 0 || v > h - 1)
		return;

	pix[(h - 1 - v)*w + u] = color;

}

void FrameBuffer::SetXYZ(int u, int v, V3 p) {

	if (u < 0 || u > w - 1 || v < 0 || v > h - 1)
		return;

	xyz[(h - 1 - v)*w + u] = p;

}



void FrameBuffer::SetZB(int u, int v, float z) {

	if (u < 0 || u > w - 1 || v < 0 || v > h - 1)
		return;

	zb[(h - 1 - v)*w + u] = z;

}

void FrameBuffer::SetBGR(unsigned int bgr) {

	for (int uv = 0; uv < w*h; uv++)
		pix[uv] = bgr;

}

// load a tiff image to pixel buffer
void FrameBuffer::LoadTiff(char* fname) {
	TIFF* in = TIFFOpen(fname, "r");
	if (in == NULL) {
		cerr << fname << " could not be opened" << endl;
		return;
	}

	int width, height;
	TIFFGetField(in, TIFFTAG_IMAGEWIDTH, &width);
	TIFFGetField(in, TIFFTAG_IMAGELENGTH, &height);
	if (w != width || h != height) {
		w = width;
		h = height;
		delete[] pix;
		if (zb)
			delete[] zb;
		pix = new unsigned int[w*h];
		zb = new float[w*h];
		size(w, h);
		glFlush();
		glFlush();
	}

	if (TIFFReadRGBAImage(in, w, h, pix, 0) == 0) {
		cerr << "failed to load " << fname << endl;
	}

	TIFFClose(in);
}

// save as tiff image
void FrameBuffer::SaveAsTiff(char *fname) {

	TIFF* out = TIFFOpen(fname, "w");

	if (out == NULL) {
		cerr << fname << " could not be opened" << endl;
		return;
	}

	TIFFSetField(out, TIFFTAG_IMAGEWIDTH, w);
	TIFFSetField(out, TIFFTAG_IMAGELENGTH, h);
	TIFFSetField(out, TIFFTAG_SAMPLESPERPIXEL, 4);
	TIFFSetField(out, TIFFTAG_BITSPERSAMPLE, 8);
	TIFFSetField(out, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
	TIFFSetField(out, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
	TIFFSetField(out, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);

	for (uint32 row = 0; row < (unsigned int)h; row++) {
		TIFFWriteScanline(out, &pix[(h - row - 1) * w], row);
	}

	TIFFClose(out);
}

void FrameBuffer::SetChecker(unsigned int col0, unsigned int col1, int csize) {

	for (int v = 0; v < h; v++) {
		for (int u = 0; u < w; u++) {
			int cv = v / csize;
			int cu = u / csize;
			if ((cu + cv) % 2) {
				Set(u, v, col0);
			}
			else {
				Set(u, v, col1);
			}
		}
	}

}

int FrameBuffer::Farther(int u, int v, float currz) {

	if (u < 0 || u > w - 1 || v < 0 || v > h - 1)
		return 1;
	int uv = (h - 1 - v)*w + u;
	if (currz < zb[uv])
		return 1;
	zb[uv] = currz;
	return 0;

}

void FrameBuffer::Draw2DSegment(V3 p0, V3 c0, V3 p1, V3 c1) {

	float du = fabsf((p0 - p1)[0]);
	float dv = fabsf((p0 - p1)[1]);
	int stepsN;
	if (du < dv) {
		stepsN = 1+(int)dv;
	}
	else {
		stepsN = 1+(int)du;
	}
	for (int i = 0; i <= stepsN; i++) {
		V3 cp, cc;
		cp = p0 + (p1 - p0) * (float)i / (float)stepsN;
		// cp[2] depth (one over w) at current pixel
		int u = (int)cp[0], v = (int)cp[1];
		if (Farther(u, v, cp[2]))
			continue;
		cc = c0 + (c1 - c0) * (float)i / (float)stepsN;
		Set(u, v, cc.GetColor());
	}

}

void FrameBuffer::DrawSquarePoint(float uf, float vf, int psize, unsigned int color) {

	int u = (int)uf;
	int v = (int)vf;
	for (int cv = v - psize / 2; cv <= v + psize / 2; cv++)
		for (int cu = u - psize / 2; cu <= u + psize / 2; cu++)
			Set(cu, cv, color);
}

void FrameBuffer::DrawDisk(float uf, float vf, float r, unsigned int color) {


	for (int cv = (int)((int)vf - r - 1); cv < (int)((int)vf + r + 1); cv++) {
		for (int cu = (int)((int)uf - r - 1); cu < (int)((int)uf + r + 1); cu++) {
			V3 pixCenter(.5f + (float)cu, .5f + (float)cv, 0.0f);
			V3 cCenter(uf, vf, 0.0f);
			if ((pixCenter - cCenter).Length() > r)
				continue;
			Set(cu, cv, color);
		}
	}
}



void FrameBuffer::Draw3DSegment(V3 P0, V3 P1, PPC *ppc, V3 c0, V3 c1) {

	V3 p0, p1;
	if (!ppc->Project(P0, p0))
		return;
	if (!ppc->Project(P1, p1))
		return;

	Draw2DSegment(p0, c0, p1, c1);

}

void FrameBuffer::Draw3DPoint(V3 P, PPC *ppc, unsigned int color, int psize) {


	V3 pP;
	if (!ppc->Project(P, pP))
		return;

	int u = (int)pP[0];
	int v = (int)pP[1];
	for (int cv = v - psize / 2; cv <= v + psize / 2; cv++) {
		for (int cu = u - psize / 2; cu <= u + psize / 2; cu++) {
			if (Farther(cu, cv, pP[2]))
				continue;
			Set(cu, cv, color);
		}
	}

}

unsigned int FrameBuffer::Get(int u, int v) {

	return pix[(h - 1 - v)*w + u];

}

float FrameBuffer::GetZ(int u, int v) {

	return zb[(h - 1 - v)*w + u];

}

unsigned int FrameBuffer::Lookup(float uf, float vf) {

	if (uf < 0.0f || uf >= (float)w ||
		vf < 0.0f || vf >= (float)h)
		return 0xFFFF00FF;
	int u = (int)uf;
	int v = (int)vf;
	return Get(u, v);
}

void FrameBuffer::GeneralizedRotation(PPC *ppc0, FrameBuffer *fb1, PPC *ppc1) {

	for (int v = 0; v < fb1->h; v++) {
		for (int u = 0; u < fb1->w; u++) {
			V3 P = ppc1->UnProject(V3(.5f + (float)u, .5f + (float)v, 1.0f));
			V3 pp;
			if (!ppc0->Project(P, pp))
				continue;
			fb1->Set(u, v, Lookup(pp[0], pp[1]));
		}
	}

}


void FrameBuffer::SetFromImage(FrameBuffer *fb) {

	for (int uv = 0; uv < w*h; uv++)
		pix[uv] = fb->pix[uv];


}

void FrameBuffer::DrawDiagonalCrossPadded(float uf, float vf, int csize,
	unsigned int color) {

	DrawDiagonalCross(uf, vf, csize, color);
	for (int v = (int)vf - csize / 2; v <= (int)vf + csize / 2; v++)
		for (int u = (int)uf - csize / 2; u <= (int)uf + csize / 2; u++)
			if (Get(u, v) != color)
				Set(u, v, 0xFF000000);

}

void FrameBuffer::DrawDiagonalCross(float uf, float vf, int csize,
	unsigned int color) {

	int u = (int)uf;
	int v = (int)vf;
	if (u - csize / 2 < 0 || u + csize / 2 > w - 1 || v - csize / 2 < 0 || v + csize / 2 > h - 1)
		return;
	for (int ci = -csize / 2; ci <= csize / 2; ci++) {
		Set(u + ci, v + ci, color);
		Set(u - ci, v + ci, color);
	}

}

void FrameBuffer::DrawStraightCrossPadded(float uf, float vf, int csize,
	unsigned int color) {

	DrawStraightCross(uf, vf, csize, color);
	for (int v = (int)vf - csize / 2; v <= (int)vf + csize / 2; v++)
		for (int u = (int)uf - csize / 2; u <= (int)uf + csize / 2; u++)
			if (Get(u, v) != color)
				Set(u, v, 0xFF000000);

}

void FrameBuffer::DrawStraightCross(float uf, float vf, int csize,
	unsigned int color) {

	int u = (int)uf;
	int v = (int)vf;
	if (u - csize / 2 < 0 || u + csize / 2 > w - 1 || v - csize / 2 < 0 || v + csize / 2 > h - 1)
		return;
	for (int ci = -csize / 2; ci <= csize / 2; ci++) {
		Set(u, v + ci, color);
		Set(u + ci, v, color);
	}

}


void FrameBuffer::SetItemBuffer(int u, int v, unsigned int tri) {

	trID[(h - 1 - v)*w + u] = tri;

}
