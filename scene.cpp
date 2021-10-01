#include "stdafx.h"

#include "scene.h"

#include "V3.h"
#include "M33.h"
#include "ppc.h"
#include "TMesh.h"
#include "AABB.h"

#include "RandomCamera.h"

#include "Path.h"

Scene *scene;

using namespace std;

#include <iostream>
#include <fstream>

Scene::Scene() {

	renderHW = 0;
	cgi = 0;
	soi = 0;

	hwPerSessionInit = 0;

	gui = new GUI();
	gui->show();

	int u0 = 20;
	int v0 = 100;
	int h = 480;
	int w = 640;

	fb = new FrameBuffer(u0, v0, w, h, 0);
	fb->label("SW 1");
	fb->show();
	fb->redraw();

	fb3 = new FrameBuffer(u0+w+30, v0, w, h, 0);
	fb3->label("SW 3");
//	fb3->show();
	fb3->redraw();


	gui->uiw->position(u0, v0 + h + 50);

	float hfov = 55.0f;
	ppc = new PPC(hfov, fb->w, fb->h);
	ppc3 = new PPC(hfov, fb3->w, fb3->h);

	tmeshesN = 13;
	tmeshes = new TMesh[tmeshesN];

	V3 cc(0.0f, 0.0f, -100.0f);
	float sideLength = 60.0f;
	tmeshes[0].SetToCube(cc, sideLength, 0xFF0000FF, 0xFF000000);
	tmeshes[0].onFlag = 0;

	tmeshes[1].LoadBin("geometry/teapot1K.bin");
//	tmeshes[1].LoadBin("geometry/teapot57K.bin");
	tmeshes[1].SetCenter(V3(0.0f, 0.0f, -140.0f));
	tmeshes[1].onFlag = 0;
	tmeshes[1].reflectorFlag = 1;

	V3 qverts[4];
	qverts[0] = V3(-20.0f, 5.0f, -20.0f);
	qverts[1] = V3(100.0f, 5.0f, -200.0f);
	qverts[2] = V3(100.0f, -5.0f, -200.0f);
	qverts[3] = V3(-20.0f, -5.0f, -20.0f);
	V3 qcolors[4] = { V3(1.0f, 0.0f, 0.0f), V3(0.0f, 0.0f, 0.0f), V3(0.0f, 0.0f, 0.0f),
		V3(1.0f, 0.0f, 0.0f) };
	tmeshes[2].SetQuad(qverts, qcolors);
	tmeshes[2].Translate(V3(0.0f, 0.0f, -20.0f));
	tmeshes[2].Translate(V3(0.0f, 6.0f, 0.0f));
	tmeshes[2].onFlag = 0;

	tmeshes[3].SetQuad(qverts, qcolors);
	tmeshes[3].Translate(V3(0.0f, 0.0f, -20.0f));
	tmeshes[3].Translate(V3(0.0f, -6.0f, 0.0f));
	tmeshes[3].msiFlag = 0;

	tmeshes[4].SetQuad(qverts, qcolors);
	tmeshes[4].Translate(V3(0.0f, 0.0f, -22.0f));
	tmeshes[4].colors[0] = tmeshes[4].colors[3] = tmeshes[4].colors[1] = 
		tmeshes[4].colors[2] = V3(0.0f, 1.0f, 0.0f);
	tmeshes[4].verts[0] = (tmeshes[4].verts[0] + tmeshes[4].verts[1]) / 2.0f;
	tmeshes[4].verts[3] = (tmeshes[4].verts[2] + tmeshes[4].verts[3]) / 2.0f;

	float qz = -1.0f;
	float qs = 400.0f;
	qverts[0] = V3(-qs, qs, qz);
	qverts[1] = V3(-qs, -qs, qz);
	qverts[2] = V3(qs, -qs, qz);
	qverts[3] = V3(qs, qs, qz);
	qcolors[0] = V3(0.0f, 0.0f, 1.0f);
	qcolors[1] = V3(0.0f, 1.0f, 1.0f);
	qcolors[2] = V3(1.0f, 1.0f, 1.0f);
	qcolors[3] = V3(1.0f, 0.0f, 1.0f);
	int texw = 1024;
	FrameBuffer *texture = new FrameBuffer(30, 30, texw, texw, 0);
//	texture->show();
	texture->label("Texture");
	texture->SetChecker(0xFF000000, 0xFFFFFFFF, 16);
	tmeshes[5].texture = texture;
	tmeshes[5].SetQuad(qverts, qcolors);
//	tmeshes[5].Translate(V3(0.0f, 0.0f, -400.0f));

	AABB aabb;
	tmeshes[1].SetAABB(aabb);
	aabb.corners[0][1] += 3.0f;

	// CONTINUE HERE
	qverts[0] = V3(aabb.corners[0][0], aabb.corners[0][1], aabb.corners[0][2]);
	qverts[1] = V3(aabb.corners[0][0], aabb.corners[0][1], aabb.corners[1][2]);
	qverts[2] = V3(aabb.corners[1][0], aabb.corners[0][1], aabb.corners[1][2]);
	qverts[3] = V3(aabb.corners[1][0], aabb.corners[0][1], aabb.corners[0][2]);
//	qcolors[0] = qcolors[1] = qcolors[2] = qcolors[3] = V3(1.0f, 0.0f, 0.0f);
	tmeshes[6].SetQuad(qverts, qcolors);
	tmeshes[6].texture = new FrameBuffer(100, 100, 64, 64, 0);
	tmeshes[6].texture->SetChecker(0xFFFFFFFF, 0xFF000000, 8);
	tmeshes[6].Scale(2.0f);

//	ppc->SetPose(ppc->C + V3(0.0f, 50.0f, 0.0f), tmeshes[1].GetCenter(), V3(0.0f, 1.0f, 0.0f));
	ppc3->SetPose(V3(0.0f, 40.0f, 50.0f), tmeshes[2].GetCenter(), V3(0.0f, 1.0f, 0.0f));

	tmeshes[3].onFlag = 0;
	tmeshes[4].onFlag = 0;
	tmeshes[5].onFlag = 0;
	tmeshes[6].onFlag = 0;
	tmeshes[7].onFlag = 0;
	tmeshes[8].onFlag = 0;

	tmeshes[9].onFlag = 1;
	qs = 20.0f;
	qz = -100.0f;
	qverts[0] = V3(-qs*2.0f, qs, qz);
	qverts[1] = V3(-qs*2.0f, -qs, qz);
	qverts[2] = V3(qs*2.0f, -qs, qz);
	qverts[3] = V3(qs*2.0f, qs, qz);
	qcolors[0] = V3(0.0f, 0.0f, 0.0f);
	qcolors[1] = V3(0.0f, 1.0f, 0.0f);
	qcolors[2] = V3(1.0f, 1.0f, 0.0f);
	qcolors[3] = V3(1.0f, 0.0f, 0.0f);
	tmeshes[9].SetQuad(qverts, qcolors);

	tmeshes[9].onFlag = 1;
	qs = 5.0f;
	qz = -50.0f;
	qverts[0] = V3(-qs*2.0f, qs, qz);
	qverts[1] = V3(-qs*2.0f, -qs, qz);
	qverts[2] = V3(qs*2.0f, -qs, qz);
	qverts[3] = V3(qs*2.0f, qs, qz);
	qcolors[0] = V3(0.0f, 0.0f, 0.0f);
	qcolors[1] = V3(0.0f, 1.0f, 0.0f);
	qcolors[2] = V3(1.0f, 1.0f, 0.0f);
	qcolors[3] = V3(1.0f, 0.0f, 0.0f);
	tmeshes[10].SetQuad(qverts, qcolors);

	tmeshes[11].LoadBin("geometry/auditorium.bin");
	tmeshes[11].Translate(tmeshes[11].GetCenter()*-1.0f);
	tmeshes[11].onFlag = false;

	vf = 20.0f;
	L = V3(ppc->C);
	ka = 0.2f;

	Render();

	hwfb = new FrameBuffer(u0 + w + 30, v0, w, h, 0);
	hwfb->label("HW fb");
	hwfb->isHW = 1;
	hwfb->show();
	hwfb->redraw();

	gpufb = new FrameBuffer(u0 + w + 30, v0 + h + 50, w, h, 0);
	gpufb->label("GPU Framebuffer");
	gpufb->isHW = 2;
//	gpufb->show();
	gpufb->redraw();

	smSetup = 0;
	EdgeVRSetup();
	return;

	SetupTeapot();

	return;

	hfov = 95.0f;
	ppc = new PPC(hfov, fb->w, fb->h);
	eri = new FrameBuffer(20 + 20 + fb->w, 100, 10, 10, 0);
	eri->LoadTiff("mydbg/ERPanorama.tiff");
//	eri->LoadTiff("mydbg/snowERPanorama.tif");
	//	eri->LoadTiff("mydbg/louvreERPanorama.tif");
	eri->label("Equirectangular panorama");
	cerr << eri->w << " " << eri->h << endl;
	eri->show();
	ppc->SetPose(ppc->C, ppc->C + V3(0.9434f, -0.3212f, 0.0823f), V3(0.0f, 1.0f, 0.0f));
	dsc = V3(34.8271f, -4.0f, 2.06099f);
	dsr = 1.0f;
	ResampleERI(eri, fb, ppc);

	TrackImageSetup();

}

void Scene::Render() {

	Render(fb, ppc);
	if (renderHW)
		hwfb->redraw();
	return;


//	fb->Draw3DPoint(L, ppc, 0xFF00FFFF, 7);
	return;

	fb3->ClearZB();
	fb3->SetBGR(0xFFFFFFFF);
	ppc->Visualize(fb3, ppc3, fb);
	fb3->redraw();
	return;

	Render(fb3, ppc3);

	ppc->Visualize(fb3, ppc3, vf);
	ppc->Visualize(fb3, ppc3, vf, fb);
	fb3->redraw();

}


void Scene::RenderHW() {

	// initializations (could be done once per session)
	glEnable(GL_DEPTH_TEST);

	// clear buffers
	glClearColor(0.0f, 0.0f, 0.5f, 1.0f);
	glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);

	// set the view desired by the application (the user)
	ppc->SetIntrinsicsHW();
	ppc->SetExtrinsicsHW();

	// draw the actual geometry
	for (int tmi = 0; tmi < tmeshesN; tmi++) {
		if (!tmeshes[tmi].onFlag)
			continue;
		tmeshes[tmi].RenderHW();
	}

}

void Scene::RenderGPU() {

	// if the first time, call per session initialization
	if (cgi == NULL) {
		cgi = new CGInterface();
		cgi->PerSessionInit();
		soi = new ShaderOneInterface();
		soi->PerSessionInit(cgi);
	}

	// clear the framebuffer
	glClearColor(0.0, 0.0f, 0.5f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// per frame initialization
	cgi->EnableProfiles();
	soi->PerFrameInit();
	soi->BindPrograms();

	// render geometry
	for (int tmi = 0; tmi < tmeshesN; tmi++) {
		if (!tmeshes[tmi].onFlag)
			continue;
		tmeshes[tmi].RenderHW();
	}

	soi->PerFrameDisable();
	cgi->DisableProfiles();



}




void Scene::Render(FrameBuffer *rfb, PPC *rppc) {

	rfb->SetBGR(0xFFFFFFFF);
	rfb->ClearZB();
	rfb->ClearItemBuffer();

	for (int tmi = 0; tmi < tmeshesN; tmi++) {
		if (!tmeshes[tmi].onFlag)
			continue;
//		tmeshes[tmi].DrawWireFrame(rfb, rppc, 0xFF000000);
		V3 C(1.0f, 0.0f, 0.0f);
//		tmeshes[tmi].Light(C, L, ka);
		tmeshes[tmi].RenderFilled(rfb, rppc, C, L, ka);
//		tmeshes[tmi].DrawWireFrame(rfb, rppc, 0xFF000000);
//		tmeshes[tmi].RayTrace(rfb, rppc);
	}

	if (path)
		path->Render(rfb, ppc);

	rfb->redraw();

}

void Scene::SetupTeapot() {

	hwfb->hide();
	gpufb->hide();
	for (int tmi = 0; tmi < tmeshesN; tmi++) {
		tmeshes[tmi].onFlag = 0;
	}
	tmeshes[1].onFlag = 1;
	tmeshes[2].onFlag = 1;
	float rw = 100.0f, rh = 50.0f;
	tmeshes[2].SetRectangle(rw, rh);
	tmeshes[2].SetCenter(tmeshes[1].GetCenter());
	tmeshes[2].Rotate(tmeshes[2].GetCenter(), V3(0.0f, 1.0f, 0.0f), 45.0f);
	Render();

}

void Scene::ResamplePPC(FrameBuffer *ppi0, PPC *ppc0, FrameBuffer *ppi1, PPC *ppc1) {

	ppi1->SetBGR(0xFFFF00FF);
	// traverse output image
	// for each pixel (u1, v1) in output image find corresponding pixel (u0, v0) in input conventional image
	for (int v1 = 0; v1 < ppi1->h; v1++) {
		for (int u1 = 0; u1 < ppi1->w; u1++) {
			// create ray through current output pixel
			V3 ray = (ppc1->c + ppc1->a*(.5f + u1) + ppc1->b*(.5f + v1)).Normalized();
			int u0, v0;
			V3 pp;
			if (!ppc0->Project(ray, pp))
				continue;
			u0 = (int)pp[0];
			v0 = (int)pp[1];
			if (u0 < 0 || u0 >= ppi0->w || v0 < 0 || v0 >= ppi0->h)
				continue;
			ppi1->Set(u1, v1, ppi0->Get(u0, v0));
		}
	}

}


int Scene::IntersectRayWithSphere(V3 O, float r, V3 C, V3 d, V3 &P) {

	V3 qabc;
	qabc[0] = 1.0f;
	qabc[1] = ((C - O) * d)*2.0f;
	qabc[2] = (C - O)*(C - O) - r*r;
	float dis = qabc[1] * qabc[1] - 4.0f * qabc[2];
	if (dis < 0.0f)
		return 0;
	float t = (-qabc[1] - sqrtf(dis)) / 2.0f;
	if (t < 0.0f)
		return 0;
	P = C + d*t;
	return 1;

}

unsigned int Scene::LookupRayIntoERI(V3 ray, FrameBuffer *eri) {

	// compute ERC projection
	float t = asinf(ray[1]) / 3.14159265f * 180.0f;
	int v0 = (int)((1.0f - (90.0f + t) / 180.0f) * (float)eri->h);

	float q = atanf(ray[2] / ray[0]) / 3.14159265f * 180.0f;
	if (ray[0] > 0.0f) {
		if (ray[2] > 0.0f)
			q = -q + 360.0f;
		else q = -q;
	}
	else {
		if (ray[2] > 0.0f)
			q = -q + 180.0f;
		else
			q = 180.0f - q;
	}
	q = 360.0f - q;
	int u0 = (int)(q / 360.0f*(float)eri->w);
	u0 = (u0 < 0) ? 0 : u0;
	u0 = (u0 >= eri->w - 1) ? eri->w - 1 : u0;
	v0 = (v0 < 0) ? 0 : v0;
	v0 = (v0 >= eri->h - 1) ? eri->h - 1 : v0;

	return eri->Get(u0, v0);

}

V3 Scene::ReflectRay(V3 ray, V3 n) {

	V3 nn = n * -1.0f;
	V3 in = nn*(ray*nn);
	V3 it = ray - in;
	V3 rr = it - in;
	return rr;

}

int Scene::IntersetRayWithPlane(V3 rayO, V3 rayDir, V3 pPoint, V3 pNormal, V3 &Q) {

	V3 P0 = pPoint, n = pNormal, rO = rayO, rD = rayDir;
	float t = ((P0 - rO)*n) / (rD*n);

	if (t < 0.0f)
		return 0;
	
	Q = rayO + rayDir*t;
	return 1;

}


void Scene::ResampleERI(FrameBuffer *eri, FrameBuffer *ppi, PPC *ppc) {

	ppi->SetBGR(0xFF000000);

	// chrome sphere parameters
	V3 sO;
	V3 vd = (ppc->a ^ ppc->b).Normalized();
	float d = 10.0f;
	float sr = 3.0f;
	sO = ppc->C + vd*d;

	// ground plane
	float h = 4.0f; // height of center of panorama above ground plane
	float gr = 10000.0f; // radius of ground plane patch
	V3 gNormal(0.0f, 1.0f, 0.0f);
	V3 gPoint(0.0f, -h, 0.0f); // G frm figure
	// assumptions: ground plane is horizontal, and equirectangular image is north south (vertical)
	// traverse output image
	// for each pixel (u1, v1) in output image find corresponding pixel (u0, v0) in equirectangular image


#if 0
	// walking sphere (diffuse, walking on wooden path)
	float ug0 = 243.0f, vg0 = 100.0f;
	V3 gP0, gP1;
	V3 ray = (ppc->c + ppc->a*(.5f + ug0) + ppc->b*(.5f + vg0)).Normalized();
	IntersetRayWithPlane(ppc->C, ray, gPoint, gNormal, gP0);
	float ug1 = 262.0f, vg1 = 200.0f;
	ray = (ppc->c + ppc->a*(.5f + ug1) + ppc->b*(.5f + vg1)).Normalized();
	IntersetRayWithPlane(ppc->C, ray, gPoint, gNormal, gP1);
	cerr << gP0 << gP1;
	return;
#endif

	V3 sunDir(-0.3798f, 0.1502f, 0.9127f);
	for (int v1 = 0; v1 < ppi->h; v1++) {
		for (int u1 = 0; u1 < ppi->w; u1++) {
			// create ray through current output pixel
			V3 ray = (ppc->c + ppc->a*(.5f + u1) + ppc->b*(.5f + v1)).Normalized();
			// compute intersection with diffuse sphere
			float rayd = FLT_MAX; // distance along ray to nearest intersection
			V3 S;
			if (IntersectRayWithSphere(dsc, dsr, ppc->C, ray, S)) {
				sunDir = sunDir.Normalized();
				V3 n = (S - dsc).Normalized(); // sphere intersection normal
				V3 col(1.0f, 1.0f, 1.0f);
				float ka = 0.2f;
				float kd = sunDir*n;
				kd = (kd < 0.0f) ? 0.0f : kd;
				V3 litCol = col*(ka + (1.0f - ka)*kd);
				ppi->Set(u1, v1, litCol.GetColor());
				rayd = (S - ppc->C).Length();
			}
			// compute intersection with ground plane
			V3 Q;
			if (IntersetRayWithPlane(ppc->C, ray, gPoint, gNormal, Q)) {
				float currrayd = (Q - ppc->C).Length();
				if (currrayd > rayd)
					continue;
				rayd = currrayd;
				float d = (Q - gPoint).Length();
				if (d < gr) {
					V3 gpRay = Q.Normalized();
					V3 S;
					V3 colv; colv.SetFromColor(LookupRayIntoERI(gpRay, eri));
					if (IntersectRayWithSphere(dsc, dsr, Q, sunDir, S)) {
						float dtob = (Q-S).Length(); // distance to blocker
						float maxd = 2.0f;
						// in shadow
						float ks = dtob / maxd + 0.2f;
						ks = (ks > 1.0f) ? 1.0f : ks;
						colv = colv * ks;
					}
					ppi->Set(u1, v1, colv.GetColor());
					continue;
				}
				else if (true && d < gr*2.0f) {
					// in transiton region
					V3 gpRay = Q.Normalized();
					float frac = (d - gr) / gr;
					V3 interRay = (gpRay*(1.0f - frac) + ray*frac).Normalized();
					V3 interColor; interColor.SetFromColor(LookupRayIntoERI(interRay, eri));
//					interColor[2] += (1.0f - interColor[2]) / 2.0f;
//					interColor[0] /= 2.0f;
//					interColor[1] /= 2.0f;
					ppi->Set(u1, v1, interColor.GetColor());
					continue;
				}
			}
			// no intersection with the ground patch
			V3 envColor; envColor.SetFromColor(LookupRayIntoERI(ray, eri));
			envColor[0] += (1.0f - envColor[0]) / 2.0f;
			envColor[1] /= 2.0f;
			envColor[2] /= 2.0f;
			ppi->Set(u1, v1, envColor.GetColor());
			continue;

			// compute intersection with sphere
			V3 P;
			if (false && IntersectRayWithSphere(sO, sr, ppc->C, ray, P)) {
				// compute normal
				V3 n = (P - sO).Normalized();
				// compute reflected ray
				V3 rr = ReflectRay(ray, n);
				V3 reflectiveColor;
				reflectiveColor.SetFromColor(LookupRayIntoERI(rr, eri));
				V3 diffuseColor(1.0f, 0.0f, 0.0f);
				float kd = 0.5f;
				V3 finalColor = diffuseColor * kd + reflectiveColor * (1.0f - kd);
				ppi->Set(u1, v1, finalColor.GetColor());
				continue;
			}

			// else lookup eye ray into equi. image
			unsigned int eriColor = LookupRayIntoERI(ray, eri);
			ppi->Set(u1, v1, eriColor);
		}
	}

	ppi->redraw();

}

void Scene::MakeDotPattern() {

	int N = 400*10; // number of dots
	int w = 640*10, h = 480*10; // resolution of image
	float R = 23.0f; // dot radius
	float extraDistance = 15.0f;
	float minD = 2*R+extraDistance; // minimum distance between dot centers
	cerr << endl;
	cerr << "INPUT: number of dots: ";
	cin >> N;
	cerr << "INPUT: image resolution in pixels w x h: ";
	cin >> w >> h;
	cerr << "INPUT: dot radius R: ";
	cin >> R;
	cerr << "INPUT: extra distance between dots in addition to 2R: ";
	cin >> extraDistance;

	float *dots = new float[2 * N];
	int currDotsN = 0;
	cerr << endl;
	while (currDotsN < N) {
		float x = R + (float)rand() / (float) RAND_MAX * (float)(w-2*(int)R-1);
		float y = R + (float)rand() / (float) RAND_MAX * (float)(h-2*(int)R-1);
		int goodDot = true;
		for (int di = 0; di < currDotsN; di++) {
			if ((x - dots[2 * di + 0])*(x - dots[2 * di + 0]) +
				(y - dots[2 * di + 1])*(y - dots[2 * di + 1]) < minD * minD) {
				goodDot = false;
				break;
			}
		}
		if (goodDot) {
			dots[2 * currDotsN + 0] = x;
			dots[2 * currDotsN + 1] = y;
			currDotsN++;
			cerr << "INFO: dots: " << currDotsN << "        \r";
		}
	}
	cerr << endl << "INFO: done generating dots" << endl;

	FrameBuffer *dotPattern = new FrameBuffer(100, 100, w, h, 0);
	dotPattern->SetBGR(0xFFFFFFFF);
	ofstream ofs("mydbg/patternCoordinates.txt");
	if (ofs.fail()) {
		cerr << "INFO: cannot open text file" << endl;
		return;
	}
	for (int di = 0; di < N; di++) {
		dotPattern->DrawDisk(dots[2 * di + 0], dots[2 * di + 1], R, 0xFF000000);
		ofs << dots[2 * di + 0] << " " << dots[2 * di + 1] << endl;
	}
//	dotPattern->show();
	dotPattern->SaveAsTiff("mydbg/pattern.tif");
	ofs.close();
}

void Scene::EdgeVRSetup() {
	
	for (int tmi = 0; tmi < tmeshesN; tmi++)
		tmeshes[tmi].onFlag = 0;
	tmeshes[12].LoadBin("geometry/Manhattan.bin");
	tmeshes[12].onFlag = 1;
//	tmeshes[12].msiFlag = 0;
	AABB aabb;
	tmeshes[12].SetAABB(aabb);
	V3 cityCenter = tmeshes[12].GetCenter();
	float cityDiagonal = aabb.GetDiagonal();
	ppc->SetPose(cityCenter + V3(-cityDiagonal / 5.0f, cityDiagonal/5.0f, 0.0f), cityCenter,
		V3(0.0f, 1.0f, 0.0f));
	ppc->Load("mydbg/view.txt");
	fb->tstep = cityDiagonal / 1000.0f;
	fb->rstep = 3.0f;
	hwfb->show();
	renderHW = 1;
	Render();
	path = new Path;
	path->speed = 50.0f;
	path->upGuidance = V3(0.0f, 1.0f, 0.0f);
	path->height = 2.5f;
}

void Scene::PlaybackPath(float fps) {

	float time = path->GetTotalTime();
	PPC ppc0(*ppc);
	for (int fi = 0; fi < (int)(time * fps); fi++) {
		path->SetCamera(ppc, ppc, (float)fi / fps);
		fb->Draw3DPoint(ppc->C, &ppc0, 0xFF00FF00, 7);
		fb->redraw();
		hwfb->redraw();
		Fl::check();
	}
	*ppc = ppc0;

}

void Scene::CollectVisibleTrianglesOnPath(float t0, float t1, float fps) {

	float time = path->GetTotalTime();
	PPC ppc0(*ppc);
	float t;
	tmeshes[12].ClearVisibleTriangles();
	for (t = t0; t < t1; t += 1.0f/fps) {
		path->SetCamera(ppc, ppc, t);
		scene->Render(fb, ppc);
		tmeshes[12].AddVisibleTriangles(fb);
		fb->redraw();
		Fl::check();
	}
	cerr << endl << "INFO: visible triangles found: " 
		<< tmeshes[12].CountVisibleTriangles() << endl;
	*ppc = ppc0;

}


void Scene::DBG() {

	{
		PlaybackPath(30.0f);
		return;
		CollectVisibleTrianglesOnPath(0.0f, 2.0f, 30.0f);
		return;
	}

	{
		TransparentDisplayPlanarProxy();
		return;

	}
	{
		TransparentDisplay();
		return;

	}
	{
		ImageTracking();
		return;
	}

	{
		MakeDotPattern();
		return;
	}

	{

		V3 gP0(34.8271f, -4.0f, 2.06099f);
		V3 gP1(6.55237f, -4, 0.985015f);
		int fN = 100;
		float va = 1.0f;
		for (int fi = 0; fi < fN; fi++) {
			V3 gP = gP0 + (gP1 - gP0)*(float)fi / (float)(fN - 1);
			V3 voff(0.0f, va*(float)fi / (float)(fN - 1), 0.0f);
			dsc = gP;
			ResampleERI(eri, fb, ppc);
			Fl::check();
		}
		return;
	}

	{
		ppc->C = V3(0.0f, 0.0f, 0.0f);
		ResampleERI(eri, fb, ppc);
		return;

	}

	{

		int fN = 3600;
		for (int fi = 0; fi < fN; fi++) {
			ResampleERI(eri, fb, ppc);
			Fl::check();
			ppc->PanLeftRight(0.3f);
//			ppc->TiltUpDown(0.1f);
//			break;
		}
		return;
		eri->hide();
		fb->label("Input PPC image");
		fb->redraw();
		FrameBuffer *ppi1 = new FrameBuffer(20 + 20 + fb->w, 100, fb->w, fb->h, 0);
		ppi1->label("User image from PPC image");
		ppi1->show();
		PPC *ppc1 = new PPC(90.0f, ppi1->w, ppi1->h);
		for (int fi = 0; fi < 900; fi++) {
			ResamplePPC(fb, ppc, ppi1, ppc1);
			ppi1->redraw();
			ppc1->PanLeftRight(0.1f);
			ppc1->TiltUpDown(0.1f);
			Fl::check();
		}
		return;
	}


	{
		RCSetup();
		return;
	}



	{

		fb->ClearZB();
		fb->SetBGR(0xFFFFFFFF);
		V3 pP;
		for (int vi = 0; vi < tmeshes[2].vertsN; vi++) {
			if (!ppc->Project(tmeshes[2].verts[vi], pP))
				continue;
			fb->DrawSquarePoint(pP[0], pP[1], 17, 0xFF00FF00);
		}
		int n = 10;
		for (int tri = 0; tri < tmeshes[2].trisN; tri++) {
			for (int ei = 0; ei < 3; ei++) {
				V3 V0 = tmeshes[2].verts[tmeshes[2].tris[tri * 3 + ei]];
				V3 V1 = tmeshes[2].verts[tmeshes[2].tris[tri * 3 + (ei+1)%3]];
				for (int pi = 0; pi < n; pi++) {
					V3 V = V0 + (V1 - V0)*(float)pi / (float)(n-1);
					if (!ppc->Project(V, pP))
						continue;
					fb->DrawSquarePoint(pP[0], pP[1], 7, 0xFF000000);
				}
			}
		}
		fb->redraw();
		return;
		int fN = 1000;
		V3 C = tmeshes[1].GetCenter();
		V3 adir(1.0f, 1.0f, 0.0f);
		adir = adir.Normalized();
		for (int fi = 0; fi < fN; fi++) {
//			tmeshes[1].Rotate(C, adir, 1.0f);
			tmeshes[1].Translate(V3(0.0f, 0.0f, -0.1f));
			Render();
			Fl::check();
		}
		return;
	}

	{
		RCSetup();
		return;
	}


	{
		ppc->TranslateRightLeft(-10.0f);
		smSetup = 1;
		hwfb->redraw();
		return;
	}

	{
		int fN = 1000;
		for (int i = 0; i < fN; i++) {
			scene->p1 = (float)i / (float) (fN-1);
			scene->morphFraction = scene->p1;
			gpufb->redraw();
			Fl::check();
		}
		return;

	}

	{
		int fN = 1000;
		PPC ppc0(*ppc);
		PPC ppc1(*ppc);
		ppc1.C = ppc1.C + V3(20.0f, 30.0f, -50.0f);
		ppc1.SetPose(ppc1.C, tmeshes[1].GetCenter(), V3(0.0f, 1.0f, 0.0f));
		for (int fi = 0; fi < fN; fi++) {
			ppc->SetInterpolated(&ppc0, &ppc1, fN, fi);
			Render(fb, ppc);
			fb->redraw();
			hwfb->redraw();
			gpufb->redraw();
			Fl::check();
		}
		*ppc = ppc0;
		return;

	}


	{

		RenderRayTracing(fb, ppc);
		fb->redraw();
		Fl::check();
		return;

	}

	{

		int fN = 10000;
		for (int fi = 0; fi < fN; fi++) {
			Render();
			Fl::check();
			tmeshes[5].Translate(V3(0.0f, 0.0f, -1.0f));
			if ((fi % 100) == 0)
				cerr << "INFO: " << (float)((int)((float)fi / (float)fN * 10000.0f)) / 100.0f << "%        \r";
		}
		return;

	}

	{
		fb3->SaveAsTiff("mydbg/ssivis.tif");
		return;

	}

	{

		V3 tcenter = tmeshes[1].GetCenter();
		V3 newC = V3(20.0f, 50.0f, -30.0f);
		ppc->SetPose(newC, tcenter, V3(0.0f, 1.0f, 0.0f));

		int fN = 1000;
		V3 L0 = L;
		V3 L1 = tmeshes[1].GetCenter();
		L = L0 + (L1 - L0)*0.75f;
		for (int fi = 0; fi < fN; fi++) {
//			L = L0 + (L1 - L0)*(float)fi / (float)(fN - 1);
			L = L.RotatePoint(L1, V3(0.0f, 1.0f, 0.0f), 360.0f*3.0f / (float)fN);
			Render();
			Fl::check();
		}
		L = L0;
		return;
	}

	{

		V3 tcenter = tmeshes[1].GetCenter();
		V3 newC = V3(20.0f, 50.0f, -30.0f);
		ppc->SetPose(newC, tcenter, V3(0.0f, 1.0f, 0.0f));
		ppc3->SetPose(V3(0.0f, 50.0f, 100.0f), tcenter, V3(0.0f, 1.0f, 0.0f));
		for (int i = 0; i < 200; i++) {
			Render();
			vf += 1.0f;
			Fl::check();
			return;
		}
		return;

		V3 aDir(0.0f, 1.0f, 0.0f);
		for (int i = 0; i < 100; i++) {
			Render();
			Fl::check();
//			tmeshes[1].Rotate(tcenter, aDir, 1.0f);
			ppc->PanLeftRight(1.0f);
		}
		return;

	}

	{
		fb->SaveAsTiff("mydbg/zb.tif");
		return;
	}

	{

		tmeshes[0].onFlag = 0;
		int fN = 300;
		float tstep = .1f;
		for (int fi = 0; fi < fN; fi++) {
			Render();
			Fl::check();
			ppc->TranslateRightLeft(-tstep);
		}
		return;
	}


	{
		int w = fb->w;
		int h = fb->h;
		float hfov = 90.0f;
		PPC ppc(hfov, w, h);
		V3 cc(0.0f, 0.0f, -100.0f);
		unsigned int color = 0xFF000000;
		float sideLength = 60.0f;
		TMesh tm;
		tm.SetToCube(cc, sideLength, 0xFF0000FF, 0xFF000000);
		int fN = 300;
		float tstep = .1f;
		for (int fi = 0; fi < fN; fi++) {
			fb->SetBGR(0xFFFFFFFF);
//			tm.DrawCubeQuadFaces(fb, &ppc, color);
			tm.DrawWireFrame(fb, &ppc, color);
			fb->redraw();
			Fl::check();
			ppc.TranslateRightLeft(-tstep);
//			ppc.TranslateFrontBack(tstep);
		}
		return;
	}



	{
		int w = fb->w;
		int h = fb->h;
		float hfov = 90.0f;
		PPC ppc(hfov, w, h);
		V3 P(0.0f, 0.0f, -100.0f);

		V3 uP, p;
		ppc.Project(P, p);
		uP = ppc.UnProject(p);
		cerr << uP;

		fb->SetBGR(0xFFF0000);
		V3 tr((float)w, 0.0f, 1.0f);
		V3 trP = ppc.UnProject(tr);
		V3 ptr;
		ppc.Project(trP, ptr);
		fb->DrawSquarePoint(ptr[0], ptr[1], 13, 0xFF00FF00);
		fb->redraw();
		return;

		return;

		V3 Q(0.0f, -10.0f, -50.0f);
		V3 q;

		fb->SetBGR(0xFFFFFFFF);
		for (int i = 0; i < 10; i++) {
			if (!ppc.Project(P, p))
				continue;
			fb->DrawSquarePoint(p[0], p[1], 5, 0xFF000000);

			if (!ppc.Project(Q, q))
				continue;
			fb->DrawSquarePoint(q[0], q[1], 5, 0xFF0000FF);

			fb->redraw();
			Fl::check();
			P = P + V3(10.0f, 0.0f, 0.0f);
			Q = Q + V3(10.0f, 0.0f, 0.0f);
		}


		if (ppc.Project(P, p)) {
			cerr << p << endl;
		}
		else {
			cerr << "INFO: point is behind the head" << endl;
		}

		return;

	}

	{

		M33 m;
		V3 r0(1.0f, 1.0f, 1.0f);
		V3 r1(-2.0f, 2.0f, 2.0f);
		V3 r2(3.0f, -3.0f, 3.0f);
		m[0] = r0;
		m[1] = r1;
		m[2] = r2;
		V3 v(1.0f, 2.0f, 3.0f);
		V3 ret = m*v;
		cerr << ret;
		M33 m1 = m.Inverted();
		cerr << m*m1.GetColumn(0) << m*m1.GetColumn(1) << m*m1.GetColumn(2);
		return;
	}


	{
		M33 m;
		V3 v0(1.0f, 3.0f, -1.0f);
		m[0] = v0;
		cerr << m[0] << endl;
		cerr << m[0][2] << endl;
		m[0][2] = 1000.0f;
		cerr << m[0][2] << endl;
		return;
	}

	{

		V3 v0(2.0f, 2.0f, 2.0f);
		V3 v1(4.0f, 3.0f, 5.0f);
		cerr << v0 + v1;
		cerr << "v0*v1 " << v0*v1 << endl;
		cerr << v0.Length() << endl;
		cerr << (v0.Normalized()).Length() << endl;
		cerr << v0;
		return;

	}

	{
		V3 v;
		v.xyz[0] = 1.0f;
		v.xyz[1] = -1.0f;
		v.xyz[2] = 0.0f;
		cerr << v[0] << endl;
		v[0] = 100.0f;
		cerr << v[0] << endl;
		return;

	}

	fb->LoadTiff("mydbg/im.tif");
	fb->redraw();
	return;
	cerr << "INFO: pressed DBG Button" << endl;

	{
		float uv0[2] = { 10.1f, 20.2f };
		float uv1[2] = { 510.1f, 420.2f };
		unsigned int col = 0xFF000000;
		int fN = 300;
		for (int fi = 0; fi < fN; fi++) {
			fb->SetBGR(0xFFFFFFFF);
//			fb->Draw2DSegment(uv0, uv1, cv, cv);
			uv0[1] += 1.0f;
			uv1[1] -= 1.0f;
			fb->redraw();
			Fl::check();
		}
		fb->SaveAsTiff("mydbg/im.tif");
	}

	return;

	{
		fb->SetBGR(0xFF0000FF);
		fb->SetChecker(0xFF000000, 0xFFFFFFFF, 40);
		fb->SetBGR(0xFFFFFFFF);
		float uv0[2] = { 20.3f, 300.45f };
		float uv1[2] = { 420.73f, 100.45f };
		unsigned int col = 0xFF000000;
//		fb->Draw2DSegment(uv0, uv1, col);
	}

}


void Scene::NewButton() {
	cerr << "INFO: pressed New Button" << endl;
}


void Scene::RenderRayTracing(FrameBuffer *rfb, PPC *rppc) {

	rfb->SetBGR(0xFFFFFFFF);
	rfb->ClearZB();
	for (int v = 0; v < rppc->h; v++) {
		for (int q = 0; q < rppc->w; q++)
			rfb->Set(q, v, 0xFF0000FF);
		rfb->redraw();
		Fl::check();
		for (int q = 0; q < rppc->w; q++)
			rfb->Set(q, v, 0xFFFFFFFF);
		for (int u = 0; u < rppc->w; u++) {
			V3 rdir = rppc->UnProject(V3(.5f + (float)u, .5f + (float)v, 1.0f)) - rppc->C;
			rdir = rdir.Normalized();
			V3 rO = rppc->C;
			V3 rc;
			float currz;
			if (!RayTrace(rO, rdir, 0, rc, currz))
				continue;
			if (rfb->Farther(u, v, currz))
				continue;
			rfb->Set(u, v, rc.GetColor());
		}
		rfb->redraw();
		Fl::check();
	}
	rfb->redraw();

}

int Scene::RayTrace(V3 rO, V3 rdir, int rayOrder, V3& rc, float &currz) {

	currz = 0.0f;
	V3 rrO, rrdir;
	for (int tmi = 0; tmi < tmeshesN; tmi++) {
		if (!tmeshes[tmi].onFlag)
			continue;
		float tmpz;
		V3 tmpc;
		V3 tmprrO, tmprrdir;
		if (tmeshes[tmi].RayTrace(rO, rdir, tmpc, tmpz, tmprrO, tmprrdir)) {
			if (currz < tmpz) {
				currz = tmpz;
				rc = tmpc;
				rrO = tmprrO;
				rrdir = tmprrdir;
			}
		}
	}

	if (currz == 0.0f)
		return 0;

	if (rayOrder == 0 && rrO[0] != FLT_MAX) {
		V3 rrc;
		float rrz;
		if (RayTrace(rrO, rrdir, 1, rrc, rrz)) {
			rc = rc*0.5f + rrc*0.5f;
		}
	}

	return 1;

}

void Scene::Render(RandomCamera *rc, FrameBuffer *rcfb) {


	Render(rcfb, rc->ppc0); // conventional rendering
	rc->ClearRays();
	for (int v = 0; v < rcfb->h; v++) {
		for (int u = 0; u < rcfb->w; u++) {
			V3 r0, r1; // perturbed ray points on near and far planes (d0 and d1)
//			rc->PerturbRay(u, v, rcfb->GetZ(u, v), r0, r1);
			rc->PerturbRay(u, v, r0, r1);
			V3 pColor = V3(1.0f, 1.0f, 1.0f);
			float pz = 0.0f;
			if (RayTrace(r0, (r1 - r0).Normalized(), 0, pColor, pz)) {
				rcfb->Set(u, v, pColor.GetColor());
			}
			rcfb->SetZB(u, v, pz);
		}
	}

}


void Scene::RCSetup() {

	gpufb->hide();
	hwfb->hide();
	for (int tmi = 0; tmi < tmeshesN; tmi++)
		tmeshes[tmi].onFlag = 0;

	float d0 = 4*10.0f;
	float d1 = 4*20.0f;
	float maxD0 = 120.0f;
	float maxD1 = 60.0f;
	float _dnear = 100;
	float _dfar = 300;
	float hfov = 55.0f;
	int w = ppc->w;
	int h = ppc->h;
	RandomCamera rcold(hfov, w, h, d0, d1, maxD0, maxD1, _dnear, _dfar);

	float df = 300.0f;
	float vcr = 1.0f;
	RandomCamera rc(hfov, w, h, vcr, df);

	V3 qcolors[4], qverts[4];
	float z0 = -100.0f;
	qverts[0] = V3(-30.0f, 20.0f, z0);
	qverts[1] = V3(-30.0f, -20.0f, z0);
	qverts[2] = V3(+30.0f, -20.0f, z0);
	qverts[3] = V3(+30.0f, 20.0f, z0);
	qcolors[0] = V3(0.0f, 0.0f, 0.0f);
	qcolors[1] = V3(0.0f, 1.0f, 0.0f);
	qcolors[2] = V3(1.0f, 1.0f, 0.0f);
	qcolors[3] = V3(1.0f, 0.0f, 0.0f);
	tmeshes[7].SetQuad(qverts, qcolors);
	tmeshes[7].texture = new FrameBuffer(10, 10, 60, 40, 0);
	tmeshes[7].texture->SetChecker(0xFF0000FF, 0xFF0000AA, 10);
	tmeshes[7].onFlag = 0;

	float z1 = -300.0f;
	qverts[0] = V3(-300.0f, 200.0f, z1);
	qverts[1] = V3(-300.0f, -200.0f, z1);
	qverts[2] = V3(+300.0f, -200.0f, z1);
	qverts[3] = V3(+300.0f, 200.0f, z1);
	tmeshes[8].SetQuad(qverts, qcolors);
	tmeshes[8].texture = new FrameBuffer(10, 10, 60, 40, 0);
	tmeshes[8].texture->SetChecker(0xFF00FF00, 0xFF00AA00, 10);
	tmeshes[8].onFlag = 0;

	V3 greenv(0.0f, 1.0f, 0.0f);
	V3 redv(1.0f, 0.0f, 0.0f);
	tmeshes[7].SetTesselatedRectangle(60.0f, 40.0f, 60, 40, redv);
	tmeshes[7].Translate(V3(0.0f, 0.0f, -100.0f));
	tmeshes[7].texture = 0;
	tmeshes[7].onFlag = 1;

	tmeshes[8].SetTesselatedRectangle(600.0f, 400.0f, 600, 400, greenv);
	tmeshes[8].Translate(V3(0.0f, 0.0f, -300.0f));
	tmeshes[8].texture = 0;
	tmeshes[8].onFlag = 1;
//	Render();

	AABB aabb;
	tmeshes[11].SetAABB(aabb);
	float diagdim = (aabb.corners[1] - aabb.corners[0]).Length();
	float newdim = 100.9f;
	tmeshes[11].Scale(newdim / diagdim);
	tmeshes[11].Rotate(tmeshes[11].GetCenter(), V3(0.0f, 1.0f, 0.0f), 90.0f);
	tmeshes[11].Rotate(tmeshes[11].GetCenter(), V3(0.0f, 0.0f, -1.0f), -90.0f);

	tmeshes[7].onFlag = tmeshes[8].onFlag = 0;
	tmeshes[11].onFlag = 1;

	fb->SetBGR(0xFFFFFFFF);
	fb->ClearZB();

	for (int tmi = 0; tmi < tmeshesN; tmi++) {
		if (!tmeshes[tmi].onFlag)
			continue;
		tmeshes[tmi].RenderRC(fb, ppc, vcr);
	}

	fb->redraw();

//	return;

//	Render(&rc, fb);
	fb->label("Random Camera Image");
	fb->redraw();

	FrameBuffer *pbfb = new FrameBuffer(30 + fb->w + 30, 100, fb->w/4, fb->h/4, 0);
	pbfb->label("Point Based Reconstruction");
	pbfb->show();
	PPC pbppc(*rc.ppc0);
	pbppc.w /= 4;
	pbppc.h /= 4;
	pbppc.a = pbppc.a * 4.0f;
	pbppc.b = pbppc.b * 4.0f;

	int fN = 1000;
	float dt = vcr * 1.0f / (float)fN;
	int pSize = 3;
	for (int fi = 0; fi < fN; fi++) {
//		rc.PointBasedRender(fb, pSize, pbfb, &pbppc);
		pbfb->SetBGR(0xFFFFFFFF);
		pbfb->ClearZB();
		for (int v = 0; v < fb->h; v++) {
			for (int u = 0; u < fb->w; u++) {
				int uv = (fb->h - v - 1)*fb->w + u;
				pbfb->Draw3DPoint(fb->xyz[uv], &pbppc, fb->pix[uv], pSize);
			}
		}
		pbfb->redraw();
		Fl::check();
		pbppc.TranslateRightLeft(dt);
	}

	pbfb->SaveAsTiff("mydbg/RightMostPointBasedReconstruction.tif");
	fb->SaveAsTiff("mydbg/RandomCamera.tif");

}

void Scene::PerSessionHWInit() {

	if (hwPerSessionInit)
		return;

	glewInit();
	GLuint texs[3];
	glGenTextures(3, texs);

	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, texs[0]);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	FrameBuffer *tex = new FrameBuffer(0, 0, 64, 64, 0);
	tex->SetChecker(0xFF000000, 0xFFFFFFFF, 8);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, 
		tex->w, tex->h, 0, GL_RGBA, GL_UNSIGNED_BYTE, tex->pix);

	glActiveTexture(GL_TEXTURE1);
	glBindTexture(GL_TEXTURE_2D, texs[1]);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	FrameBuffer *tex2 = new FrameBuffer(0, 0, 64, 64, 0);
	tex2->SetChecker(0xFF0000FF, 0xFFFFFFFF, 8);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA,
		tex2->w, tex2->h, 0, GL_RGBA, GL_UNSIGNED_BYTE, tex2->pix);


	glActiveTexture(GL_TEXTURE2);
	glBindTexture(GL_TEXTURE_2D, texs[2]);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

	FrameBuffer *tex3 = new FrameBuffer(0, 0, 64, 64, 0);
	tex3->LoadTiff("mydbg/sm.tif");
//	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA,
//		tex3->w, tex3->h, 0, GL_RGBA, GL_UNSIGNED_BYTE, tex3->pix);

	ifstream zb("mydbg/smzb");
	zb.read((char*)tex3->zb, sizeof(float)*tex3->w*tex3->h);
	zb.close();
	glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT,
		tex3->w, tex3->h, 0, GL_DEPTH_COMPONENT, GL_FLOAT, tex3->zb);

	hwPerSessionInit = 1;

}

void Scene::SetupSM() {

	smSetup = 0;
	glReadPixels(0, 0, hwfb->w, hwfb->h, GL_RGBA, GL_UNSIGNED_BYTE, hwfb->pix);
	hwfb->SaveAsTiff("mydbg/sm.tif");
	glReadPixels(0, 0, hwfb->w, hwfb->h, GL_DEPTH_COMPONENT, GL_FLOAT, hwfb->zb);

	float minz = 1.0f;
	for (int v = 0; v < hwfb->h; v++)
		for (int u = 0; u < hwfb->w; u++) {
			float z = hwfb->GetZ(u, v);
			if (z < minz)
				minz = z;
		}
	for (int v = 0; v < hwfb->h; v++)
		for (int u = 0; u < hwfb->w; u++) {
			float z = hwfb->GetZ(u, v);
			z = (z - minz) / (1.0f - minz);
			V3 zv(z, z, z);
			hwfb->Set(u, v, zv.GetColor());
		}
	hwfb->SaveAsTiff("mydbg/smz.tif");
	ofstream zb("mydbg/smzb");
	zb.write((char*)hwfb->zb, sizeof(float)*hwfb->w*hwfb->h);
	zb.close();

}


void Scene::TrackImage(FrameBuffer *itt, V3 *corners, PPC *ppc0, FrameBuffer *frame, PPC *ppc) {



}


void Scene::TrackImageSetup() {


	eri->hide();
	fb->hide();

	imageToTrack = new FrameBuffer(20, 20, 10, 10, 0);
	imageToTrack->LoadTiff("mydbg/photos/textSmall.tif");
	cerr << "Image to track res: " << imageToTrack->w << " " << imageToTrack->h << endl;
	imageToTrack->label("image to track");
	imageToTrack->show();

	currFrame = new FrameBuffer(20 + imageToTrack->w + 20, 20, 10, 10, 0);
	currFrame->LoadTiff("mydbg/photos/frame0.tif");
	currFrame->label("frame");
	
//	currFrame->SetFromImage(imageToTrack);

	visfb = new FrameBuffer(20 + imageToTrack->w + 20, 20, 10, 10, 0);
	visfb->LoadTiff("mydbg/photos/frame0.tif");
	visfb->SetFromImage(currFrame);
	visfb->label("tracking vis");
	visfb->show();
	gui->uiw->position(20 + imageToTrack->w + 20 + visfb->w + 20, 20);

	V3 bbl0(225, 222, 0);
	V3 bbl1(221, 619, 0);
	float pix2cm = (bbl0 - bbl1).Length() / 9.0f / 5.0f;

//	V3 ittCorners[4]; // bottom left, bottom right, top right, top left
	ittCorners[0] = V3(0.0f, 0.0f, 0.0f);
	ittCorners[1] = V3((float)imageToTrack->w / pix2cm, 0.0f, 0.0f);
	ittCorners[2] = V3((float)imageToTrack->w / pix2cm, 0.0f, -(float)imageToTrack->h / pix2cm);
	ittCorners[3] = V3(0.0f, 0.0f, -(float)imageToTrack->h / pix2cm);

	float hfov = 51.0f;
	currPPC = new PPC(hfov, currFrame->w, currFrame->h);
	V3 ic = (ittCorners[3] + ittCorners[1]) / 2.0f;
	float h = (ittCorners[1] - ittCorners[0]).Length() / 2.0f / tan(hfov/2.0f/180.0f*3.14159265f);
	h += 10.0f;
	V3 eyeGuess(0.0f, h, 0.0f); eyeGuess = eyeGuess + ic;
	eyeGuess = ic + V3(0.0f, 30.0f, 30.0f);
	currPPC->SetPose(eyeGuess, ic, V3(0.0f, 1.0f, 0.0f));

//	currPPC->Load("mydbg/trackedPPC.txt");

	cerr << "INFO: initial guess color difference: " 
		<< ColorDifference(imageToTrack, ittCorners, currFrame, currPPC) << endl;
	currPPC->Save("mydbg/trackingGuess.txt");
	visfb->SetFromImage(currFrame);
	DrawImageRectangle(visfb, currPPC, ittCorners, V3(0.0f, 1.0f, 0.0f));
	visfb->redraw();

}


float Scene::ColorDifference(FrameBuffer *itt, V3 *corners, FrameBuffer *frame, 
	PPC *ppc) {

	V3 TL = corners[3];
	V3 du = (corners[1] - corners[0]) / (float) itt->w;
	V3 dv = (corners[0] - corners[3]) / (float) itt->h;

	float ret = 0.0f;
	int pixN = 0;
	for (int v = 0; v < itt->h; v++) {
		for (int u = 0; u < itt->w; u++) {
			V3 P = TL + du * (.5f + (float)u) + dv * (.5f + (float)v);
			V3 PP;
			if (!ppc->Project(P, PP))
				continue;
			V3 c0; c0.SetFromColor(itt->Get(u, v));
			int u1 = (int)PP[0];
			int v1 = (int)PP[1];
			if (u1 < 0 || u1 > frame->w - 1 || v1 < 0 || v1 > frame->h - 1)
				continue;
			V3 c1; c1.SetFromColor(frame->Get(u1, v1));
			V3 c10 = c1 - c0;
			ret += c10*c10;
			pixN++;
		}
	}
	if (pixN != itt->h * itt->w)
		return FLT_MAX;
	ret /= (float)pixN;
	ret /= 3;
	ret = sqrtf(ret);
	return ret;

}

void Scene::DrawImageRectangle(FrameBuffer *frame, PPC *ppc, V3 *corners, V3 color) {

	for (int i = 0; i < 4; i++)
		frame->Draw3DSegment(corners[i], corners[(i + 1) % 4], ppc, color, color);

}


void Scene::ImageTracking() {

	PPC *trackedPPC = new PPC(*currPPC);
	float minError = ColorDifference(imageToTrack, ittCorners, currFrame, trackedPPC);
	V3 panTiltRoll; // rotational dofs
	V3 trans; // translational degrees of freedom (dofs)
	V3 panTiltRollRange(20.0f, 20.0f, 5.0f);
	V3 transRange(20.0f, 20.0f, 20.0f);
	int panTiltRollStepsN[3] = { 7, 7, 3 };
	int transStepsN[3] = { 7, 7, 7 };
	int stepsN = panTiltRollStepsN[0] * panTiltRollStepsN[1] * panTiltRollStepsN[2] *
		transStepsN[0] * transStepsN[1] * transStepsN[2];
	int stepi = 0;
	PPC *ppc0 = new PPC(*currPPC);

	for (int txi = 0; txi < transStepsN[0]; txi++) {
		float txs = transRange[0] / (float) (transStepsN[0] - 1);
		float tx = -transRange[0] / 2.0f + txs*(float)txi;
		for (int tyi = 0; tyi < transStepsN[1]; tyi++) {
			float tys = transRange[1] / (float)(transStepsN[1] - 1);
			float ty = -transRange[1] / 2.0f + tys*(float)tyi;
			for (int tzi = 0; tzi < transStepsN[2]; tzi++) {
				float tzs = transRange[2] / (float)(transStepsN[2] - 1);
				float tz = -transRange[2] / 2.0f + tzs*(float)tzi;
				for (int rxi = 0; rxi < panTiltRollStepsN[0]; rxi++) {
					float pans = panTiltRollRange[0] / (float)(panTiltRollStepsN[0] - 1);
					float pan = -panTiltRollRange[0] / 2.0f + pans*(float)rxi;
					for (int ryi = 0; ryi < panTiltRollStepsN[1]; ryi++) {
						float tilts = panTiltRollRange[1] / (float)(panTiltRollStepsN[1] - 1);
						float tilt = -panTiltRollRange[1] / 2.0f + tilts*(float)ryi;
						for (int rzi = 0; rzi < panTiltRollStepsN[2]; rzi++) {
							float rolls = panTiltRollRange[2] / (float)(panTiltRollStepsN[2] - 1);
							float roll = -panTiltRollRange[2] / 2.0f + rolls*(float)rzi;
							*currPPC = *ppc0;
							currPPC->PanLeftRight(pan);
							currPPC->TiltUpDown(tilt);
							currPPC->RollLeftRight(roll);
							currPPC->TranslateFrontBack(tz);
							currPPC->TranslateUpDown(ty);
							currPPC->TranslateRightLeft(tx);
							float currColorDiff = ColorDifference(imageToTrack, ittCorners, currFrame, currPPC);
							if (currColorDiff < minError) {
								cerr << endl << "INFO: smaller error found: " << currColorDiff << endl;
								*trackedPPC = *currPPC;
								minError = currColorDiff;
								visfb->SetFromImage(currFrame);
								visfb->ClearZB();
								DrawImageRectangle(visfb, currPPC, ittCorners, V3(0.0f, 1.0f, 0.0f));
								visfb->redraw();
								Fl::check();
							}
							stepi++;
							if ((stepi % 101) == 0) {
								float prgr = (float)(int)((float)stepi / (float)stepsN * 10000.0f);
								prgr /= 100.0f;
								cerr << "INFO: " << prgr << "%    \r";    
#if 1
								visfb->SetFromImage(currFrame);
								visfb->ClearZB();
								DrawImageRectangle(visfb, currPPC, ittCorners, V3(1.0f, 0.0f, 0.0f));
								DrawImageRectangle(visfb, trackedPPC, ittCorners, V3(0.0f, 1.0f, 0.0f));
								visfb->redraw();
								Fl::check();
#endif
							}
						}
					}
				}
			}
		}
	}
	trackedPPC->Save("mydbg/trackedPPC.txt");
	cerr << "INFO: done" << endl;

}

void Scene::TransparentDisplay() {

	piPhi = piDelta = 0.0f;

	fb->hide();

	float displayW = 140.f; // mm
	float displayH = 230.f; // mm
	float tabletCamVFOV = 80.0f;
	float tabletCamF = displayH / 2.0f / tanf(tabletCamVFOV / 2.0f / 180.0f * 3.14159265f);
	float tabletCamHFOV = 2.0f * atanf(displayW / 2.0f / tabletCamF) * 180.0f / 3.14159265f;
	cerr << "INFO: tablet cam hfov = " << tabletCamHFOV << endl;
	PPC *tabletCam = new PPC(tabletCamHFOV, (int)displayW, (int)displayH);
	tabletCam->C = V3(displayW / 2.0f, displayH / 2.0f, 0.0f);

	float userCamF = 390.f; // mm
	float userCamHFOV = 2.0f * atanf(displayW / 2.0f / userCamF)*180.0f / 3.14159265f;
	cerr << "INFO: user cam hfov: " << userCamHFOV << endl;
	PPC *userCam = new PPC(userCamHFOV, (int)displayW, (int)displayH);
	userCam->C = V3(0.0f, 0.0f, userCamF);

	V3 p((float)userCam->w * 3.0f / 4.0f, (float)userCam->h/2.0f, 1);
	float dn = 1000.0f;
	float df = 20000.0f;
	int stepsN = 191;
	FrameBuffer *conv, *transp, *truth;
	int visw = 400, vish = 600;
	conv = new FrameBuffer(20, 20, visw, vish, 0);
	conv->label("Conventional");
	conv->show();
	transp = new FrameBuffer(20+conv->w+20, 20, visw, vish, 0);
	transp->label("Transparent");
	transp->show();
	truth = new FrameBuffer(20 + conv->w + 20 + transp->w + 20, 20, visw, vish, 0);
	truth->label("truth");
	truth->show();
	float vishfov = 2.0f * 180.0f / 3.14159265f * atan((float)visw / (float)vish);
	cerr << "vishfov " << vishfov << endl;
	PPC *visppc = new PPC(50.0f, conv->w, conv->h);
	visppc->C = userCam->C;
	int maxGotBetter = 0;
	int avgGotBetter = 0;
	FrameBuffer *tex = new FrameBuffer(10, 10, 1, 1, 1);
	tex->LoadTiff("mydbg/munchthedance.tif");
//	tex->SaveAsTiff("mydbg/munchthedance.tif");
	//tex->show();
	FrameBuffer *dispr, *dispg, *dispb;
	dispr = new FrameBuffer(20, 450, userCam->w, userCam->h, 0);
	dispr->label("conv frame");
	dispr->show();
	dispg = new FrameBuffer(20 + conv->w + 20, 450, userCam->w, userCam->h, 0);
	dispg->label("transp frame");
	dispg->show();
	dispb = new FrameBuffer(20 + conv->w + 20 + conv->w + 20, 450, userCam->w, userCam->h, 0);
	dispb->label("truth frame");
	dispb->show();
	char folder[1000];
	strcpy_s(folder, "C:/Users/Voicu Popescu/Dropbox/Research/TransparentDisplay/transparentDisplaySimulator");
	char fname[1000];
	sprintf_s(fname, "%s/data.txt", folder);
	ofstream ofs(fname);
	ofs << "Distance RepErrAvgTransp  RepErrAvgConv ";
	ofs << "RepErrMaxTransp RepErrMaxConv ";
	ofs << "MissTransp MissConv ";
	ofs << "RedunTransp RedunConv ";
	ofs << endl;
	for (int si = 0; si < stepsN; si++) {
		float d = dn + (df - dn)*(float)si / (float)(stepsN - 1);
		float ph = 2.0f*(d + userCamF);
		VisTransp(transp, tex, d, ph, userCam, tabletCam, visppc, 0);
		VisTransp(conv, tex, d, ph, userCam, tabletCam, visppc, 1);
		VisTransp(truth, tex, d, ph, userCam, tabletCam, visppc, 2);
		VisDisp(dispg, tex, d, ph, userCam, tabletCam, 0);
		VisDisp(dispr, tex, d, ph, userCam, tabletCam, 1);
		VisDisp(dispb, tex, d, ph, userCam, tabletCam, 2);
		float gerrAvg = 0.0f;
		float rerrAvg = 0.0f;
		float gerrMax = 0.0f;
		float rerrMax = 0.0f;
		float gerrMiss = 0.0f;
		float rerrMiss = 0.0f;
		float gerrRedun = 0.0f;
		float rerrRedun = 0.0f;
		float fu = userCam->GetF();
		float ft = tabletCam->GetF();
		for (int v = 0; v < userCam->h; v++) {
			for (int u = 0; u < userCam->w; u++) {
				// compute transparency error on transparent display
				// uv is on transparent display
				V3 g;
				g[0] = (float)u + .5f;
				g[1] = (float)v + .5f;
				V3 rayg = userCam->GetImagePlanePoint(g[0], g[1]) - userCam->C;
				V3 gt;
				tabletCam->Project(rayg + tabletCam->C, gt);
				gt[2] = ft / d;
				V3 gP = tabletCam->UnProject(gt);
				V3 gb; userCam->Project(gP, gb);
				if (gb[0] < 0.0f || gb[0] > (float) userCam->w
					|| gb[1] < 0.0f || gb[1] > (float) userCam->h) {
					gerrRedun += 1.0f;
					dispg->Set(u, v, 0xFF00FFFF);
				}
				float gerr = (V3(g[0], g[1], 0.0f) - V3(gb[0], gb[1], 0.0f)).Length();
				V3 errc(0.0f, 0.0f, 0.0f);
				errc[0] = gerr / 20.0f;
//				dispg->Set(u, v, errc.GetColor());
				gerrAvg += gerr;
				if (gerrMax < gerr)
					gerrMax = gerr;

				// compute transparency error on convnetional display
				// uv is on convnetional display
				V3 r;
				r[0] = (float)u + .5f;
				r[1] = (float)v + .5f;
				r[2] = ft / d;
				V3 rP = tabletCam->UnProject(r);
				V3 rb; userCam->Project(rP, rb);
				if (rb[0] < 0.0f || rb[0] > (float) userCam->w
					|| rb[1] < 0.0f || rb[1] > (float) userCam->h) {
					rerrRedun += 1.0f;
					if (dispr->Get(u, v) != 0xFF0000FF 
						&& dispr->Get(u, v) != 0xFF00FF00
						&& dispr->Get(u, v) != 0xFF000000)
						dispr->Set(u, v, 0xFF00FFFF);
				}
				float rerr = (V3(r[0], r[1], 0.0f) - V3(rb[0], rb[1], 0.0f)).Length();
				rerrAvg += rerr;
				if (rerrMax < rerr)
					rerrMax = rerr;

				if (u == 5 * userCam->w / 8 && v == userCam->h / 2 && 0) {
					dispg->DrawDiagonalCrossPadded(g[0], g[1], 17, 0xFF00FF00);
					dispb->DrawStraightCrossPadded(gb[0], gb[1], 17, 0xFF00FF00);
					dispg->DrawStraightCrossPadded(gb[0], gb[1], 17, 0xFF00FF00);
					dispr->DrawDiagonalCrossPadded(r[0], r[1], 17, 0xFF0000FF);
					dispb->DrawStraightCrossPadded(rb[0], rb[1], 17, 0xFF0000FF);
					dispr->DrawStraightCrossPadded(rb[0], rb[1], 17, 0xFF0000FF);
				}

				V3 M = userCam->UnProject(V3(.5f + (float) u, 
					.5f + (float) v, fu / (d + fu)));
				V3 m;  tabletCam->Project(M, m);
				if (m[0] < 0.0f || m[0] > (float) tabletCam->w ||
					m[1] < 0.0f || m[1] > (float) tabletCam->h) {
					rerrMiss += 1.0f;
					gerrMiss += 1.0f;
					if (dispb->Get(u, v) != 0xFF0000FF && dispb->Get(u, v) != 0xFF00FF00)
						dispb->Set(u, v, 0xFFFFFF00);
				}
				else {
					userCam->Project(userCam->C + 
						(tabletCam->GetImagePlanePoint(m[0], m[1])-
							tabletCam->C), g);
					if (g[0] < 0.0f || g[0] > (float) userCam->w ||
						g[1] < 0.0f || g[1] > (float) userCam->h) {
						gerrMiss += 1.0f;
						if (dispb->Get(u, v) != 0xFF0000FF 
							&& dispb->Get(u, v) != 0xFF00FF00
							&& dispb->Get(u, v) != 0xFF000000)
							dispb->Set(u, v, 0xFFFF0000);
					}
				}
			}
		}
		dispb->redraw();

		gerrAvg /= (float)(userCam->w * userCam->h);
		rerrAvg /= (float)(userCam->w * userCam->h);
		cerr << "d = " << d << " avg: " << gerrAvg << " " << rerrAvg << " Max: " << gerrMax << " " << rerrMax << endl;
		ofs << d << " " << gerrAvg << " " << rerrAvg << " ";
		ofs << gerrMax << " " << rerrMax << " ";
		if (gerrAvg < rerrAvg && !avgGotBetter) {
			cerr << "INFO: average got better" << endl;
			avgGotBetter = 1;
		}
		if (gerrMax < rerrMax && !maxGotBetter) {
			cerr << "INFO: max got better" << endl;
			maxGotBetter = 1;
		}

		gerrMiss /= (float)(userCam->w * userCam->h);
		rerrMiss /= (float)(userCam->w * userCam->h);
		cerr << "d = " << d << " missing: " << gerrMiss * 100.0f << "% "
			<< rerrMiss * 100.0f << "%" << endl;
		ofs << gerrMiss << " " << rerrMiss << " ";

		gerrRedun /= (float)(userCam->w * userCam->h);
		rerrRedun /= (float)(userCam->w * userCam->h);
		cerr << "d = " << d << " Redundant: " << gerrRedun * 100.0f << "% "
			<< rerrRedun * 100.0f << "%" << endl;
		ofs << gerrRedun << " " << rerrRedun << endl;

		Fl::check();
		if ((si % 10) == 0) {
			sprintf_s(fname, "%s/VisConventional_d%d.tif", folder, (int)d);
			conv->SaveAsTiff(fname);
			sprintf_s(fname, "%s/VisTransparent_d%d.tif", folder, (int)d);
			transp->SaveAsTiff(fname);
			sprintf_s(fname, "%s/VisTruth_d%d.tif", folder, (int)d);
			truth->SaveAsTiff(fname);

			sprintf_s(fname, "%s/DisplayConventional_d%d.tif", folder, (int)d);
			dispr->SaveAsTiff(fname);
			sprintf_s(fname, "%s/DisplayTransparent_d%d.tif", folder, (int)d);
			dispg->SaveAsTiff(fname);
			sprintf_s(fname, "%s/DisplayTruth_d%d.tif", folder, (int)d);
			dispb->SaveAsTiff(fname);

		}
	}
	ofs.close();

	return;

}

void Scene::VisTransp(FrameBuffer *visfb, FrameBuffer *tex, float d, float ph,
	PPC *userCamera, PPC* tabletCamera, PPC *visppc, int cond) {

	visfb->SetBGR(0xFFFFFFFF);
	for (int v = 0; v < visppc->h; v++) {
		for (int u = 0; u < visppc->w; u++) {
			V3 ray = visppc->GetImagePlanePoint(.5f + (float)u,
				.5f + (float)v) - visppc->C;
			V3 rayO = visppc->C;
			V3 ip;
			// intersect with display
			if (IntersetRayWithPlane(visppc->C, ray, V3(0.0f, 0.0f, 0.0f), V3(0.0f, 0.0f, 1.0f), ip)) {
				V3 dpp;
				if (userCamera->Project(ip, dpp)) {
					if (!(dpp[0] < 0.0f || dpp[0] > (float) userCamera->w ||
						dpp[1] < 0.0f || dpp[1] > (float) userCamera->h)) {
						float framew = 3.0f;
						if (dpp[0] < framew || dpp[0] > (float) userCamera->w - framew ||
							dpp[1] < framew || dpp[1] > (float)userCamera->h - framew) {
							visfb->Set(u, v, 0xFF888888);
							continue;
						}
						if (cond == 0) {
							// transp distant geometry
							rayO = tabletCamera->C;
						}
						else if (cond == 1) {
							// conv
							ray = tabletCamera->GetImagePlanePoint(dpp[0], dpp[1]) -
								tabletCamera->C;
						}
						else if (cond == 2) {
							// no change;
							// truth
						}
						else if (cond == 3) {
							// transp planar proxy
							V3 proxyPoint(0.0f, 0.0f, -d);
							V3 proxyNormal(0.0f, 0.0f, 1.0f);
							proxyNormal = proxyNormal.RotateVector(V3(0.0f, 1.0f, 0.0f), piPhi);
							V3 proxyInt;
							if (!IntersetRayWithPlane(rayO, ray, proxyPoint, proxyNormal, proxyInt)) {
								cerr << "INTERNAL ERROR: no intersection with proxy plane" << endl;
								exit(0);
							}
							rayO = tabletCamera->C;
							ray = (proxyInt - rayO).Normalized();
						}
					}
				}
			}
			unsigned int col;
			if (IntWall(rayO, ray, d, ph, tex, col))
				visfb->Set(u, v, col);
		}
	}
	visfb->redraw();
}

int Scene::IntWall(V3 rayO, V3 ray, float d, float ph, FrameBuffer *tex,
	unsigned int &col) {

	V3 wallc = V3(0.0f, 0.0f, -d);
	V3 walln = V3(0.0f, 0.0f, 1.0f);
	walln = (walln.RotateVector(V3(0.0f, 1.0f, 0.0f), piPhi)).Normalized();
	wallc = wallc + walln * piDelta;
	float pw = ph * (float)tex->w / (float)tex->h;

	// intersect with wall
	V3 ip;
	if (!IntersetRayWithPlane(rayO, ray, wallc, walln, ip))
		return 0;
	float x = (ip[0] + pw / 2.0f) / pw;
	float y = (-ip[1] + ph / 2.0f) / ph;
	int pu = (int)(x * (float)tex->w);
	int pv = (int)(y * (float)tex->h);
	if (pu < 0 || pu > tex->w - 1 || pv < 0 || pv > tex->h - 1)
		return 0;
	col = tex->Get(pu, pv);
	return 1;
}


void Scene::VisDisp(FrameBuffer *disp, FrameBuffer *tex, float d, float ph,
	PPC *userCam, PPC *tabletCam, int cond) {

	for (int v = 0; v < disp->h; v++) {
		for (int u = 0; u < disp->w; u++) {
			float uf = .5f + (float)u;
			float vf = .5f + (float)v;
			disp->Set(u, v, 0xFFFFFFFF);
			V3 rayO, ray;
			if (cond == 0) {
				rayO = tabletCam->C;
				ray = userCam->GetImagePlanePoint(uf, vf) - userCam->C;
			}
			else if (cond == 1) {
				rayO = tabletCam->C;
				ray = tabletCam->GetImagePlanePoint(uf, vf) - rayO;
			}
			else if (cond == 2) {
				rayO = userCam->C;
				ray = userCam->GetImagePlanePoint(uf, vf) - rayO;
			}
			else if (cond == 3) {
				// transp planar proxy
				rayO = userCam->C;
				ray = userCam->GetImagePlanePoint(uf, vf) - rayO;
				V3 proxyPoint(0.0f, 0.0f, -d);
				V3 proxyNormal(0.0f, 0.0f, 1.0f);
				proxyNormal = proxyNormal.RotateVector(V3(0.0f, 1.0f, 0.0f), piPhi);
				V3 proxyInt;
				if (!IntersetRayWithPlane(rayO, ray, proxyPoint, proxyNormal, proxyInt)) {
					cerr << "INTERNAL ERROR: no intersection with proxy plane" << endl;
					exit(0);
				}
				rayO = tabletCam->C;
				ray = (proxyInt - rayO).Normalized();
			}
			unsigned int col;
			if (IntWall(rayO, ray, d, ph, tex, col))
				disp->Set(u, v, col);
		}
	}
	disp->redraw();

}


void Scene::TransparentDisplayPlanarProxy() {

	fb->hide();

	float displayW = 140.f; // mm
	float displayH = 230.f; // mm
	float tabletCamVFOV = 80.0f;
	float tabletCamF = displayH / 2.0f / tanf(tabletCamVFOV / 2.0f / 180.0f * 3.14159265f);
	float tabletCamHFOV = 2.0f * atanf(displayW / 2.0f / tabletCamF) * 180.0f / 3.14159265f;
	cerr << "INFO: tablet cam hfov = " << tabletCamHFOV << endl;
	PPC *tabletCam = new PPC(tabletCamHFOV, (int)displayW, (int)displayH);
	tabletCam->C = V3(displayW / 2.0f, displayH / 2.0f, 0.0f);

	float userCamF = 390.f; // mm
	float userCamHFOV = 2.0f * atanf(displayW / 2.0f / userCamF)*180.0f / 3.14159265f;
	cerr << "INFO: user cam hfov: " << userCamHFOV << endl;
	PPC *userCam = new PPC(userCamHFOV, (int)displayW, (int)displayH);
	userCam->C = V3(0.0f, 0.0f, userCamF);

	V3 p((float)userCam->w * 3.0f / 4.0f, (float)userCam->h / 2.0f, 1);
	FrameBuffer *conv, *transp, *truth;
	int visw = 400, vish = 600;
	conv = new FrameBuffer(20, 20, visw, vish, 0);
	conv->label("Conventional");
	conv->show();
	transp = new FrameBuffer(20 + conv->w + 20, 20, visw, vish, 0);
	transp->label("Transparent");
	transp->show();
	truth = new FrameBuffer(20 + conv->w + 20 + transp->w + 20, 20, visw, vish, 0);
	truth->label("truth");
	truth->show();
	float vishfov = 2.0f * 180.0f / 3.14159265f * atan((float)visw / (float)vish);
	cerr << "vishfov " << vishfov << endl;
	PPC *visppc = new PPC(50.0f, conv->w, conv->h);
	visppc->C = userCam->C;
	int maxGotBetter = 0;
	int avgGotBetter = 0;
	FrameBuffer *tex = new FrameBuffer(10, 10, 1, 1, 1);
	tex->LoadTiff("mydbg/munchthedance.tif");
	//	tex->SaveAsTiff("mydbg/munchthedance.tif");
	//tex->show();
	FrameBuffer *dispr, *dispg, *dispb;
	dispr = new FrameBuffer(20, 450, userCam->w, userCam->h, 0);
	dispr->label("conv frame");
	dispr->show();
	dispg = new FrameBuffer(20 + conv->w + 20, 450, userCam->w, userCam->h, 0);
	dispg->label("transp frame");
	dispg->show();
	dispb = new FrameBuffer(20 + conv->w + 20 + conv->w + 20, 450, userCam->w, userCam->h, 0);
	dispb->label("truth frame");
	dispb->show();
	char folder[1000];
	strcpy_s(folder, "C:/Users/Voicu Popescu/Dropbox/Research/TransparentDisplay/transparentDisplaySimulator");
	char fname[1000];
	sprintf_s(fname, "%s/data.txt", folder);
	ofstream ofs(fname);
	ofs << "Distance RepErrAvgTranspProxy  RepErrAvgConv ";
	ofs << "RepErrMaxTranspProxy RepErrMaxConv ";
	ofs << "MissTranspProxy MissConv ";
	ofs << "RedunTranspProxy RedunConv ";
	ofs << endl;
	piPhi = 45.0f;
	float piDeltaMax = 500.0f;
	float d = 1000.0f;
	int stepsN = 51;
	for (int si = 0; si < stepsN; si++) {
		piDelta = -piDeltaMax / 2.0f + piDeltaMax * (float)si / (float)(stepsN - 1);
		float ph = 2.0f*(d + userCamF);
		VisTransp(transp, tex, d, ph, userCam, tabletCam, visppc, 3);
		VisTransp(conv, tex, d, ph, userCam, tabletCam, visppc, 1);
		VisTransp(truth, tex, d, ph, userCam, tabletCam, visppc, 2);
		VisDisp(dispg, tex, d, ph, userCam, tabletCam, 3);
		VisDisp(dispr, tex, d, ph, userCam, tabletCam, 1);
		VisDisp(dispb, tex, d, ph, userCam, tabletCam, 2);
		float gerrAvg = 0.0f;
		float rerrAvg = 0.0f;
		float gerrMax = 0.0f;
		float rerrMax = 0.0f;
		float gerrMiss = 0.0f;
		float rerrMiss = 0.0f;
		float gerrRedun = 0.0f;
		float rerrRedun = 0.0f;
		float fu = userCam->GetF();
		float ft = tabletCam->GetF();
		for (int v = 0; v < userCam->h; v++) {
			for (int u = 0; u < userCam->w; u++) {
				// compute transparency error on transparent display, planar proxy
				// uv is on transparent display
				V3 g;
				g[0] = (float)u + .5f;
				g[1] = (float)v + .5f;
				V3 rayg = userCam->GetImagePlanePoint(g[0], g[1]) - userCam->C;

				V3 proxyPoint(0.0f, 0.0f, -d);
				V3 proxyNormal(0.0f, 0.0f, 1.0f);
				proxyNormal = proxyNormal.RotateVector(V3(0.0f, 1.0f, 0.0f), piPhi);
				V3 proxyInt;
				if (!IntersetRayWithPlane(userCam->C, rayg, proxyPoint, proxyNormal, proxyInt)) {
					cerr << "INTERNAL ERROR: no intersection with proxy plane" << endl;
					exit(0);
				}
				V3 sceneInt;
				if (!IntersetRayWithPlane(tabletCam->C, proxyInt - tabletCam->C,
					proxyPoint + proxyNormal*piDelta, proxyNormal, sceneInt)) {
					cerr << "INTERNAL ERROR: no intersection with scene plane transp" << endl;
					exit(0);
				}

				V3 gb; userCam->Project(sceneInt, gb);
				if (gb[0] < 0.0f || gb[0] > (float) userCam->w
					|| gb[1] < 0.0f || gb[1] > (float) userCam->h) {
					gerrRedun += 1.0f;
					dispg->Set(u, v, 0xFF00FFFF);
				}
				float gerr = (V3(g[0], g[1], 0.0f) - V3(gb[0], gb[1], 0.0f)).Length();
				gerrAvg += gerr;
				if (gerrMax < gerr)
					gerrMax = gerr;

				// compute transparency error on convnetional display
				// uv is on convnetional display
				V3 r;
				r[0] = (float)u + .5f;
				r[1] = (float)v + .5f;

				V3 sceneIntConv;
				if (!IntersetRayWithPlane(tabletCam->C, 
					tabletCam->GetImagePlanePoint(r[0], r[1]) - tabletCam->C,
					proxyPoint + proxyNormal*piDelta, proxyNormal, sceneIntConv)) {
					cerr << "INTERNAL ERROR: no intersection with scene plane conv" << endl;
					exit(0);
				}

				V3 rb; userCam->Project(sceneIntConv, rb);
				if (rb[0] < 0.0f || rb[0] > (float) userCam->w
					|| rb[1] < 0.0f || rb[1] > (float) userCam->h) {
					rerrRedun += 1.0f;
					dispr->Set(u, v, 0xFF00FFFF);
				}
				float rerr = (V3(r[0], r[1], 0.0f) - V3(rb[0], rb[1], 0.0f)).Length();
				rerrAvg += rerr;
				if (rerrMax < rerr)
					rerrMax = rerr;

				// uv on truth display
				V3 b;
				b[0] = (float)u + .5f;
				b[1] = (float)v + .5f;
				V3 sceneIntTruth;
				if (!IntersetRayWithPlane(userCam->C,
					userCam->GetImagePlanePoint(b[0], b[1]) - userCam->C,
					proxyPoint + proxyNormal*piDelta, proxyNormal, sceneIntTruth)) {
					cerr << "INTERNAL ERROR: no intersection with scene plane truth" << endl;
					exit(0);
				}
				V3 m;  tabletCam->Project(sceneIntTruth, m);
				if (m[0] < 0.0f || m[0] > (float) tabletCam->w ||
					m[1] < 0.0f || m[1] > (float) tabletCam->h) {
					rerrMiss += 1.0f;
					gerrMiss += 1.0f;
					if (dispb->Get(u, v) != 0xFF0000FF && dispb->Get(u, v) != 0xFF00FF00)
						dispb->Set(u, v, 0xFFFFFF00);
				}
				else {
					V3 proxyInt;
					if (!IntersetRayWithPlane(tabletCam->C,
						sceneIntTruth - tabletCam->C,
						proxyPoint, proxyNormal, proxyInt)) {
						cerr << "INTERNAL ERROR: no intersection with proxy plane" << endl;
						exit(0);
					}
					userCam->Project(proxyInt, g);
					if (g[0] < 0.0f || g[0] > (float) userCam->w ||
						g[1] < 0.0f || g[1] > (float) userCam->h) {
						gerrMiss += 1.0f;
						dispb->Set(u, v, 0xFFFF0000);
					}
				}
			}
		}
		dispb->redraw();

		gerrAvg /= (float)(userCam->w * userCam->h);
		rerrAvg /= (float)(userCam->w * userCam->h);
		cerr << "delta = " << piDelta << " avg transpProxy conv: " << gerrAvg << " " << rerrAvg << " Max transpProxy conv: " << gerrMax << " " << rerrMax << endl;
		ofs << piDelta << " " << gerrAvg << " " << rerrAvg << " ";
		ofs << gerrMax << " " << rerrMax << " ";
		if (gerrAvg < rerrAvg && !avgGotBetter) {
			cerr << "INFO: average got better" << endl;
			avgGotBetter = 1;
		}
		if (gerrMax < rerrMax && !maxGotBetter) {
			cerr << "INFO: max got better" << endl;
			maxGotBetter = 1;
		}

		gerrMiss /= (float)(userCam->w * userCam->h);
		rerrMiss /= (float)(userCam->w * userCam->h);
		cerr << "delta = " << piDelta << " missing transpProxy conv: " << gerrMiss * 100.0f << "% "
			<< rerrMiss * 100.0f << "%" << endl;
		ofs << gerrMiss << " " << rerrMiss << " ";

		gerrRedun /= (float)(userCam->w * userCam->h);
		rerrRedun /= (float)(userCam->w * userCam->h);
		cerr << "d = " << d << " redundant transpProxy conv: " << gerrRedun * 100.0f << "% "
			<< rerrRedun * 100.0f << "%" << endl;
		ofs << gerrRedun << " " << rerrRedun << endl;

		Fl::check();
		if ((si % 10) == 0) {
			sprintf_s(fname, "%s/VisConventional_delta%d.tif", folder, (int)piDelta);
			conv->SaveAsTiff(fname);
			sprintf_s(fname, "%s/VisTransparent_delta%d.tif", folder, (int)piDelta);
			transp->SaveAsTiff(fname);
			sprintf_s(fname, "%s/VisTruth_delta%d.tif", folder, (int)piDelta);
			truth->SaveAsTiff(fname);

			sprintf_s(fname, "%s/DisplayConventional_delta%d.tif", folder, (int)piDelta);
			dispr->SaveAsTiff(fname);
			sprintf_s(fname, "%s/DisplayTransparent_delta%d.tif", folder, (int)piDelta);
			dispg->SaveAsTiff(fname);
			sprintf_s(fname, "%s/DisplayTruth_delta%d.tif", folder, (int)piDelta);
			dispb->SaveAsTiff(fname);

		}
	}
	ofs.close();

	return;

}

