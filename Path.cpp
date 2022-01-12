#include "stdafx.h"

#include "Path.h"

#include <iostream>
#include <fstream>
#include "TMesh.h"

void Path::AddGroundNode(V3 newGroundNode) {

	V3 newNode = newGroundNode + upGuidance * height;
	AddNode(newNode);

}

void Path::AddNode(V3 newNode) {

	if (!nodes) {
		nodes = new V3[1000];
		nodesN = 0;
		allocatedNodesN = 1000;
	}

	if (nodesN == allocatedNodesN) {
		cerr << "INFO: out of nodes" << endl;
		return;
	}

	nodes[nodesN] = newNode;
	nodesN++;

}


V3 Path::GetSegmentDirection(float timeStep) {

	float prevTime = 0.0f; // time needed to cover all previous path segments
	for (int ni = 0; ni < nodesN - 1; ni++) {
		float segmentTime = (nodes[ni + 1] - nodes[ni]).Length() / speed; // time needed to cover current path segment
		if (prevTime + segmentTime > timeStep) {
			V3 ret = (nodes[ni + 1] - nodes[ni]).Normalized();
			return ret;
		}
		prevTime += segmentTime; // move to next time
	}
	cerr << "INFO: requested time step is beyond end of path";
	return V3(1.0f, 0.0f, 0.0f);

}

void Path::Save(char* fname)
{
	ofstream ofs(fname);
	ofs << this->nodesN << endl;
	for (int i = 0; i < this->nodesN; i++) {
		ofs << this->nodes[i];
	}
	ofs << this->speed << endl;
	ofs << this->upGuidance;
	ofs << this->height << endl;
	ofs.flush();
	ofs.close();
}

void Path::Load(char* fname)
{
	ifstream ifs(fname);
	if (ifs.fail()) {
		cerr << "INFO: cannot open file: " << fname << endl;
		return;
	}

	// node count must be read first to allocate for array
	ifs >> this->nodesN;
	this->allocatedNodesN = 1000;
	this->nodes = new V3[this->allocatedNodesN];

	for (int i = 0; i < this->nodesN; i++) {
		ifs >> this->nodes[i];
	}
	ifs >> this->speed;
	ifs >> this->upGuidance;
	ifs >> this->height;
	ifs.close();
}

// returns position and view direction along path at total play back time "timeStep"
void Path::GetCurrentPositionAndVD(float timeStep, V3& currentPosition,
	V3& currentVD) {

	float prevTime = 0.0f; // time needed to cover all previous path segments
	for (int ni = 0; ni < nodesN - 1; ni++) {
		float segmentTime = (nodes[ni + 1] - nodes[ni]).Length() / speed; // time needed to cover current path segment
		if (prevTime + segmentTime > timeStep) {
			float segmentFractionTime = timeStep - prevTime; // desired time point on current path segment
			float frac = segmentFractionTime / segmentTime;
			currentPosition = nodes[ni] + (nodes[ni + 1] - nodes[ni]) * frac;
			currentVD = (nodes[ni + 1] - nodes[ni]).Normalized();
			float turnTime = 1.0f;
			if (segmentTime - segmentFractionTime < turnTime && ni < nodesN - 2) {
				float turnFrac = 1.0f - (segmentTime - segmentFractionTime) / turnTime;
				if (segmentTime < turnTime)
					turnFrac = segmentFractionTime / segmentTime;
				V3 nextVD = (nodes[ni + 2] - nodes[ni + 1]).Normalized();
				currentVD = currentVD + (nextVD - currentVD) * turnFrac;
				currentVD = currentVD.Normalized();
			}
			return;
		}
		prevTime += segmentTime; // move to next time
	}
	cerr << "INFO: requested time step is beyond end of path";
	currentPosition = nodes[nodesN - 1];

}

void Path::SetCamera(PPC* ppc0, PPC* ppc, float timeStep) {

	*ppc = *ppc0;
	if (speed < 0.0f) {
		cerr << "INFO: set path speed first";
		return;
	}
	V3 currentPosition, currentVD;
	GetCurrentPositionAndVD(timeStep, currentPosition, currentVD);
	*ppc = *ppc0;
	ppc->SetPose(currentPosition, currentPosition + currentVD, upGuidance);

}

float Path::GetTotalTime() {

	float ret = 0;
	for (int ni = 0; ni < nodesN-1; ni++) {
		ret += (nodes[ni + 1] - nodes[ni]).Length()/speed;
	}
	return ret;
}

void Path::Render(FrameBuffer* fb, PPC* ppc) {

	V3 c(0.0f, 0.0f, 1.0f);
	for (int ni = 0; ni < nodesN - 1; ni++) {
		fb->Draw3DSegment(nodes[ni], nodes[ni + 1], ppc, c, c);
	}
	for (int ni = 0; ni < nodesN; ni++) {
		fb->Draw3DPoint(nodes[ni], ppc, 0xFF0000FF, 13);
	}

}


void Path::accumulateVisTrisOnPath(FrameBuffer* hwfb, PPC* ppc, float fps, TMesh* tMesh)
{
	tMesh->ClearVisibleTriangles();

	for (int i = 0; i < this->nodesN - 1; i++) {
		accumulateVisTrisOnSegment(hwfb, ppc, fps, tMesh, i);
	}
}

void Path::accumulateVisTrisOnSegment(FrameBuffer* hwfb, PPC* ppc, float fps, TMesh* tMesh, int nodeStart)
{
	// node start can not be larger than the second-to-last node
	if (nodeStart > this->nodesN - 2 || this->nodesN < 2 || nodeStart < 0) {
		cerr << "ERROR: nodeStart or nodesN size invalid" << endl;
		return;
	}

	// use accumulation of node lengths to find the frame to start rendering and the frame to stop
	auto timeStart = 0.0f;
	for (int i = 0; i < nodeStart; i++) {
		timeStart += (nodes[i + 1] - nodes[i]).Length() / speed;
	}
	auto timeEnd = timeStart + (nodes[nodeStart + 1] - nodes[nodeStart]).Length() / speed;

	auto frameStart = (int)(timeStart * fps);
	auto frameEnd = (int)(timeEnd * fps);

	// rendering and triangle accumulation
	PPC ppc0(*ppc);
	for (int fi = frameStart; fi < frameEnd; fi++) {
		this->SetCamera(ppc, ppc, (float)fi / fps);
		hwfb->redraw();
		Fl::check();
		tMesh->addVisibleTrianglesHWFrameBuffer(hwfb);
	}
	*ppc = ppc0;
	hwfb->redraw();
	Fl::check();
}

void Path::accumulateVisTrisOnInterval(FrameBuffer* hwfb, PPC* ppc, float fps, TMesh* tMesh, int frameStart, int frameEnd)
{
	PPC ppc0(*ppc);
	for (int fi = frameStart; fi < frameEnd; fi++) {
		this->SetCamera(ppc, ppc, (float)fi / fps);
		hwfb->redraw();
		Fl::check();
		tMesh->addVisibleTrianglesHWFrameBuffer(hwfb);
	}
	*ppc = ppc0;
}



// frame points that may be the start or end of an interval
std::vector<int> Path::pathIntervalCriticalFramePoints(float tPeriod, float fps)
{
	auto totalTime = this->GetTotalTime();
	auto intervalCount = (int)ceil(totalTime / tPeriod);

	auto criticalPoints = std::vector<int>();

	for (int ci = 0; ci < intervalCount; ci++) {
		criticalPoints.push_back((int)((float)ci * tPeriod * fps));
	}
	criticalPoints.push_back((int)(totalTime * fps));
	return criticalPoints;
}

#include <vector>

std::vector<TMesh> Path::accumulateVisTrisOnConstantIntervals(FrameBuffer* hwfb, PPC* ppc, float fps, TMesh* tMesh, float tPeriod)
{
	auto totalTime = this->GetTotalTime();
	auto intervalCount = (int)ceil(totalTime / tPeriod);

	auto criticalPoints = pathIntervalCriticalFramePoints(tPeriod, fps);

	auto tMeshChunks = std::vector<TMesh>();

	// iteration for each interval
	for (int i = 0; i < intervalCount; i++) {
		tMesh->ClearVisibleTriangles();
		accumulateVisTrisOnInterval(hwfb, ppc, fps, tMesh, criticalPoints[i], criticalPoints[i + 1]);
		tMeshChunks.push_back(tMesh->constructVisibleMesh());
	}
	tMesh->ClearVisibleTriangles();

	return tMeshChunks;
}

std::vector<TMesh> Path::accumulateNewlyVisTrisOnConstantIntervals(FrameBuffer* hwfb, PPC* ppc, float fps, TMesh* tMesh, float tPeriod)
{
	auto totalTime = this->GetTotalTime();
	auto intervalCount = (int)ceil(totalTime / tPeriod);

	auto criticalPoints = pathIntervalCriticalFramePoints(tPeriod, fps);

	auto tMeshChunks = std::vector<TMesh>();

	tMesh->ClearVisibleTriangles();

	// iteration for each interval
	for (int i = 0; i < intervalCount; i++) {
		auto visTrisBeforeAcc = tMesh->visTrisCopy();

		accumulateVisTrisOnInterval(hwfb, ppc, fps, tMesh, criticalPoints[i], criticalPoints[i + 1]);
		auto visTrisAfterAcc = tMesh->visTrisCopy();
		
		
		// vistris that have been marked visible already must be deleted
		tMesh->removeVisTris(visTrisBeforeAcc);
		tMeshChunks.push_back(tMesh->constructVisibleMesh());
		tMesh->setVisTris(visTrisAfterAcc);
	}
	tMesh->ClearVisibleTriangles();

	return tMeshChunks;
}
