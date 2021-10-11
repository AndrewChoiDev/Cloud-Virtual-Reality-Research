#pragma once

#include "V3.h"
#include "ppc.h"
#include "TMesh.h"

class Path {
public:
	int nodesN;
	int allocatedNodesN;
	V3 *nodes;
	float speed; // constant speed along path

	V3 upGuidance; // up direction, perpendicular to path
	float height; // path node height above the clicked point

	Path() : nodesN(0), nodes(0), allocatedNodesN(0), speed(-1.0f),
		upGuidance(0.0f, 1.0f, 0.0f), height(1.0f) {};
	void AddNode(V3 newNode);
	void AddGroundNode(V3 newGroundNode); // raises node then adds it
	void SetCamera(PPC *ppc0, PPC *ppc, float timeStep);
	void GetCurrentPositionAndVD(float timeStep, V3& currentPosition, V3& currentVD);
	V3 GetSegmentDirection(float timeStep);
	void Save(char *fname); // save all path info in a text file, look at PPC save
	void Load(char *fname); // load path from text file
	float GetTotalTime();
	void Render(FrameBuffer* fb, PPC* ppc);
	void accumulateVisTrisOnSegment(FrameBuffer* hwfb, PPC* ppc, float fps, TMesh* tMesh);
};