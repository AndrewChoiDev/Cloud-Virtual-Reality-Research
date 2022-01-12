#pragma once
#include <vector>
#include "Line.h"
#include "framebuffer.h"
#include "ppc.h"
#include <optional>

#include "PosGraphMatrix.h"

struct IndexedIntersection {
	V3 intersect;
	int i;
	int j;
};

class StreetLines
{
private:
	std::vector<Line> lines;

	std::optional<V3> lines2DGroundIntersect(Line a, Line b);
	float crossMag(float aX, float aZ, float bX, float bZ);
	std::vector<IndexedIntersection> getIntersections();
	std::vector<Line> getEdges(std::vector<IndexedIntersection> intersections);
public:
	void AddLine(V3 a, V3 b);
	void SWRender(FrameBuffer* fb, PPC* ppc);

	void SetGraphStructure(PosGraphMatrix& graphMat);

};