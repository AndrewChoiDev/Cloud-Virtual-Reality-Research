#include "stdafx.h"
#include "StreetLines.h"
#include <algorithm>
#include "V3.h"
#include <iostream>
std::optional<V3> StreetLines::lines2DGroundIntersect(Line a, Line b)
{
	auto aDeltaX = a.end[0] - a.start[0];
	auto aDeltaZ = a.end[2] - a.start[2];
	auto bDeltaX = b.end[0] - b.start[0];
	auto bDeltaZ = b.end[2] - b.start[2];

	auto deltaCross = crossMag(aDeltaX, aDeltaZ, bDeltaX, bDeltaZ);

	if (abs(deltaCross) < 0.00001f) {
		return {};
	}

	auto startDiffX = b.start[0] - a.start[0];
	auto startDiffZ = b.start[2] - a.start[2];

	auto tA = crossMag(startDiffX, startDiffZ, bDeltaX, bDeltaZ) / deltaCross;
	auto tB = crossMag(startDiffX, startDiffZ, aDeltaX, aDeltaZ) / deltaCross;

	if (tA < 0 || tA > 1 || tB < 0 || tB > 1) {
		return {};
	}

	return a.start + (a.end - a.start) * tA;
}

float StreetLines::crossMag(float aX, float aZ, float bX, float bZ)
{
	return aX * bZ - aZ * bX;
}

std::vector<IndexedIntersection> StreetLines::getIntersections()
{
	auto intersections = std::vector<IndexedIntersection>();


	auto len = this->lines.size();

	V3 c(0.0f, 0.0f, 1.0f);
	V3 up(0.0f, 5.0f, 0.0f);

	for (int i = 0; i < len; i++) {
		for (int j = i + 1; j < len; j++) {
			auto hit = lines2DGroundIntersect(lines[i], lines[j]);
			if (hit.has_value()) {
				intersections.push_back(
					IndexedIntersection{ hit.value(), i, j });
			}
		}
	}

	return intersections;
}

std::vector<Line> StreetLines::getEdges(std::vector<IndexedIntersection> intersections)
{

	auto edges = std::vector<Line>();

	// each entry is a list of intersections
	// that occur on the line i, the index of list's entry
	auto intersectionsPerLineList = 
		std::vector<std::vector<IndexedIntersection>>();

	//initialize each list
	for (int i = 0; i < this->lines.size(); i++) {
		intersectionsPerLineList.push_back(
			std::vector<IndexedIntersection>());
	}

	// add each intersection to lists
	for each (auto inter in intersections) {
		intersectionsPerLineList[inter.i].push_back(inter);
		intersectionsPerLineList[inter.j].push_back(inter);
	}

	for (int i = 0; i < intersectionsPerLineList.size(); i++) {
		auto lineStart = this->lines[i].start;
		std::sort(intersectionsPerLineList[i].begin(),
			intersectionsPerLineList[i].end(),
			[&lineStart](IndexedIntersection& a, IndexedIntersection& b)
			-> bool {
				return (a.intersect - lineStart).Length() < (b.intersect - lineStart).Length();
			});

		for (int j = 0; j < intersectionsPerLineList[i].size() - 1; j++) {
			edges.push_back(
				Line(
					intersectionsPerLineList[i][j].intersect,
					intersectionsPerLineList[i][j+1].intersect
				)
			);
		}
	}

	return edges;
}

void StreetLines::AddLine(V3 a, V3 b)
{
	lines.push_back(Line(a, b));
}

void StreetLines::SWRender(FrameBuffer* fb, PPC* ppc)
{
	V3 c(0.0f, 0.0f, 1.0f);

	V3 up(0.0f, 5.0f, 0.0f);

	for each (auto line in this->lines) {
		fb->Draw3DSegment(line.start + up, line.end + up, ppc, c, c);
	}
}

void StreetLines::SetGraphStructure(PosGraphMatrix& graphMat)
{
	graphMat.Clear();

	auto intersections = getIntersections();
	auto edges = getEdges(intersections);

	std::cerr << "INFO: intersection count: " << intersections.size()
		<< endl;
	std::cerr << "INFO: edge count: " << edges.size()
		<< endl;

	for each (auto inter in intersections)
	{
		graphMat.AddNode(inter.intersect);
	}

	for each (auto edge in edges)
	{
		graphMat.AddEdge(edge.start, edge.end);
		//std::cerr << "INFO: edge added: " << edge.start << ", " << edge.end << std::endl;
	}
}
