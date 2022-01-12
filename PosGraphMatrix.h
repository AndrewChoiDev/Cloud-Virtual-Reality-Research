#pragma once


#include <vector>
#include "V3.h"
#include <memory>
#include "framebuffer.h"
#include "ppc.h"
#include <optional>

#include "Path.h"

const char CONNECTED = 1;
const char NOT_CONNECTED = 0;

class PosGraphMatrix
{


public:
	void Clear();
	void AddNode(V3 pos);
	void AddEdge(V3 posA, V3 posB);
	void SWRender(FrameBuffer* fb, PPC* ppc);
	bool ContainsNode(V3 pos);
	
	std::optional<V3> nearestNode(V3 pos);

	std::optional<Path> DijkstraShortestPath(V3 source, V3 target);


private:
	std::vector<V3> nodes;
	std::unique_ptr<char[]> adjacencyMatrix;

	int nodeIndex(V3 pos);
	bool connected(int u, int v);

};

