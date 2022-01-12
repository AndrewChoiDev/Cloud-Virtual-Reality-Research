#include "stdafx.h"
#include "PosGraphMatrix.h"
#include <algorithm>

void PosGraphMatrix::Clear()
{
	this->nodes.clear();
	this->adjacencyMatrix.reset(new char[0]);
}

void PosGraphMatrix::AddNode(V3 pos)
{
	// no duplicates allowed in graph
	for each (auto node in nodes)
	{
		if (node == pos) {
			return;
		}
	}

	nodes.push_back(pos);
	auto oldWidth = nodes.size() - 1;
	auto newWidth = nodes.size();

	auto newMatrix = new char[newWidth * newWidth];

	for (int i = 0; i < oldWidth; i++) {
		for (int j = 0; j < oldWidth; j++)
		{
			newMatrix[i * newWidth + j] = adjacencyMatrix[i * oldWidth + j];
		}
	}
	// set not connected to new entries in array
	for (int l = 0; l < newWidth; l++) {
		newMatrix[l * newWidth + (newWidth - 1)] = NOT_CONNECTED;
		newMatrix[(newWidth - 1) * newWidth + l] = NOT_CONNECTED;
	}

	adjacencyMatrix.reset(newMatrix);
}

void PosGraphMatrix::AddEdge(V3 posA, V3 posB)
{
	// get the two indices of the elements that match posA and B
	auto findA = std::find(nodes.begin(), nodes.end(), posA);
	auto findB = std::find(nodes.begin(), nodes.end(), posB);

	// elements could not be found
	if (findA == nodes.cend() || findB == nodes.cend()) {
		return;
	}
	auto indexA = std::distance(nodes.begin(), findA);
	auto indexB = std::distance(nodes.begin(), findB);

	auto width = this->nodes.size();

	adjacencyMatrix[indexA * width + indexB] = CONNECTED;
	//adjacencyMatrix[indexB * width + indexA] = CONNECTED;
}

void PosGraphMatrix::SWRender(FrameBuffer* fb, PPC* ppc)
{
	for each (auto node in this->nodes) {
		fb->Draw3DPoint(node, ppc, 0xFF0000FF, 5);
	}

	auto len = this->nodes.size();

	V3 c(3.0f, 0.0f, 0.0f);
	V3 up(0.0f, 10.0f, 0.0f);

	for (int i = 0; i < len; i++) {
		for (int j = i + 1; j < len; j++) {
			if (connected(i, j)) {
				fb->Draw3DSegment(this->nodes[i] + up, 
					this->nodes[j] + up, ppc, c, c);
			}
		}
	}
}

bool PosGraphMatrix::ContainsNode(V3 pos)
{
	for each (auto node in this->nodes) {
		if (node == pos) {
			return true;
		}
	}
	return false;
}



std::optional<V3> PosGraphMatrix::nearestNode(V3 pos)
{
	if (this->nodes.size() == 0) {
		return {};
	}

	auto minNode = this->nodes[0];

	for each (auto node in this->nodes)
	{
		if ((node - pos).Length() < (minNode - pos).Length()) {
			minNode = node;
		}
	}

	return minNode;
}

#include <limits>

std::optional<Path> PosGraphMatrix::DijkstraShortestPath(V3 source, V3 target)
{
	auto sourceIndex = nodeIndex(source);
	auto targetIndex = nodeIndex(target);

	// positions dont even exist in this graph
	// or positions are same
	if (targetIndex == -1 || sourceIndex == -1 || 
		targetIndex == sourceIndex) {
		return {};
	}


	auto size = this->nodes.size();
	auto dist = new float[size];
	auto prev = new int[size];

	dist[sourceIndex] = 0.0f;
	prev[sourceIndex] = -1;

	// contains indices of vertices that have yet to be checked
	std::vector<int> remainingVertices = {};

	for (int i = 0; i < this->nodes.size(); i++) {
		if (!(nodes[i] == source)) {
			dist[i] = std::numeric_limits<float>::infinity();
			prev[i] = -1;
		}

		remainingVertices.push_back(i);
	}

	while (remainingVertices.size() > 0) {
		
		// find index of minimum value in dist[]
		auto min = 0;
		{
			auto minIndex = 0;
			for (int i = 1; i < remainingVertices.size(); i++) {
				if (dist[remainingVertices[i]] < dist[remainingVertices[minIndex]]) {
					minIndex = i;
				}
			}

			min = remainingVertices[minIndex];

			// removes element from vector
			remainingVertices.erase(remainingVertices.begin() + minIndex);
		}

		for (int i = 0; i < this->nodes.size(); i++) {
			if (this->connected(min, i)) {
				auto candidateDist = dist[min] + (nodes[min] - nodes[i]).Length();
				if (candidateDist < dist[i]) {
					dist[i] = candidateDist;
					prev[i] = min;
				}
			}
		}
	}

	if (prev[targetIndex] == -1) {
		return {};
	}

	std::vector<int> pathNodes = {};
	int insertValue = targetIndex;
	while (insertValue != -1) {
		pathNodes.push_back(insertValue);
		insertValue = prev[insertValue];
	}


	auto path = Path();
	path.speed = 10.0f;
	while (pathNodes.size() > 0) {
		int nodeIndex = pathNodes[pathNodes.size() - 1];
		pathNodes.pop_back();
		path.AddNode(this->nodes[nodeIndex]);
	}

	return path;
}

int PosGraphMatrix::nodeIndex(V3 pos)
{
	for (int i = 0; i < this->nodes.size(); i++) {
		if (nodes[i] == pos) {
			return i;
		}
	}

	return -1;
}

bool PosGraphMatrix::connected(int u, int v)
{
	auto width = this->nodes.size();
	return adjacencyMatrix[u * width + v] == CONNECTED
		|| adjacencyMatrix[v * width + u] == CONNECTED;
}
