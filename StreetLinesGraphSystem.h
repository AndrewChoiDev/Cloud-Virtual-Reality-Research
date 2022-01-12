#pragma once
#include <optional>
#include "StreetLines.h"
#include "PosGraphMatrix.h"
#include "Path.h"

const char STREET_LINES_MODE = 0;
const char GRAPH_STRUCTURE_MODE = 1;
const char PATH_STRUCTURE_MODE = 2;
class StreetLinesGraphSystem
{

private:
	StreetLines streetLines;
	PosGraphMatrix posGraphMatrix;
	std::optional<V3> constructionLineStart;
	std::optional<V3> constructionShortestPathStart;
	char mode = STREET_LINES_MODE;
public:
	StreetLinesGraphSystem() {
		this->path.speed = 10.0f;
	}
	void initConstructionLine(V3 pos);
	void endConstructionLine(V3 pos);

	bool constructionLineExists();
	void SWRender(FrameBuffer* fb, PPC* ppc);
	void handleKey(char key, int& needRender);
	void handleMousePressPoint(V3 p);
	Path path = {};
};

