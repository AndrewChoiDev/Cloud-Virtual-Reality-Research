#include "stdafx.h"
#include "StreetLinesGraphSystem.h"
#include <iostream>

void StreetLinesGraphSystem::initConstructionLine(V3 pos)
{
	this->constructionLineStart = pos;
}

void StreetLinesGraphSystem::endConstructionLine(V3 pos)
{
	if (this->constructionLineStart.has_value()) {
		streetLines.AddLine(
			this->constructionLineStart.value(), pos);

		constructionLineStart.reset();
	}
}

bool StreetLinesGraphSystem::constructionLineExists()
{
	return constructionLineStart.has_value();
}

void StreetLinesGraphSystem::SWRender(FrameBuffer* fb, PPC* ppc)
{
	if (this->mode == STREET_LINES_MODE) {
		this->streetLines.SWRender(fb, ppc);
		if (constructionLineExists()) {
			fb->Draw3DPoint(this->constructionLineStart.value(), ppc, 0xFFFF0000, 7);
		}
	}
	if (this->mode == GRAPH_STRUCTURE_MODE || this->mode == PATH_STRUCTURE_MODE) {
		this->posGraphMatrix.SWRender(fb, ppc);
	}
	if (this->mode == PATH_STRUCTURE_MODE) {
		this->path.Render(fb, ppc);

		if (this->constructionShortestPathStart.has_value()) {
			fb->Draw3DPoint(this->constructionShortestPathStart.value(), 
				ppc, 0xFFFF0000, 7);
		}
	}
	
}

void StreetLinesGraphSystem::handleKey(char key, int& needRender)
{
	switch (key) {
		// toggle mode
		case '7': {
			if (this->mode == STREET_LINES_MODE) {
				this->mode = GRAPH_STRUCTURE_MODE;
			}
			else if (this->mode == GRAPH_STRUCTURE_MODE) {
				this->mode = PATH_STRUCTURE_MODE;
			}
			else {
				this->mode = STREET_LINES_MODE;
			}
			std::cerr << "INFO: Graph System Mode Toggled" << endl;
			needRender = 1;
			break;
		}
		case '8': {
			if (this->mode == GRAPH_STRUCTURE_MODE) {
				this->streetLines.SetGraphStructure(this->posGraphMatrix);
			}
			else if (this->mode == PATH_STRUCTURE_MODE) {
				this->path = {};
				path.speed = 10.0f;
			}

			needRender = 1;
			std::cerr << "INFO: Graph System Action Performed" << endl;
			break;
		}
		//default: break;
	}
}

void StreetLinesGraphSystem::handleMousePressPoint(V3 p)
{
	p = p + V3(0.0f, 10.f, 0.0f);
	switch (this->mode) {
		case STREET_LINES_MODE: {
			if (constructionLineExists()) {
				endConstructionLine(p);
			}
			else {
				initConstructionLine(p);
			}
			break;
		}
		case GRAPH_STRUCTURE_MODE: {
			break;
		}
		case PATH_STRUCTURE_MODE: {
			auto nearestNode = posGraphMatrix.nearestNode(p);
			if (nearestNode.has_value()) {
				this->path.AddGroundNode(nearestNode.value());
				/*if (this->constructionShortestPathStart.has_value()) {

					auto potentialNewPath = this->posGraphMatrix.DijkstraShortestPath(
						this->constructionShortestPathStart.value(),
						nearestNode.value()
					);

					if (potentialNewPath.has_value()) {
						this->path = potentialNewPath.value();
					}

					this->constructionShortestPathStart.reset();
				}
				else {
					this->constructionShortestPathStart = nearestNode.value();
				}*/
			}
			break;
		}
							  
	}

		
}

