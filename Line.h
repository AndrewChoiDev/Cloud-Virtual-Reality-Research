#pragma once

#include "V3.h"
class Line
{
public:
	V3 start;
	V3 end;
	Line(V3 start, V3 end) { this->start = start; this->end = end; };
};

