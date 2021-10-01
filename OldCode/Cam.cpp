#include "stdafx.h"

#include "cam.h"
#include "m33.h"


PPC2::PPC2(float hfov, int _w, int _h) {


	w = _w;
	h = _h;
	a = V3(1.0f, 0.0f, 0.0f);
	b = V3(0.0f, -1.0f, 0.0f);
	C = V3(0.0f, 0.0f, 0.0f);
	c = V3(-(float)w / 2.0f, (float)h / 2.0f, -(float)w / 2.0f / tan(hfov /180.0f * 3.1415926f / 2.0f));

}

int PPC2::Project(V3 P, V3& Pp) {

	M33 M;
	M.SetColumn(0, a);
	M.SetColumn(1, b);
	M.SetColumn(2, c);

	V3 q = M.Inverted()*(P - C);
	if (q[2] <= 0.0f)
		return 0;

	q[0] /= q[2];
	q[1] /= q[2];

	Pp = q;
	return 1;
}

void PPC2::Pan(float pan) {

	a = a.RotateDirectionAboutDirection(b.UnitVector()*-1.0f, pan);
	c = c.RotateDirectionAboutDirection(b.UnitVector()*-1.0f, pan);

}


void PPC2::Tilt(float tilt) {

	b = b.RotateDirectionAboutDirection(a.UnitVector(), tilt);
	c = c.RotateDirectionAboutDirection(a.UnitVector(), tilt);

}

int PPC2::InsideImage(V3 pp) {

	if (pp[0] < 0.0f || pp[0] >= (float)w || pp[1] < 0.0f ||
		pp[1] >= (float)h)
		return 0;

	return 1;

}

V3 PPC2::GetRay(float uf, float vf) {

	return a*uf + b*vf + c;

}