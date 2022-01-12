
#pragma once

#include <ostream>
#include <istream>

using namespace std;

class V3 {
public:
	V3() {};
	V3(float x, float y, float z) { xyz[0] = x; xyz[1] = y; xyz[2] = z; };
	float xyz[3];
	float& operator[](int i);
	V3 operator+(V3 v1);
	V3 operator-(V3 v1);
	friend bool operator== (const V3& a, const V3& b);
	// cerr << v;
	friend ostream& operator<<(ostream& ostr, V3 v);
	friend istream& operator>>(istream& istr, V3 &v);
	float operator*(V3 v1);
	V3 operator*(float scf);
	float Length();
	V3 Normalized();
	V3 operator/(float scf);
	V3 operator^(V3 v2);
	void SetFromColor(unsigned int color);
	unsigned int GetColor();
	V3 RotatePoint(V3 aO, V3 adir, float theta);
	V3 RotateVector(V3 adir, float theta);
	V3 Light(V3 lv, V3 nv, float ka);
	V3 Reflect(V3 ray);
	bool areValuesfinite();

	static V3 indexToRGB(int index);
	static int RGBToIndex(V3 rgb);
};
