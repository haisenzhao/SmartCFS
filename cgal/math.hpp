#ifndef math_hpp
#define math_hpp

#include "stdafx.h"
#include <vector>
#include <set>
#include <stdexcept>
#include <cstring>
#include <cstdlib>
#include <ostream>
#include <functional>
#include <queue>
#include <sstream>
#include <math.h>
#include <fstream>
#include <io.h>
#include <direct.h>

typedef glm::highp_vec2 Vector2d;
typedef std::vector<Vector2d> Vector2d1;
typedef std::vector<std::vector<Vector2d>> Vector2d2;
typedef std::vector<std::vector<std::vector<Vector2d>>> Vector2d3;

typedef glm::highp_vec3 Vector3d;
typedef std::vector<Vector3d> Vector3d1;
typedef std::vector<std::vector<Vector3d>> Vector3d2;
typedef std::vector<std::vector<std::vector<Vector3d>>> Vector3d3;

typedef glm::highp_ivec2 Vector2i;
typedef glm::highp_ivec3 Vector3i;

template <typename datum>
using Vector1 = std::vector<datum>;

template <typename datum>
using Vector2 = std::vector<std::vector<datum>>;

template <typename datum>
using Vector3 = std::vector<std::vector<std::vector<datum>>>;

//typedef glm::highp_dvec2 Vector2d;
//typedef std::vector<Vector2d> Vector2d1;
//typedef std::vector<std::vector<Vector2d>> Vector2d2;
//typedef std::vector<std::vector<std::vector<Vector2d>>> Vector2d3;

//typedef glm::highp_dvec3 Vector3d;
//typedef std::vector<Vector3d> Vector3d1;
//typedef std::vector<std::vector<Vector3d>> Vector3d2;
//typedef std::vector<std::vector<std::vector<Vector3d>>> Vector3d3;

typedef std::vector<bool> Vector1b1;
typedef std::vector<std::vector<bool>> Vector1b2;
typedef std::vector<std::vector<std::vector<bool>>> Vector1b3;

typedef std::vector<int> Vector1i1;
typedef std::vector<std::vector<int>> Vector1i2;
typedef std::vector<std::vector<std::vector<int>>> Vector1i3;

typedef std::vector<double> Vector1d1;
typedef std::vector<std::vector<double>> Vector1d2;
typedef std::vector<std::vector<std::vector<double>>> Vector1d3;

typedef std::vector<std::string> VectorStr1;
typedef std::vector<std::vector<std::string>> VectorStr2;
typedef std::vector<std::vector<std::vector<std::string>>> VectorStr3;

//typedef glm::highp_ivec2 Vector2i;
//typedef glm::highp_ivec3 Vector3i;

typedef std::vector<Vector2i> Vector2i1;
typedef std::vector<std::vector<Vector2i>> Vector2i2;
typedef std::vector<std::vector<std::vector<Vector2i>>> Vector2i3;

typedef std::vector<Vector3i> Vector3i1;
typedef std::vector<std::vector<Vector3i>> Vector3i2;
typedef std::vector<std::vector<std::vector<Vector3i>>> Vector3i3;

typedef std::vector<std::pair<int, int>> VectorPI1;
typedef std::vector<std::vector<std::pair<int, int>>> VectorPI2;

typedef std::tuple<int, int, int> TI3;
typedef std::vector<std::tuple<int, int, int>> VectorTI3;



namespace Math {

	static const double DOUBLE_EPSILON = 1.0E-05;
	static const double Math_PI = 3.14159265359;

	inline std::string IntString(int i)
	{
		std::stringstream ss;
		std::string str;
		ss << i;
		ss >> str;
		return str;
	}
	inline std::string DoubleString(double d)
	{
		std::stringstream ss;
		ss.precision(5);
		std::string str;
		ss << d;
		ss >> str;
		return str;
	}

	inline double GetLength(Vector3d v){
		return glm::length(v);
	}

	inline double GetLength(Vector2d v){
		return glm::length(v);
	}

	inline double GetLength(const Vector2d v0, const Vector2d v1) {
		return GetLength(v0 - v1);
	}

	inline Vector2d Vector3d2d(Vector3d v)
	{
		return Vector2d(v[0], v[1]);
	}

	inline Vector2d1 Vector3d2d(Vector3d1 vecs_3d)
	{
		Vector2d1 vecs_2d;
		for (auto v : vecs_3d)
			vecs_2d.emplace_back(Vector3d2d(v));
		return vecs_2d;
	}

	inline Vector2d2 Vector3d2d(Vector3d2 vecs_3ds)
	{
		Vector2d2 vecs_2d;
		for (auto vecs_3d : vecs_3ds)
			vecs_2d.emplace_back(Vector3d2d(vecs_3d));
		return vecs_2d;
	}


	inline Vector3d Vector2d3d(Vector2d v, double z = 0.0)
	{
		return Vector3d(v[0], v[1], z);
	}

	inline Vector3d1 Vector2d3d(Vector2d1 vecs_2d, double z = 0.0)
	{
		Vector3d1 vecs_3d;
		for (auto v : vecs_2d)
			vecs_3d.emplace_back(Vector2d3d(v, z));
		return vecs_3d;
	}

	inline Vector3d2 Vector2d3d(Vector2d2 vecs_2d, double z = 0.0)
	{
		Vector3d2 vecs_3d;
		for (auto v : vecs_2d)
			vecs_3d.emplace_back(Vector2d3d(v, z));
		return vecs_3d;
	}
	


	inline Vector3d RotationAxis(Vector3d p, double angle, Vector3d n)
	{
		glm::mat4 inputMatrix(0.0);
		inputMatrix[0][0] = p[0];
		inputMatrix[1][0] = p[1];
		inputMatrix[2][0] = p[2];
		inputMatrix[3][0] = 1.0;
		double u = n[0];
		double v = n[1];
		double w = n[2];

		glm::mat4  rotationMatrix;

		double L = (u*u + v*v + w*w);

		//angle = angle * M_PI / 180.0; //converting to radian value
		double u2 = u * u;
		double v2 = v * v;
		double w2 = w * w;

		rotationMatrix[0][0] = (u2 + (v2 + w2) * cos(angle)) / L;
		rotationMatrix[0][1] = (u * v * (1 - cos(angle)) - w * sqrt(L) * sin(angle)) / L;
		rotationMatrix[0][2] = (u * w * (1 - cos(angle)) + v * sqrt(L) * sin(angle)) / L;
		rotationMatrix[0][3] = 0.0;

		rotationMatrix[1][0] = (u * v * (1 - cos(angle)) + w * sqrt(L) * sin(angle)) / L;
		rotationMatrix[1][1] = (v2 + (u2 + w2) * cos(angle)) / L;
		rotationMatrix[1][2] = (v * w * (1 - cos(angle)) - u * sqrt(L) * sin(angle)) / L;
		rotationMatrix[1][3] = 0.0;

		rotationMatrix[2][0] = (u * w * (1 - cos(angle)) - v * sqrt(L) * sin(angle)) / L;
		rotationMatrix[2][1] = (v * w * (1 - cos(angle)) + u * sqrt(L) * sin(angle)) / L;
		rotationMatrix[2][2] = (w2 + (u2 + v2) * cos(angle)) / L;
		rotationMatrix[2][3] = 0.0;

		rotationMatrix[3][0] = 0.0;
		rotationMatrix[3][1] = 0.0;
		rotationMatrix[3][2] = 0.0;
		rotationMatrix[3][3] = 1.0;

		double outputMatrix[4][1] = { 0.0, 0.0, 0.0, 0.0 };

		for (int i = 0; i < 4; i++){
			for (int j = 0; j < 1; j++){
				outputMatrix[i][j] = 0;
				for (int k = 0; k < 4; k++){
					outputMatrix[i][j] += rotationMatrix[i][k] * inputMatrix[k][j];
				}
			}
		}
		return Vector3d(outputMatrix[0][0], outputMatrix[0][1], outputMatrix[0][2]);
	}

	inline Vector3d GetCrossproduct(Vector3d& v1, Vector3d& v2) {
		return glm::cross(v1, v2);
	}

	/*double GetDotProduct(Vector3d& v1, Vector3d& v2){
		return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
	}
*/
	inline Vector3d SetVectorLength(Vector3d &v, double length)
	{
		double l = GetLength(v);

		v[0] = v[0] / l*length;
		v[1] = v[1] / l*length;
		v[2] = v[2] / l*length;

		return v;
	}

	inline bool IsAlmostZero(double value) {
		return value < 10.0 * DOUBLE_EPSILON && value > -10.0 * DOUBLE_EPSILON;
	}

	inline bool IsAlmostZero_Double(double value, double EPSILON) {
		return value < 10.0 * DOUBLE_EPSILON && value > -10.0 *EPSILON;
	}

	inline void ClearVector3d(Vector3d &v)
	{
		if (IsAlmostZero(v[0]))v[0] = 0.0;
		if (IsAlmostZero(v[1]))v[1] = 0.0;
		if (IsAlmostZero(v[2]))v[2] = 0.0;
	}

	inline double GetAngleBetween(Vector3d v1, Vector3d v2) {

		double d = glm::dot(v1, v2) / (glm::length(v1) * glm::length(v2));

		if (IsAlmostZero(d - 1.0))
			return 0.0;

		if (IsAlmostZero(d + 1.0))
			return Math_PI;

		return glm::acos(glm::dot(v1, v2) / (glm::length(v1) * glm::length(v2)));
	}


	inline double GetAngleBetween(Vector2d v1, Vector2d v2) {
		double d = glm::dot(v1, v2) / (glm::length(v1) * glm::length(v2));
		if (IsAlmostZero(d - 1.0))
			return 0.0;
		if (IsAlmostZero(d + 1.0))
			return Math_PI;
		return glm::acos(d);
	}

	inline void ColorMapping(double isolevel, double &output_c_0, double &output_c_1, double &output_c_2)
	{
		Vector3d v;
		if (isolevel >= 0 && isolevel <= 0.25)
		{
			v[0] = 0;
			v[1] = isolevel / 0.25;
			v[2] = 1;
		}

		if (isolevel>0.25&&isolevel <= 0.50)
		{
			v[0] = 0;
			v[1] = 1;
			v[2] = 1 - (isolevel - 0.25) / 0.25;
		}

		if (isolevel>0.50&&isolevel <= 0.75)
		{
			v[0] = (isolevel - 0.50) / 0.25;
			v[1] = 1;
			v[2] = 0;
		}

		if (isolevel>0.75&&isolevel <= 1.0)
		{
			v[0] = 1;
			v[1] = 1 - (isolevel - 0.75) / 0.25;
			v[2] = 0;
		}

		if (isolevel < 0.0)
		{
			v[0] = 0.0;
			v[1] = 0.0;
			v[2] = 0.0;
		}

		if (isolevel > 1.0)
		{
			v[0] = 0.5;
			v[1] = 0.0;
			v[2] = 0.0;
		}
		output_c_0 = v[0];
		output_c_1 = v[1];
		output_c_2 = v[2];
	}

	inline void ClearFolder(const std::string& path)
	{
		if (_access(path.c_str(), 0) == -1)
		{
			_mkdir(path.c_str());
		}
		else
		{
			std::string cmd = "rmdir " + path + ". / s / q";
			system(cmd.c_str());
			_mkdir(path.c_str());
		}
	}

	static void MAssert(std::string& str)
	{
		std::cerr << str << std::endl;
		system("pause");
	}

	static void MAssert(char* str)
	{
		std::cerr << str << std::endl;
		system("pause");
	}
}
#endif 


