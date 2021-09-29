#pragma once

#include <cmath>
#include <math.h>

using namespace std;

namespace hpcg {
	
	static const double DOUBLE_EPSILON = 1.0E-05;
	static const float SINGLE_EPSILON = 1.0E-05f;
    static const double DEG_TO_RAD = 0.0174532925;
    static const double RAD_TO_DEG = 57.2957795;
	static const float  PRACTICAL_SINGLE_EPSILON = 1.0E-05f;
	static const double PRACTICAL_DOUBLE_EPSILON = 1.0E-07;
	static const double MM_PI = 3.14159265358979323846;

	static const double MAXDOUBLE = 1000000000.0;
	static const double TAU = 6.28318530;




	inline std::string MyGetUserName()
	{
		return "nodebug";
	}


	inline std::vector<string> Split(std::string& s, std::string delim)
	{
		std::vector<string> ret;

		s.erase(0, s.find_first_not_of(" "));
		s.erase(s.find_last_not_of(" ") + 1);

		size_t last = 0;
		size_t index = s.find_first_of(delim, last);
		while (index != std::string::npos)
		{
			if (index - last > 0)
				ret.push_back(s.substr(last, index - last));

			last = index + 1;
			index = s.find_first_of(delim, last);
		}
		if (index - last > 0)
		{
			if (index - last > 0) ret.push_back(s.substr(last, index - last));
		}

		return ret;
	}


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

	template <class Type>
	Type StringToNum(const string& str)
	{
		istringstream iss(str);
		Type num;
		iss >> num;
		return num;
	}

	inline bool isAlmostZero(double value) {
		return value < 10.0 * DOUBLE_EPSILON && value > -10.0 * DOUBLE_EPSILON;
	}

	inline bool isAlmostZero(float value) {
		return value < 10.0f * SINGLE_EPSILON && value > -10.0f * SINGLE_EPSILON;
	}
	
	//Practical zero limits
	inline bool isPracticalZero(double value) {
		return value < 10.0 * PRACTICAL_DOUBLE_EPSILON && value > -10.0 * PRACTICAL_DOUBLE_EPSILON;
	}

	inline bool isPracticalZero(float value) {
		return value < 10.0f * PRACTICAL_SINGLE_EPSILON && value > -10.0f * PRACTICAL_SINGLE_EPSILON;
	}

	/// Returns true if two given floating point numbers are epsilon-equal.
	/// Method automatically adjust the epsilon to the absolute size of given numbers.
	inline bool areAlmostEqual(double value1, double value2) {
		// in case they are Infinities (then epsilon check does not work)
		if (value1 == value2) {
			return true;
		}

		// computes (|value1-value2| / (|value1| + |value2| + 10.0)) < DOUBLE_EPSILON
		double eps = (glm::abs(value1) + glm::abs(value2) + 10.0) * DOUBLE_EPSILON;
		double delta = value1 - value2;
		return (-eps < delta) && (eps > delta);

	}
	
	/// Returns true if two given floating point numbers are epsilon-equal.
	/// Method automatically adjust the epsilon to the absolute size of given numbers.
	inline bool areAlmostEqual(float value1, float value2) {
		// in case they are Infinities (then epsilon check does not work)
		if (value1 == value2) {
			return true;
		}

		// computes (|value1-value2| / (|value1| + |value2| + 10.0)) < SINGLE_EPSILON
		double eps = (glm::abs(value1) + glm::abs(value2) + 10.0) * SINGLE_EPSILON;
		double delta = value1 - value2;
		return (-eps < delta) && (eps > delta);

	}

	inline bool isEpsilonGreaterThanZero(float value) {
		return value > 10 * SINGLE_EPSILON;
	}

	inline bool isEpsilonGreaterThanZero(double value) {
		return value > 10 * DOUBLE_EPSILON;
	}

	inline bool isInfinity(float value) {
		return value == std::numeric_limits<float>::infinity() || value == -std::numeric_limits<float>::infinity();
	}
	
	inline bool isInfinity(double value) {
		return value == std::numeric_limits<double>::infinity() || value == -std::numeric_limits<double>::infinity();
	}
	
	inline bool isFiniteNumber(float value) {
		return value == value && !isInfinity(value);
	}
	
	inline bool isFiniteNumber(double value) {
		return value == value && !isInfinity(value);
	}

    inline void SnapToGrid(double *p, double dist) {
        double fd = fmod(*p, dist);
        (fd < 0.5*dist) ? *p -= fd : *p += (dist - fd);
    }

	template<typename GlmVector>
	inline double getAngleBetween(GlmVector v1, GlmVector v2) {

		double d = glm::dot(v1, v2) / (glm::length(v1) * glm::length(v2));

		if (isAlmostZero(d - 1.0))
			return 0.0;

		if (isAlmostZero(d + 1.0))
			return MM_PI;

		return glm::acos(glm::dot(v1, v2) / (glm::length(v1) * glm::length(v2)));
	}

	template<typename GlmVector>
	inline double getLength(GlmVector v){
		return glm::length(v);
	}

	template<typename GlmVector>
	inline double getLength(GlmVector v1, GlmVector v2){
		return glm::length(v1-v2);
	}

    template<typename GlmVector>
    inline double getTriangleArea(GlmVector& v1, GlmVector& v2, GlmVector& v3) {
        return glm::length(glm::cross(v2-v1, v3-v1))/2.0;
    }

	inline double getConvexHullArea(Vector2d1 hull)
	{
		double area = 0.0;
		Vector2d center(0.0,0.0);

		for (int i = 0; i < hull.size(); i++)
		{
			center[0] += hull[i][0];
			center[1] += hull[i][1];
		}

		center[0] = center[0] / hull.size();
		center[1] = center[1] / hull.size();

		for (int i = 0; i < hull.size(); i++)
		{
			Vector3d v0(hull[i][0], hull[i][1],0.0);
			Vector3d v1(hull[(i + 1) % hull.size()][0], hull[(i + 1) % hull.size()][1],0.0);
			Vector3d v2(center[0], center[1],0.0);

			area += getTriangleArea(v0,v1,v2);
			//area += getTriangleArea(hull[i], hull[(i + 1) % hull.size()],center);
		}
		
		return area;
	}


	template<typename GlmVector>
	inline GlmVector getCrossproduct(GlmVector& v1, GlmVector& v2) {
		return glm::cross(v1, v2);
	}

	template<typename GlmVector>
	double getDotProduct(GlmVector& v1, GlmVector& v2)
	{
		return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
	}


    template<typename GlmVector>
    bool byHeight(GlmVector& a, GlmVector& b) {
        return (a.y > b.y);
    }

	//Corrected a terrible bug in Juraj code!!, Yell at him next time you see it!!
    template<typename GlmVector>
    bool IsEqual(GlmVector& a, GlmVector& b) {
        return isAlmostZero(glm::distance2(a, b));
    }

	//A not so precise version of the former to fight against numerical impresition
	template<typename GlmVector>
	bool IsAlmostEqual(GlmVector& a, GlmVector& b) {
		return isPracticalZero(glm::distance2(a, b));
	}


    template<typename GlmVector>
    inline double Barycenter(GlmVector& a, GlmVector& b, GlmVector& c) {
        return (b.x-a.x)*(c.y-a.y) - (b.y-a.y)*(c.x-a.x);
    }

	template<typename GlmVector>
	inline bool CrossTriangle(const GlmVector& v_0, const GlmVector& v_1, const GlmVector& v_2, const GlmVector& point, const GlmVector& ray) {
		glm::mat3 M = glm::mat3(-ray, v_1 - v_0, v_2 - v_0);
		glm::vec3 solution = glm::inverse(M) * (point - v_0);

		if (solution[0] >= 0.0f && 
			solution[1] >= 0.0f && solution[2] >= 0.0f && solution[1] <= 1.0f && solution[2] <= 1.0f && 
			(solution[1] + solution[2]) <= 1.0f) {
			return true;
		}
		return false;
	}

	template<typename GlmVector>
	inline bool TriangleCrossPlane(const GlmVector& v_0, const GlmVector& v_1, const GlmVector& v_2, const double& height) {
		if (v_0.y < height && v_1.y < height && v_2.y < height) {
			return false;
		}
		else if (v_0.y > height && v_1.y > height && v_2.y > height) {
			return false;
		}
		return true;
	}

	inline Vector2d GetPointOnPlane(const Vector3d& p_0, const Vector3d& p_1, const double& height) {
		double t = (height - p_1.y) / (p_0.y - p_1.y);
		return Vector2d(p_0.x * t + (1.0 - t) * p_1.x, p_0.z * t + (1.0 - t) * p_1.z);
	}

	inline double VolumeUnderTriangle(const Vector3d& p1, const Vector3d& p2, const Vector3d& p3) {
		//Sort the point to the Y coordinate
		Vector3d points[3];
		if (p1.y < p2.y) {
			if (p1.y < p3.y) {
				if (p2.y < p3.y) {
					points[0] = p1;
					points[1] = p2;
					points[2] = p3;
				} else {
					points[0] = p1;
					points[1] = p3;
					points[2] = p2;
				}
			} else { //
				points[0] = p3;
				points[1] = p1;
				points[2] = p2;
			}
		} else {
			if (p2.y < p3.y) {
				if (p1.y < p3.y) {
					points[0] = p2;
					points[1] = p1;
					points[2] = p3;
				}
				else {
					points[0] = p2;
					points[1] = p3;
					points[2] = p1;
				}
			}
			else {
				points[0] = p3;
				points[1] = p2;
				points[2] = p1;
			}
		}
		//We have the points sorted by Y direction
		double projection = points[0].x * points[1].z - points[1].x * points[0].z
			              + points[1].x * points[2].z - points[2].x * points[1].z
			              + points[2].x * points[0].z - points[0].x * points[2].z;

		return (points[0].y + points[1].y + points[2].y) * projection / 6.0;
	}

	inline double RadiantoAngle(double r)
	{
		return r / MM_PI*180.0;
	}
	/*	inline glm::mat4 RotationAxis(glm::mat4 inputMatrix, double angle, Vector3d n)
	{

		double u = n[0];
		double v = n[1];
		double w = n[2];

		glm::mat4  rotationMatrix;

		double L = (u*u + v*v + w*w);

		//angle = angle * MM_PI / 180.0; //converting to radian value
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

		glm::mat4 outputMatrix(0.0);

		for (int i = 0; i < 4; i++){
			for (int j = 0; j < 4; j++){
				for (int k = 0; k < 4; k++){
					outputMatrix[i][j] += rotationMatrix[i][k] * inputMatrix[k][j];
				}
			}
		}

		return outputMatrix;
	}*/

	inline Vector3d SetVectorLength(Vector3d &v, double length)
	{
		double l = getLength(v);

		v[0] = v[0] / l*length;
		v[1] = v[1] / l*length;
		v[2] = v[2] / l*length;

		return v;
	}

	inline Vector2d SetVectorLength(Vector2d &v, double length)
	{
		double l = getLength(v);
		v[0] = v[0] / l*length;
		v[1] = v[1] / l*length;
		return v;
	}


	inline Vector3d GetCenter(Vector3d v0, Vector3d v1)
	{
		Vector3d center = v0 + v1;
		
		center[0] = center[0] / 2.0;
		center[1] = center[1] / 2.0;
		center[2] = center[2] / 2.0;

		return center;
	}

	inline Vector3d GetCenter(Vector3d v0, Vector3d v1,Vector3d v2)
	{
		Vector3d center = v0 + v1+v2;

		center[0] = center[0] / 3.0;
		center[1] = center[1] / 3.0;
		center[2] = center[2] / 3.0;

		return center;
	}

	inline Vector3d GetCenter(const std::vector<Vector3d> &vecs)
	{
		Vector3d center(0.0, 0.0, 0.0);
		for (int i = 0; i < vecs.size(); i++) center += vecs[i];
		center = center / (float)vecs.size();
		return center;
	}


	inline Vector2d GetCenter(const std::vector<Vector2d> &vecs)
	{
		Vector2d center(0.0, 0.0);
		for (int i = 0; i < vecs.size(); i++) center += vecs[i];
		center = center / (float)vecs.size();
		return center;
	}

	inline Vector3d GetCrossproduct(const Vector3d& v1, const Vector3d& v2) {
		return glm::cross(v1, v2);
	}

	template <class Type>
	inline double GetAngleBetween(Type v1, Type v2) {
		double d = glm::dot(v1, v2) / (glm::length(v1) * glm::length(v2));
		if (IsAlmostZero(d - 1.0))
			return 0.0;
		if (IsAlmostZero(d + 1.0))
			return MM_PI;
		return glm::acos(d);
	}

	inline bool IsAlmostZero(double value) {
		return (value < DOUBLE_EPSILON) && (value > -DOUBLE_EPSILON);
	}


	inline Vector3d PosApplyM(const Vector3d& v, const glm::dmat4& M)
	{
		return Vector3d(M * glm::vec4(v, 1.0));
	}

	inline std::vector<Vector3d> PosApplyM(const std::vector<Vector3d>& vecs, const glm::dmat4& M)
	{
		std::vector<Vector3d> ps;
		for (auto p : vecs)
			ps.emplace_back(PosApplyM(p, M));
		return ps;
	}

	inline std::pair<Vector3d, Vector3d> PosApplyM(const std::pair<Vector3d, Vector3d> &vecs, const glm::dmat4& M)
	{
		return std::pair<Vector3d, Vector3d>(PosApplyM(vecs.first, M), PosApplyM(vecs.second, M));
	}


	inline Vector3d2 PosApplyM(const Vector3d2& veces, const glm::dmat4& M)
	{
		Vector3d2 pses;
		for (auto vecs : veces)
			pses.emplace_back(PosApplyM(vecs, M));
		return pses;
	}

	inline Vector3d3 PosApplyM(const Vector3d3& veces, const glm::dmat4& M)
	{
		Vector3d3 pses;
		for (auto vecs : veces)
			pses.emplace_back(PosApplyM(vecs, M));
		return pses;
	}


	inline glm::dmat4 RotationMatrix(const Vector3d& n, double angle)
	{
		//return glm::rotate(angle, n);
		double u = n[0];
		double v = n[1];
		double w = n[2];

		glm::dmat4  rotationMatrix;

		double L = (u * u + v * v + w * w);

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

		return rotationMatrix;
	}

	inline void AAA(const Vector3d& v,Vector3d &n)
	{
		auto a = v[0];
		auto b = v[1];
		auto c = v[2];
		bool bx = IsAlmostZero(a);
		bool by = IsAlmostZero(b);
		bool bz = IsAlmostZero(c);

		if (bx&&by&&bz)
		{
			std::cerr << "if (bx&&by&&bz)" << std::endl;
			system("pause");
		}

		if (bx&&by&&!bz)
		{
			n[0] = 1.0;
			n[1] = 1.0;
			n[2] = 0.0;
		}
		if (bx&&!by&&bz)
		{
			n[0] = 1.0;
			n[1] = 0.0;
			n[2] = 1.0;
		}

		if (bx&&!by&&!bz)
		{
			n[0] = 1.0;
			n[1] = 1.0;
			n[2] = -b/c;
		}

		if (!bx&&by&&bz)
		{
			n[0] = 0.0;
			n[1] = 1.0;
			n[2] = 1.0;
		}

		if (!bx&&by&&!bz)
		{
			n[0] = 1.0;
			n[1] = 1.0;
			n[2] = -a / c;
		}
		if (!bx&&!by&&bz)
		{
			n[0] = 1.0;
			n[1] = -a/b;
			n[2] = 1.0;
		}

		if (!bx&&!by&&!bz)
		{
			n[0] = 1.0;
			n[1] = 1.0;
			n[2] = -(a+b) / c;
		}
	
	}

	inline glm::dmat4 RotationMatrixXYZ(const Vector3d& xx, const Vector3d &yy, const Vector3d &zz)
	{
		auto x = xx;
		auto y = yy;
		auto z = zz;

		Math::ClearVector3d(x);
		Math::ClearVector3d(y);
		Math::ClearVector3d(z);
		x = x / (float)Math::GetLength(x);
		y = y / (float)Math::GetLength(y);
		z = z / (float)Math::GetLength(z);

		glm::dmat4  rotationMatrix;

		rotationMatrix[0][0] = x[0];
		rotationMatrix[0][1] = y[0];
		rotationMatrix[0][2] = z[0];
		rotationMatrix[0][3] = 0.0;

		rotationMatrix[1][0] = x[1];
		rotationMatrix[1][1] = y[1];
		rotationMatrix[1][2] = z[1];
		rotationMatrix[1][3] = 0.0;

		rotationMatrix[2][0] = x[2];
		rotationMatrix[2][1] = y[2];
		rotationMatrix[2][2] = z[2];
		rotationMatrix[2][3] = 0.0;

		rotationMatrix[3][0] = 0.0;
		rotationMatrix[3][1] = 0.0;
		rotationMatrix[3][2] = 0.0;
		rotationMatrix[3][3] = 1.0;

		return rotationMatrix;

	}

	inline glm::dmat4 RotationMatrix(const Vector3d& o, const Vector3d &t, const Vector3d &n)
	{
		double angle = GetAngleBetween(o, t);

		if (IsAlmostZero(angle))
		{
			glm::dmat4  rotationMatrix;

			rotationMatrix[0][0] = 1.0;
			rotationMatrix[0][1] = 0.0;
			rotationMatrix[0][2] = 0.0;
			rotationMatrix[0][3] = 0.0;
			rotationMatrix[1][0] = 0.0;
			rotationMatrix[1][1] = 1.0;
			rotationMatrix[1][2] = 0.0;
			rotationMatrix[1][3] = 0.0;
			rotationMatrix[2][0] = 0.0;
			rotationMatrix[2][1] = 0.0;
			rotationMatrix[2][2] = 1.0;
			rotationMatrix[2][3] = 0.0;
			rotationMatrix[3][0] = 0.0;
			rotationMatrix[3][1] = 0.0;
			rotationMatrix[3][2] = 0.0;
			rotationMatrix[3][3] = 1.0;

			return rotationMatrix;
		}
		else
		{
			return RotationMatrix(n, angle);
		}
	}



	inline glm::dmat4 RotationMatrix(const Vector3d& o, const Vector3d &t)
	{
		Vector3d n = GetCrossproduct(o, t);
		double angle = GetAngleBetween(o, t);

		if (IsAlmostZero(angle - MM_PI))
		{
			AAA(o,n);
		}

		if (IsAlmostZero(angle))
		{
			glm::dmat4  rotationMatrix;

			rotationMatrix[0][0] = 1.0;
			rotationMatrix[0][1] = 0.0;
			rotationMatrix[0][2] = 0.0;
			rotationMatrix[0][3] = 0.0;
			rotationMatrix[1][0] = 0.0;
			rotationMatrix[1][1] = 1.0;
			rotationMatrix[1][2] = 0.0;
			rotationMatrix[1][3] = 0.0;
			rotationMatrix[2][0] = 0.0;
			rotationMatrix[2][1] = 0.0;
			rotationMatrix[2][2] = 1.0;
			rotationMatrix[2][3] = 0.0;
			rotationMatrix[3][0] = 0.0;
			rotationMatrix[3][1] = 0.0;
			rotationMatrix[3][2] = 0.0;
			rotationMatrix[3][3] = 1.0;

			return rotationMatrix;
		}
		else
		{
			return RotationMatrix(n, angle);
		}
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

		//angle = angle * MM_PI / 180.0; //converting to radian value
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


	inline std::vector<int> SetUnion(std::vector<int> &a, std::vector<int> &b)
	{
		std::vector<int> v(a.size() + b.size());
		std::vector<int>::iterator it;

		std::sort(a.begin(), a.end());
		std::sort(b.begin(), b.end());

		it = std::set_union(a.begin(), a.end(), b.begin(), b.end(), v.begin());

		v.resize(it - v.begin());

		return v;
	}


	inline std::vector<int> SetUnion(std::vector<std::vector<int>> &sets)
	{
		std::vector<int> start = sets[0];

		for (int i = 1; i < sets.size(); i++)
		{
			start = SetUnion(start,sets[i]);
		}

		return start;
	}

	inline bool VectorContain(std::vector<int> vecs, int element)
	{
		for (int i = 0; i < vecs.size(); i++)
		{
			if (vecs[i] == element)
				return true;
		}

		return false;
	}

	inline int VectorContainReturnIndex(std::vector<int> vecs, int element)
	{
		for (int i = 0; i < vecs.size(); i++)
		{
			if (vecs[i] == element)
				return i;
		}

		return -1;
	}

	inline bool VectorContainForSpecialCase(std::vector<std::vector<int>> vecs, std::vector<int> element)
	{
		for (int i = 0; i < vecs.size(); i++)
		{
			if (vecs[i][0] == element[0] && vecs[i][1] == element[1])
				return true;
		}

		return false;
	}
	inline bool VectorContainForSpecialCase1(std::vector<std::vector<int>> vecs, std::vector<int> element)
	{
		for (int i = 0; i < vecs.size(); i++)
		{
			if (vecs[i][0] == element[0] && vecs[i][1] == element[1]) return true;
			if (vecs[i][0] == element[1] && vecs[i][1] == element[0]) return true;
		}

		return false;
	}

	inline int VectorContainForSpecialCase2(std::vector<std::vector<int>> vecs,int element_0, int element_1)
	{
		for (int i = 0; i < vecs.size(); i++)
		{
			if (vecs[i][0] == element_0 && vecs[i][1] == element_1) return i;
			if (vecs[i][0] == element_1 && vecs[i][1] == element_0) return i;
		}
		return -1;
	}

	inline int VectorContainForSpecialCase3(std::vector<std::vector<int>> vecs, int element_0, int element_1)
	{
		for (int i = 0; i < vecs.size(); i++)
		{
			if (vecs[i][0] == element_0 && vecs[i][1] == element_1) return i;
		}
		return -1;
	}



	inline int VectorIndex(std::vector<int> vecs, int element)
	{
		for (int i = 0; i < vecs.size(); i++)
		{
			if (vecs[i] == element)
				return i;
		}

		return -1;
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

	inline int VectorIndex(std::vector<std::vector<int>> vecs, std::vector<int> element)
	{
		for (int i = 0; i < vecs.size(); i++)
		{
			if (vecs[i][0] == element[0] && vecs[i][1] == element[1]) return i;
			if (vecs[i][0] == element[1] && vecs[i][1] == element[0]) return i;
		}
		return -1;
	}

	inline Vector3d RelatedNormal(double angle_x, double angle_z)
	{
		Vector3d n(0.0, 1.0, 0.0);
		n= RotationAxis(n, angle_x, Vector3d(1.0, 0.0, 0.0));
		n = RotationAxis(n, angle_z, Vector3d(0.0, 0.0, 1.0));
		return n;
	}

	inline void RelatedSphereCoordinate(Vector3d n, double &angle_0, double &angle_1)
	{
		if (isAlmostZero(n[0]) && isAlmostZero(n[2]))
		{
			if (n[1] > 0)
			{
				angle_0 = 0.0;
				angle_1 = 0.0;
			}
			else
			{
				angle_0 = MM_PI;
				angle_1 = 0.0;
			}
			return;
		}
		angle_0 = getAngleBetween(n, Vector3d(0.0, 1.0, 0.0));
		angle_1 = getAngleBetween(Vector3d(0.0, 0.0, 1.0), Vector3d(n[0], 0.0, n[2]));
		if (n[0] <= 0)
		{
		//	angle_1 = 2 * MM_PI - angle_1;
		}
	}


	inline void OutputVectors(std::string out_path, const Vector3d3 &vecs)
	{
		std::ofstream file(out_path);
		file << vecs.size() << std::endl;

		for (int i = 0; i < vecs.size(); i++){
			file << vecs[i].size() << std::endl;
			for (int j = 0; j < vecs[i].size(); j++){
				file << vecs[i][j].size() << std::endl;
				for (int k = 0; k < vecs[i][j].size(); k++)
					file<< vecs[i][j][k][0] << " " << vecs[i][j][k][1] << " " << vecs[i][j][k][2] << " ";
				file << "" << std::endl;
			}
		}
	
		file.clear();
		file.close();
	}

	inline void OutputVectors(std::string out_path, const Vector3d1 &vecs)
	{
		std::ofstream file(out_path);
		file << vecs.size() << std::endl;
		for (int i = 0; i < vecs.size(); i++)
			file << vecs[i][0] << " " << vecs[i][1] << " " << vecs[i][2] <<std::endl;
		file.clear();
		file.close();
	}

	inline bool LoadExisting(std::string path)
	{
		std::ifstream file(path, std::ios::in);
		if (!file) return false;
		return true;
	}

	inline bool LoadVectors(std::string path, Vector3d3 &vec_3)
	{
		//zigzag_final_path
		int nb_0,nb_1,nb_2;
		std::ifstream file(path, std::ios::in);

		if (!file) return false;
		
		file >> nb_0;
		for (int i = 0; i < nb_0; i++)
		{
			file >> nb_1;
			Vector3d2 vec_2;
			for (int j = 0; j < nb_1; j++)
			{
				file >> nb_2;
				Vector3d1 vec_1(nb_2,Vector3d(0.0,0.0,0.0));
				for (int k = 0; k < nb_2; k++)
					file >> vec_1[k][0] >> vec_1[k][1] >> vec_1[k][2];
				vec_2.emplace_back(vec_1);
			}
			vec_3.emplace_back(vec_2);
		}
		file.clear();
		file.close();

		return true;
	}


	inline bool LoadVectors(std::string path, Vector3d1 &vec_3)
	{
		//zigzag_final_path
		std::ifstream file(path, std::ios::in);

		if (!file) return false;

		int nb;
		file >> nb;
		for (int i = 0; i < nb; i++)
		{
			vec_3.emplace_back(Vector3d());
			file >> vec_3.back()[0] >> vec_3.back()[1] >> vec_3.back()[2];
		}
		file.clear();
		file.close();

		return true;
	}




	inline void RelatedRotation(Vector3d n, double &angle_x, double &angle_z)
	{
		if (isAlmostZero(n[0]) && isAlmostZero(n[1]))
		{
			if (n[2] > 0)
			{
				angle_x = MM_PI / 2.0;
			}
			else
			{
				angle_x = -MM_PI / 2.0;
			}
			angle_z = 0.0;
			return;
		}

		angle_z = getAngleBetween(Vector3d(0.0, 1.0, 0.0), Vector3d(n[0], n[1], 0.0));

		if (n[0] > 0)
			angle_z = -angle_z;

		double l = getLength(n);
		double angle = getAngleBetween(n, Vector3d(0.0, 0.0, 1.0));

		Vector3d p(0.0, l*glm::sin(angle), l*glm::cos(angle));
		angle_x = getAngleBetween(Vector3d(0.0, 1.0, 0.0), p);

		if (n[2] < 0)
			angle_x = -angle_x;
	}

	inline void Output_tree(int nodes_nb, std::vector<int> &edges, std::string path)
	{
		std::ofstream file(path);

		file << "Mark Newman on Sat Jul 22 05:32:16 2006" << std::endl;
		file << "graph" << std::endl;
		file << "[" << std::endl;
		file << "  directed 0" << std::endl;

		for (int i = 0; i < nodes_nb; i++)
		{
			file << "node" << std::endl;
			file << "[" << std::endl;
			file << "id " << i << std::endl;
			file << "label " << i << std::endl;

			file << "]" << std::endl;
		}

		for (int i = 0; i < edges.size(); i = i + 2)
		{
			file << "edge" << std::endl;
			file << "[" << std::endl;

			file << "source " << edges[i] << std::endl;
			file << "target " << edges[i + 1] << std::endl;

			file << "]" << std::endl;
		}

		file << "]" << std::endl;

		file.clear();
		file.close();
	}
	
}

