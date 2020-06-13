#include "stdafx.h"
#include "cgalpackage.h"
#include <iostream>
#include <fstream>
#include <cassert>

#include <vector>
#include <string>
#include <algorithm>
#include <list>

#define _USE_MATH_DEFINES  
#include <math.h>

double Strip_Get_Total_length(Vector2d1 &input_points)
{
	double total_length = 0.0;
	if (input_points.size() >= 2){
		for (int i = 0; i < input_points.size() - 1; i++){
			total_length += CGAL_2D_Distance_Point_Point(input_points[i], input_points[i + 1]);
		}
	}
	return total_length;
}

double Circuit_Get_Total_length(Vector2d1 &input_points)
{
	double total_length = 0.0;
	for (int i = 0; i < input_points.size(); i++){
		total_length += CGAL_2D_Distance_Point_Point(input_points[i], input_points[(i + 1) % input_points.size()]);
	}
	return total_length;
}

Vector2d Strip_Get_One_Point_From_Strip(double d, Vector2d1 &input_points, int &index)
{
	Vector2d v;

	double length = 0.0;
	for (int i = 0; i < input_points.size() - 1; i++)
	{
		length += sqrt((double)CGAL_2D_Distance_Point_Point(input_points[i], input_points[(i + 1) % input_points.size()]));
	}
	double total_length = length;

	length = 0.0;

	for (int i = 0; i < input_points.size() - 1; i++)
	{
		double l = sqrt((double)CGAL_2D_Distance_Point_Point(input_points[i], input_points[(i + 1) % input_points.size()]));

		if (d*total_length >= length&&d*total_length <= length + l)
		{
			double ll = (d - length / total_length)*total_length / l;
			v[0] = input_points[i].x + (input_points[(i + 1) % input_points.size()].x - input_points[i].x)*ll;
			v[1] = input_points[i].y + (input_points[(i + 1) % input_points.size()].y - input_points[i].y)*ll;
			index = i;
			break;
		}
		length += l;
	}
	return v;
}

void GenerateACircle(int divided_nb, Vector2d v, double distance, Vector2d1 &circle_points)
{
	for (int i = 0; i < divided_nb; i++)
	{
		Vector2d p0(v[0] + abs(distance)*sin(i * 2 * M_PI / (double)divided_nb), v[1] + abs(distance)*cos(i * 2 * M_PI / (double)divided_nb));
		circle_points.push_back(p0);
	}
}

Vector2d CircleIntersectWithSegment(Vector2d1 &circle_points, Vector2d v0, Vector2d v1)
{
	for (int i = 0; i < circle_points.size(); i++)
	{
		Vector2d inter;
		if (CGAL_2D_Intersection_Segment_Segment(circle_points[i], circle_points[(i + 1) % circle_points.size()], v0, v1, inter))
		{
			return inter;
		}
	}
	Vector2d v(0.0, 0.0);
	return v;
}

void FindLastIndexPoint(Vector2d1 &input_points, Vector2d center, double radius, Vector2d1 &circle_points, int start_index, Vector2d1 &output_points)
{
	for (int i = start_index; i >= 0; i--)
	{
		if (i >= 0 && i < input_points.size())
		{
			double distance = CGAL_2D_Distance_Point_Point(input_points[i], center);

			if (CGAL_2D_Location_Point_Polygon(input_points[i], circle_points))
			{
				output_points.push_back(input_points[i]);
			}
			else
			{
				Vector2d v;
				if (output_points.size() > 0)
				{
					v = CircleIntersectWithSegment(circle_points, input_points[i], output_points[output_points.size() - 1]);
				}
				else
				{
					v = CircleIntersectWithSegment(circle_points, input_points[i], center);
				}
				output_points.push_back(v);
				break;
			}
		}
	}
}


void FindNextIndexPoint(Vector2d1 &input_points, Vector2d center, double radius, Vector2d1 &circle_points, int start_index, Vector2d1 &output_points)
{
	for (int i = start_index; i < input_points.size(); i++)
	{
		if (i >= 0 && i < input_points.size())
		{
			if (CGAL_2D_Location_Point_Polygon(input_points[i], circle_points))
			{
				output_points.push_back(input_points[i]);
			}
			else
			{
				Vector2d v;

				if (output_points.size() > 0)
				{
					v = CircleIntersectWithSegment(circle_points, input_points[i], output_points[output_points.size() - 1]);
				}
				else
				{
					v = CircleIntersectWithSegment(circle_points, input_points[i], center);
				}
				output_points.push_back(v);
				break;
			}
		}
	}
}

//par system
double Circuit_Find_Nearest_Point_Par(Vector2d v, Vector2d1 &input_points)
{
	Vector2d n_p;

	double total_length = Circuit_Get_Total_length(input_points);

	double min_d = DBL_MAX;
	int min_i = -1;
	for (int i = 0; i < input_points.size(); i++)
	{
		Vector2d p0(input_points[i].x, input_points[i].y);
		Vector2d p1(input_points[(i + 1) % input_points.size()].x, input_points[(i + 1) % input_points.size()].y);

		double l = CGAL_2D_Distance_Point_Segment(v,p0,p1);

		if (l < min_d)
		{
			min_d = l;
			min_i = i;
		}
	}

	if (min_i >= 0)
	{
		double length = 0.0;
		for (int i = 0; i < min_i; i++)
		{
			length += CGAL_2D_Distance_Point_Point(input_points[i], input_points[(i + 1) % input_points.size()]);
		}

		Vector2d p0(input_points[min_i].x, input_points[min_i].y);
		Vector2d p1(input_points[(min_i + 1) % input_points.size()].x, input_points[(min_i + 1) % input_points.size()].y);

		double l0 = CGAL_2D_Distance_Point_Point(v, p0);
		double l1 = CGAL_2D_Distance_Point_Point(v, p1);

		if (min_d<l0 &&min_d<l1&&abs(min_d - l0)>0.0001&&abs(min_d - l1)>0.0001)
		{
			double l = CGAL_2D_Distance_Point_Point(p0, p1);
			if (l < 0.00001)
			{
				v[0] = p0[0];
				v[1] = p0[1];
			}
			else
			{
				Vector2d vec(p1[0] - p0[0], p1[1] - p0[1]);
				Vector2d r_vec(vec[1], -vec[0]);
				if (vec[0] < 0.00001)
				{
					r_vec[0] = -vec[1];
					r_vec[1] = vec[0];
				}

				if (vec[1] < 0.00001)
				{
					r_vec[0] = vec[1];
					r_vec[1] = -vec[0];
				}

				if (CGAL_2D_Intersection_Segment_Line(p0, p1, v, v + r_vec, n_p))
				{

				}
				else
				{
					assert(false);
				}
			}
		}
		else
		{
			if (l0 < l1)
			{
				n_p[0] = p0[0];
				n_p[1] = p0[1];
			}
			else
			{
				n_p[0] = p1[0];
				n_p[1] = p1[1];
			}

		}

		length += CGAL_2D_Distance_Point_Point(input_points[min_i], n_p);

		return length / total_length;
	}
	else
	{
		assert(false);
	}

	return -1.0;
}


Vector2d Circuit_Get_One_Point_From_Offset(double d, Vector2d1 &input_points)
{
	Vector2d v;
	double length = 0.0;
	for (int i = 0; i < input_points.size(); i++)
	{
		length += CGAL_2D_Distance_Point_Point(input_points[i], input_points[(i + 1) % input_points.size()]);
	}

	double total_length = length;
	length = 0.0;

	for (int i = 0; i < input_points.size(); i++)
	{
		double l = CGAL_2D_Distance_Point_Point(input_points[i], input_points[(i + 1) % input_points.size()]);

		if (d*total_length >= length&&d*total_length <= length + l)
		{
			double ll = (d - length / total_length)*total_length / l;
			v[0] = input_points[i].x + (input_points[(i + 1) % input_points.size()].x - input_points[i].x)*ll;
			v[1] = input_points[i].y + (input_points[(i + 1) % input_points.size()].y - input_points[i].y)*ll;
			break;
		}
		length += l;
	}

	return v;
}


void Circuit_Select_One_Part_Offset(Vector2d1 &input_points, double d0, double d1, Vector2d1 &vecs)
{
	if (abs(d0 - d1) < 0.0000001)
	{
		Vector2d v = Circuit_Get_One_Point_From_Offset(d0, input_points);
		vecs.push_back(v);
	}
	else
	{
		Vector2d v = Circuit_Get_One_Point_From_Offset(d0, input_points);
		vecs.push_back(v);

		if (d1 >= 0 && d1 <= 1.0)
		{
			double total_length = Circuit_Get_Total_length(input_points);

			std::vector<double> vec_ds;

			vec_ds.push_back(0.0);
			double length = 0.0;
			for (int i = 0; i < input_points.size(); i++)
			{
				length += CGAL_2D_Distance_Point_Point(input_points[i], input_points[(i + 1) % input_points.size()]);
				vec_ds.push_back(length / total_length);
			}


			if (d0 > d1)
			{
				for (int i = vec_ds.size() - 1; i >= 0; i--)
				{
					if (vec_ds[i]<d0&&vec_ds[i]>d1)
					{
						v = Circuit_Get_One_Point_From_Offset(vec_ds[i], input_points);

						if (!(abs(v[0] - vecs[vecs.size() - 1][0]) < 0.000001&&abs(v[1] - vecs[vecs.size() - 1][1]) < 0.000001))
						{
							vecs.push_back(v);
						}
					}

					if (vec_ds[i] < d1)
					{
						break;
					}
				}
			}
			else
			{
				for (int i = vec_ds.size() - 1; i >0; i--)
				{
					if (vec_ds[i] < d0)
					{
						v = Circuit_Get_One_Point_From_Offset(vec_ds[i], input_points);

						if (!(abs(v[0] - vecs[vecs.size() - 1][0]) < 0.000001&&abs(v[1] - vecs[vecs.size() - 1][1]) < 0.000001))
						{
							vecs.push_back(v);
						}
					}
				}

				for (int i = vec_ds.size() - 1; i >0; i--)
				{
					if (vec_ds[i] > d1)
					{
						v = Circuit_Get_One_Point_From_Offset(vec_ds[i], input_points);
						if (!(abs(v[0] - vecs[vecs.size() - 1][0]) < 0.000001&&abs(v[1] - vecs[vecs.size() - 1][1]) < 0.000001))
						{
							vecs.push_back(v);
						}
					}
				}
			}

			v = Circuit_Get_One_Point_From_Offset(d1, input_points);
			if (!(abs(v[0] - vecs[vecs.size() - 1][0]) < 0.000001&&abs(v[1] - vecs[vecs.size() - 1][1]) < 0.000001))
			{
				vecs.push_back(v);
			}

			if (abs(vecs[1][0] - vecs[0][0]) < 0.000001&&abs(vecs[1][1] - vecs[0][1]) < 0.000001)
			{
				vecs.erase(vecs.begin());
			}

			std::vector<double>().swap(vec_ds);
		}
	}
}


double GetHalfCircleArea(Vector2d1 &last_points, Vector2d1 &next_points, Vector2d center, double radius, Vector2d1 &circle_points)
{
	Vector2d1 polygon_points;
	polygon_points.push_back(center);

	for (int i = 0; i < last_points.size() - 1; i++)
		polygon_points.push_back(last_points[i]);

	double par_0 = Circuit_Find_Nearest_Point_Par(last_points[last_points.size() - 1], circle_points);
	double par_1 = Circuit_Find_Nearest_Point_Par(next_points[next_points.size() - 1], circle_points);
	Circuit_Select_One_Part_Offset(circle_points, par_0, par_1, polygon_points);

	for (int i = next_points.size() - 2; i >= 0; i--)
		polygon_points.push_back(next_points[i]);

	return CGAL_2D_Polygon_Area(polygon_points);
}


void CGAL_Intergral_Curvature(Vector2d1 &input_points, int sampling_points_nb, double radius, double thresholder,
	Vector2d1 &output_points, std::vector<double> &output_rates)
{
	double total_length = Strip_Get_Total_length(input_points);

	for (int i = 0; i < sampling_points_nb; i++)
	{
		if (i % (int)(sampling_points_nb/100.0) == 0)
		{
			std::cout << "Intergral: " << i << " / " << sampling_points_nb << std::endl;
		}

		double par = i*1.0 / (double)sampling_points_nb;

		if (par > radius / total_length && par<1.0 - radius / total_length)
		{

			int last_index;
			Vector2d center = Strip_Get_One_Point_From_Strip(par, input_points, last_index);

			Vector2d1 circle_points;
			GenerateACircle(30, center, radius, circle_points);

			Vector2d1 last_points;
			FindLastIndexPoint(input_points, center, radius, circle_points, last_index, last_points);
			Vector2d1 next_points;
			FindNextIndexPoint(input_points, center, radius, circle_points, last_index + 1, next_points);

			if (last_points.size() > 0 && next_points.size() > 0)
			{
				double area = GetHalfCircleArea(last_points, next_points, center, radius, circle_points);

				if (area > M_PI*radius*radius*0.5)
				{
					area = M_PI*radius*radius - area;
				}

				area = area / (M_PI*radius*radius);

				if (area < thresholder)
				{
					output_points.push_back(center);
					output_rates.push_back(area);
				}
			}
		}

	}
}