#ifndef CIRCUIT_ONCE
#define CIRCUIT_ONCE
#pragma once

#include <stdafx.h>
#include <MathHelper.h>
#include <Strip.h>

#include <Constant.h>

using namespace hpcg;
using namespace std;

namespace cnc {

	class Circuit
	{
	public:

		//get center
		/***************************************************************************************************/
		static Vector3d GetCenter(Vector3d1 &input_points)
		{
			Vector3d center(0.0, 0.0, 0.0);

			for (int i = 0; i < input_points.size(); i++)
			{
				center[0] += input_points[i][0];
				center[1] += input_points[i][1];
				center[2] += input_points[i][2];
			}

			center[0] = center[0] / input_points.size();
			center[1] = center[1] / input_points.size();
			center[2] = center[2] / input_points.size();

			return center;
		}


		//get total distance
		/***************************************************************************************************/
		static double GetTotalLength(Vector3d1 &input_points)
		{
			double length = 0.0;
			for (int i = 0; i < input_points.size(); i++)
			{
				length += CGAL_3D_Distance_Point_Point(input_points[i], input_points[(i + 1) % input_points.size()]);
			}
			return length;
		}
		/***************************************************************************************************/

		//many kinds of distance computing in 2D
		/***************************************************************************************************/
		static double Distance(Vector2d v, Vector2d1 &input_points)
		{
			if (input_points.size() <= 1)
			{
				if (input_points.size() == 0)
				{
					return MAXDOUBLE;
				}
				else
				{
					return Strip::Distance(v, input_points[0]);
				}
			}
			else
			{
				double m_d = MAXDOUBLE;
				for (int i = 0; i < input_points.size(); i++)
				{
					double d = CGAL_2D_Distance_Point_Segment(v, input_points[i], input_points[(i + 1) % input_points.size()]);
					m_d = glm::min(m_d, d);
				}
				return m_d;
			}
		}

		//many kinds of distance computing in 3D
		/***************************************************************************************************/
		static double Distance(Vector3d v, Vector3d1 &input_points)
		{
			if (input_points.size() <= 1)
			{
				if (input_points.size() == 0)
				{
					return MAXDOUBLE;
				}
				else
				{
					return Strip::Distance(v, input_points[0]);
				}
			}
			else
			{
				double m_d = MAXDOUBLE;
				for (int i = 0; i < input_points.size(); i++)
				{
					double d = CGAL_3D_Distance_Point_Segment(v, input_points[i],input_points[(i + 1) % input_points.size()]);
					m_d = glm::min(m_d, d);
				}
				return m_d;
			}
		}

		static double GeodesicDistance(Vector3d v, Vector3d1 &input_points)
		{
			if (input_points.size() <= 1)
			{
				if (input_points.size() == 0)
				{
					return MAXDOUBLE;
				}
				else
				{
					return Strip::Distance(v, input_points[0]);
				}
			}
			else
			{
				double m_d = MAXDOUBLE;
				for (int i = 0; i < input_points.size(); i++)
				{
					double d = CGAL_Geodesic_Distance(INPUTPATH, v, input_points[i]);
					m_d = glm::min(m_d, d);
				}
				return m_d;
			}
		}

		static double Distance(Vector3d v, std::vector<Vector3d1 > &input_pointsese)
		{
			double min_d = MAXDOUBLE;
			for (int i = 0; i < input_pointsese.size(); i++)
			{
				double d = Distance(v, input_pointsese[i]);
				min_d = glm::min(min_d, d);
			}
			return min_d;
		}

		static double Distance(Vector3d v0, Vector3d v1, Vector3d1 &input_points)
		{
			double m_d = MAXDOUBLE;
			for (int i = 0; i < input_points.size(); i++)
			{
				double d = CGAL_3D_Distance_Segment_Segment(v0, v1, input_points[i],input_points[(i + 1) % input_points.size()]);
				m_d = glm::min(m_d, d);
			}
			return m_d;
		}

		static double Distance(Vector3d1 &input_points_0, Vector3d1 &input_points_1)
		{
			double m_d = MAXDOUBLE;

			for (int i = 0; i < input_points_0.size(); i++)
			{
				double d = Distance(input_points_0[i], input_points_0[(i + 1) % input_points_0.size()], input_points_1);
				m_d = glm::min(m_d, d);
			}
			for (int i = 0; i < input_points_1.size(); i++)
			{
				double d = Distance(input_points_1[i], input_points_1[(i + 1) % input_points_1.size()], input_points_0);
				m_d = glm::min(m_d, d);
			}
			return m_d;
		}
		/***************************************************************************************************/

		//GetOffsetPart
		/***************************************************************************************************/
		static double GetOffsetPartsLength(Vector3d2 &parts)
		{
			double length = 0.0;
			for (int i = 0; i < parts.size(); i++)
				length += Strip::GetTotalLength(parts[i]);
			return length;
		}
		/***************************************************************************************************/


		//parameter
		/***************************************************************************************************/
		static Vector3d FindNearestPoint(Vector3d v, Vector3d1 &input_points)
		{
			Vector3d n_p;

			double total_length = GetTotalLength(input_points);

			double min_d = MAXDOUBLE;
			int min_i = -1;
			for (int i = 0; i < input_points.size(); i++)
			{
				double l = CGAL_3D_Distance_Point_Segment(v, input_points[i], input_points[(i + 1) % input_points.size()]);

				if (l < min_d)
				{
					min_d = l;
					min_i = i;
				}
			}

			if (min_i >= 0)
			{
				Vector3d p0 = input_points[min_i];
				Vector3d p1 = input_points[(min_i + 1) % input_points.size()];

				double l0 = CGAL_3D_Distance_Point_Point(v, p0);
				double l1 = CGAL_3D_Distance_Point_Point(v, p1);

				if (min_d<l0 &&min_d<l1&&abs(min_d - l0)>0.0001&&abs(min_d - l1)>0.0001)
				{
					double length = std::sqrt(l0*l0 - min_d*min_d);
					double l = CGAL_3D_Distance_Point_Point(p0, p1);

					double w0 = (l - length) / l;
					double w1 = length / l;

					n_p[0] = w0*p0[0] + w1*p1[0];
					n_p[1] = w0*p0[1] + w1*p1[1];
					n_p[2] = w0*p0[2] + w1*p1[2];
				}
				else
				{
					if (l0 < l1)
						n_p = p0;
					else
						n_p = p1;
				}

			}
			else
			{
				assert(false);
			}

			return n_p;
		}

		static void FindNearestPoint(Vector3d v, Vector3d1 &input_points, int &location_index, double &par)
		{
			Vector3d n_p;

			double total_length = GetTotalLength(input_points);

			double min_d = MAXDOUBLE;
			int min_i = -1;
			for (int i = 0; i < input_points.size(); i++)
			{
				double l = CGAL_3D_Distance_Point_Segment(v, input_points[i], input_points[(i + 1) % input_points.size()]);

				if (l < min_d)
				{
					min_d = l;
					min_i = i;
				}
			}

			location_index = min_i;

			if (min_i >= 0)
			{
				Vector3d p0 = input_points[min_i];
				Vector3d p1 = input_points[(min_i + 1) % input_points.size()];

				double l0 = CGAL_3D_Distance_Point_Point(v, p0);
				double l1 = CGAL_3D_Distance_Point_Point(v, p1);

				if (min_d<l0 &&min_d<l1&&abs(min_d - l0)>0.0001&&abs(min_d - l1)>0.0001)
				{
					double length = std::sqrt(l0*l0 - min_d*min_d);
					double l = CGAL_3D_Distance_Point_Point(p0, p1);

					double w0 = (l - length) / l;
					double w1 = length / l;

					n_p[0] = w0*p0[0] + w1*p1[0];
					n_p[1] = w0*p0[1] + w1*p1[1];
					n_p[2] = w0*p0[2] + w1*p1[2];
					par = w1;
				}
				else
				{
					if (l0 < l1)
					{
						n_p = p0;
						par = 0.0;
					}
					else
					{
						n_p = p1;
						par = 1.0;
					}
				}
			}
			else
			{
				assert(false);
			}
		}

		//parameter
		/***************************************************************************************************/
		static double FindNearestPointPar(Vector3d v, Vector3d1 &input_points)
		{
			Vector3d n_p;

			double total_length = GetTotalLength(input_points);

			double min_d = MAXDOUBLE;
			int min_i = -1;
			for (int i = 0; i < input_points.size(); i++)
			{
				double l = CGAL_3D_Distance_Point_Segment(v,input_points[i], input_points[(i + 1) % input_points.size()]);

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
					length += CGAL_3D_Distance_Point_Point(input_points[i],input_points[(i + 1) % input_points.size()]);
				}

				Vector3d p0 = input_points[min_i];
				Vector3d p1 = input_points[(min_i + 1) % input_points.size()];

				double l0 = CGAL_3D_Distance_Point_Point(v, p0);
				double l1 = CGAL_3D_Distance_Point_Point(v, p1);

				if (min_d<l0 &&min_d<l1&&abs(min_d - l0)>0.0001&&abs(min_d - l1)>0.0001)
				{
					length += std::sqrt(l0*l0 - min_d*min_d);;
					return length / total_length;
				}
				else
				{
					if (l0 < l1)
						n_p= p0;
					else
						n_p = p1;

					length += CGAL_3D_Distance_Point_Point(input_points[min_i], n_p);
					return length / total_length;
				}

			}
			else
			{
				assert(false);
			}

			return -1.0;
		}

		static Vector3d GetOnePointFromOffset(double d, Vector3d1 &input_points)
		{
			Vector3d v;
			double length = 0.0;
			for (int i = 0; i < input_points.size(); i++)
			{
				length += CGAL_3D_Distance_Point_Point(input_points[i], input_points[(i + 1) % input_points.size()]);
			}

			double total_length = length;
			length = 0.0;

			for (int i = 0; i < input_points.size(); i++)
			{
				double l = CGAL_3D_Distance_Point_Point(input_points[i], input_points[(i + 1) % input_points.size()]);

				if (d*total_length >= length&&d*total_length <= length + l)
				{
					double ll = (d - length / total_length)*total_length / l;
					v[0] = input_points[i].x + (input_points[(i + 1) % input_points.size()].x - input_points[i].x)*ll;
					v[1] = input_points[i].y + (input_points[(i + 1) % input_points.size()].y - input_points[i].y)*ll;
					v[2] = input_points[i].z + (input_points[(i + 1) % input_points.size()].z - input_points[i].z)*ll;
					break;
				}
				length += l;
			}

			return v;
		}

		/***************************************************************************************************/
		static Vector3d1 SelectOnePartOffset(Vector3d1 &input_points, double d0, double d1)
		{
			Vector3d1 vecs;
			if (abs(d0 - d1) < 0.0000001)
			{
				Vector3d v = GetOnePointFromOffset(d0, input_points);
				vecs.push_back(v);
			}
			else
			{
				Vector3d v = GetOnePointFromOffset(d0, input_points);
				vecs.push_back(v);

				if (d1 >= 0 && d1 <= 1.0)
				{
					double total_length = GetTotalLength(input_points);

					std::vector<double> vec_ds;

					vec_ds.push_back(0.0);
					double length = 0.0;
					for (int i = 0; i < input_points.size(); i++)
					{
						length += CGAL_3D_Distance_Point_Point(input_points[i], input_points[(i + 1) % input_points.size()]);
						vec_ds.push_back(length / total_length);
					}


					if (d0 > d1)
					{
						for (int i = vec_ds.size() - 1; i >= 0; i--)
						{
							if (vec_ds[i]<d0&&vec_ds[i]>d1)
							{
								v = GetOnePointFromOffset(vec_ds[i], input_points);

								if (!(abs(v[0] - vecs[vecs.size() - 1][0]) < 0.000001&&abs(v[1] - vecs[vecs.size() - 1][1]) < 0.000001&&abs(v[2] - vecs[vecs.size() - 1][2]) < 0.000001))
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
								v = GetOnePointFromOffset(vec_ds[i], input_points);

								if (!(abs(v[0] - vecs[vecs.size() - 1][0]) < 0.000001&&abs(v[1] - vecs[vecs.size() - 1][1]) < 0.000001&&abs(v[2] - vecs[vecs.size() - 1][2]) < 0.000001))
								{
									vecs.push_back(v);
								}
							}
						}

						for (int i = vec_ds.size() - 1; i >0; i--)
						{
							if (vec_ds[i] > d1)
							{
								v = GetOnePointFromOffset(vec_ds[i], input_points);
								if (!(abs(v[0] - vecs[vecs.size() - 1][0]) < 0.000001&&abs(v[1] - vecs[vecs.size() - 1][1]) < 0.000001&&abs(v[2] - vecs[vecs.size() - 1][2]) < 0.000001))
								{
									vecs.push_back(v);
								}
							}
						}
					}

					v = GetOnePointFromOffset(d1, input_points);
					//if (!(abs(v[0] - vecs[vecs.size() - 1][0]) < 0.000001&&abs(v[1] - vecs[vecs.size() - 1][1]) < 0.000001&&abs(v[2] - vecs[vecs.size() - 1][2]) < 0.000001))
					{
						vecs.push_back(v);
					}

					if (abs(vecs[1][0] - vecs[0][0]) < 0.000001&&abs(vecs[1][1] - vecs[0][1]) < 0.000001&&abs(vecs[1][2] - vecs[0][2]) < 0.000001)
					{
						//vecs.erase(vecs.begin());
					}

					std::vector<double>().swap(vec_ds);
				}
			}

			return 	 vecs;
		}

		static void MakeIndexAsFirst(Vector3d1 &input_points,int index)
		{
			Vector3d1 temp;

			for (int i = index; i < input_points.size(); i++)
				temp.push_back(input_points[i]);
			for (int i = 0; i < index;i++)
				temp.push_back(input_points[i]);

			Vector3d1().swap(input_points);

			input_points.assign(temp.begin(), temp.end());
		}

		static void NearestPoint(Vector3d1 &input_points_0, Vector3d1 &input_points_1, int &nearest_int_0, int &nearest_int_1)
		{
			nearest_int_0 = -1;
			nearest_int_1 = -1;

			double m_d = MAXDOUBLE;
			for (int i = 0; i < input_points_0.size(); i++)
			{
				double d = Distance(input_points_0[i], input_points_1);
				if (d < m_d)
				{
					nearest_int_0 = i;
					m_d = d;
				}
			}

			m_d = MAXDOUBLE;
			for (int i = 0; i < input_points_1.size(); i++)
			{
				double d = Distance(input_points_1[i],  input_points_0);
				if (d < m_d)
				{
					nearest_int_1 = i;
					m_d = d;
				}
			}
		}

		static bool CheckSameDirection(Vector3d1 &input_points_0, Vector3d1 &input_points_1, double toolpath_size)
		{
			double delta = toolpath_size*TOLERANCE;
			
			int turn_vote = 0;
			int no_turn_vote = 0;

			for (int i = 0; i < input_points_0.size(); i++)
			{
				double d = Distance(input_points_0[i], input_points_1);
				if (abs(d - toolpath_size) < delta)
				{
					Vector3d upper_v_0 = input_points_0[i];
					double par_0 = Circuit::FindNearestPointPar(upper_v_0, input_points_0);
					double par_1 = Circuit::DeltaDEuclideanDistance(par_0, delta, input_points_0);
					Vector3d upper_v_1 = Circuit::GetOnePointFromOffset(par_1, input_points_0);

					par_0 = Circuit::FindNearestPointPar(upper_v_0, input_points_1);
					par_1 = Circuit::DeltaDEuclideanDistance(par_0, delta, input_points_1);
					Vector3d lower_v_0 = Circuit::GetOnePointFromOffset(par_0, input_points_1);
					Vector3d lower_v_1 = Circuit::GetOnePointFromOffset(par_1, input_points_1);
					double angle = getAngleBetween(upper_v_1 - upper_v_0, lower_v_1 - lower_v_0);
					if (angle > 3.141592653 / 2.0)
					{
						turn_vote++;
					}
					else
					{
						no_turn_vote++;
					}
				}
			}

			if (turn_vote + no_turn_vote > 0)
			{
				if (turn_vote > no_turn_vote)
					return false;
				else
					return true;
			}
			else
			{
				return true;
			}

		}

		static bool DistanceDouble(Vector3d1 &input_points_0, Vector3d1 &input_points_1, double toolpath_size)
		{
			bool b = false;

			for (int i = 0; i < input_points_0.size() && !b; i++)
			{
				double d = Distance(input_points_0[i], input_points_0[(i + 1) % input_points_0.size()], input_points_1);

				if (abs(d - toolpath_size) < toolpath_size*TOLERANCE)
				{
					b = true;
					return true;
				}
			}
			for (int i = 0; i < input_points_1.size() && !b; i++)
			{
				double d = Distance(input_points_1[i], input_points_1[(i + 1) % input_points_1.size()], input_points_0);

				if (abs(d - toolpath_size) < toolpath_size*TOLERANCE)
				{
					b = true;
					return true;
				}
			}
			return b;
		}

		//points moving on circuit
		static double DeltaDEuclideanDistance(double d, double distance, Vector3d1 &input_points)
		{
			double delta = distance / (GetTotalLength(input_points));
			double return_d=d+delta;

			if (return_d > 1.0)
				return_d = return_d - 1.0;

			if (return_d < 0.0)
				return_d = return_d + 1.0;

			return return_d;

		}

		//not nessary for now
		//points moving on circuit
		static Vector3d2 DeltaDEuclideanDistance(Vector3d v, double distance, Vector3d1 &input_points)
		{
			Vector3d2 output_vecs;

			if (GetTotalLength(input_points) <= 2 * distance)
				return output_vecs;

			///////////////////////////////////////////////////////////////////
			//int location_index;
			//double location_par;
			//FindNearestPoint(v, input_points, location_index, location_par);
			//par_0 = Circuit::FindNearestPointPar(parts_0[i][j], node_0->points);
			//par_1 = Circuit::DeltaDEuclideanDistance(par_0, toolpath_size, node_0->points);
			//cutting_point_0 = Circuit::GetOnePointFromOffset(par_0, node_0->points);
			//cutting_point_1 = Circuit::GetOnePointFromOffset(par_1, node_0->points);

			double par = FindNearestPointPar(v, input_points);
			double par_0 = DeltaDEuclideanDistance(par, distance, input_points);
			double par_1 = DeltaDEuclideanDistance(par, -distance, input_points);

			Vector3d1 vec_0 = Circuit::SelectOnePartOffset(input_points, par_0, par);
			std::reverse(vec_0.begin(),vec_0.end());
			output_vecs.push_back(vec_0);

			output_vecs.push_back(Circuit::SelectOnePartOffset(input_points, par, par_1));

			Strip::RemovingShortLines(output_vecs[0], 0.001);
			Strip::RemovingShortLines(output_vecs[1], 0.001);

			///////////////////////////////////////////////////////////////////
			//forward direction
			return output_vecs;
		}

		static void Laplace_Smoothing(Vector3d1 &input_points, int nb)
		{
			for (int j = 0; j < nb;j++)
			{
				Vector3d1 points;
				for (int i = 0; i < input_points.size(); i++)
				{
					Vector3d v0 = input_points[(i + input_points.size() - 1) % input_points.size()];
					Vector3d v1 = input_points[i];
					Vector3d v2 = input_points[(i + input_points.size() + 1) % input_points.size()];
					points.push_back((v0 + v1 + v2) / (float)3.0);
				}
				input_points = points;
			}
		}

		static void Laplace_Smoothing(Vector3d2 &input_points, int nb)
		{
			for (auto & ips : input_points)
			{
				Laplace_Smoothing(ips,nb);
			}
		}

		static void Tangent_line(Vector3d1 &input_points, Vector3d1 &tangents)
		{
			for (int i = 0; i < input_points.size(); i++)
			{
				Vector3d v0 = input_points[(i + input_points.size() - 1) % input_points.size()];
				Vector3d v1 = input_points[i];
				Vector3d v2 = input_points[(i + input_points.size() + 1) % input_points.size()];
			
				double length_0 = getLength(v1 - v0);
				double length_1 = getLength(v2 - v1);

				Vector3d tangent = ((float)length_0*SetVectorLength((v1 - v0), 1.0) + 
					                (float)length_1*SetVectorLength((v2 - v1), 1.0)) / (float)(length_0 + length_1);
				tangents.push_back(SetVectorLength(tangent, 1.0));
			}
		}


		static void KeepingOriginalSampling(Vector3d1 &input_points, double delta_length, Vector3d1 &output_points)
		{


			for (int i = 0; i < input_points.size(); i++)
			{
				auto s = input_points[i];
				auto e = input_points[(i + 1) % input_points.size()];
				output_points.emplace_back(s);

				double l = CGAL_3D_Distance_Point_Point(s,e);
				double d = delta_length;
				while (true)
				{
					if (d < l)
					{
						Vector3d v = (float)(1.0 - d / l)*s + (float)(d / l)*e;
						output_points.emplace_back(v);
						d = d + delta_length;
					}
					else
						break;
				}
			}
		}


		//remove duplicate points
		//remove collinear points
		static Vector2d1 ClearCircuit(const Vector2d1 &vecs)
		{
			//remove duplicate points
			Vector2d1 vecs_0;
			for (int i = 0; i < vecs.size(); i++)
			{
				if (i != vecs.size() - 1)
				{
					if (vecs_0.empty()) vecs_0.emplace_back(vecs[i]);
					else
						if (!Math::IsAlmostZero_Double(Math::GetLength(vecs_0.back(), vecs[i]),0.1))
							vecs_0.emplace_back(vecs[i]);
				}
				else
				{
					if (vecs_0.empty()) vecs_0.emplace_back(vecs[i]);
					else
					if (!Math::IsAlmostZero_Double(Math::GetLength(vecs_0.back(), vecs[i]), 0.1) &&
						!Math::IsAlmostZero_Double(Math::GetLength(vecs_0.front(), vecs[i]), 0.1))
							vecs_0.emplace_back(vecs[i]);
				}
			}

			if (Math::IsAlmostZero_Double(Math::GetLength(vecs_0.front(), vecs_0.back()),0.1)) vecs_0.erase(vecs_0.begin());

			//remove collinear points
			Vector2d1 vecs_1;
			for (int i = 0; i < vecs_0.size(); i++)
			{
				auto pre_v = vecs_0[(i + vecs_0.size() - 1) % vecs_0.size()];
				auto cur_v = vecs_0[i];
				auto next_v = vecs_0[(i + 1) % vecs_0.size()];
				double angle = Math::GetAngleBetween(cur_v - pre_v, next_v - cur_v);
				if (!Math::IsAlmostZero(angle))
					vecs_1.emplace_back(cur_v);
			}
			return vecs_1;
		}

		static void OutputPoints(const std::string path, const Vector3d1 &offset, const std::vector<std::string> names = std::vector<std::string>())
		{
			std::ofstream debug_file(path);
			int debug_index = 1;

			for (int j = 0; j < offset.size(); j++)
			{
				auto s = offset[j];
				if (names.size() == offset.size())
					CGAL_Export_Path_Point(debug_file, debug_index, names[j], 1.0, 0.0, 0.0, s, 0.1);
				else
					CGAL_Export_Path_Point(debug_file, debug_index, "name", 1.0, 0.0, 0.0, s, 0.1);
			}

			debug_file.clear();
			debug_file.close();
		}

		static void OutputPoints(const std::string path, const Vector3d2 &offset)
		{
			Vector3d1 points;
			for (const auto py : offset)for (const auto p : py)points.emplace_back(p);
			OutputPoints(path,points);
		}

		static void OutputOffsets(const std::string path, const Vector3d1 &offset)
		{
			std::ofstream debug_file(path);
			int debug_index = 1;

			for (int j = 0; j < offset.size(); j++)
			{
				auto s = offset[j];
				auto e = offset[(j + 1) % offset.size()];
				CGAL_Export_Path_Segment(debug_file, debug_index, "name", 1.0, 0.0, 0.0, s, e, 0.1);
			}

			debug_file.clear();
			debug_file.close();
		}

		static void OutputOffsets1(const std::string path, const Vector3d1 &offset)
		{
			std::ofstream debug_file(path);
			int debug_index = 1;

			for (int j = 0; j < offset.size(); j++)
			{
				auto s = offset[j];
				auto e = offset[(j + 1) % offset.size()];
				CGAL_Export_Path_Segment(debug_file, debug_index, "name", 1.0, 0.0, 0.0, s, e, 0.1);
			}

			debug_file.clear();
			debug_file.close();
		}

		static void OutputOffsets(const std::string path, const Vector3d2 &offsets, const std::string name)
		{
			std::ofstream debug_file(path);
			int debug_index = 1;
			for (int i = 0; i < offsets.size(); i++)
			{
				for (int j = 0; j < offsets[i].size(); j++)
				{
					auto s = offsets[i][j];
					auto e = offsets[i][(j + 1) % offsets[i].size()];
					CGAL_Export_Path_Segment(debug_file, debug_index, name+"_"+std::to_string(i), 1.0, 0.0, 0.0, s, e, 0.1);
				}
			}
			debug_file.clear();
			debug_file.close();
		}

		static void OutputSegments(const std::string path, const Vector3d2 &offsets, const std::string name)
		{
			std::ofstream debug_file(path);
			int debug_index = 1;
			for (int i = 0; i < offsets.size(); i++)
			{
				for (int j = 0; j < offsets[i].size()-1; j++)
				{
					auto s = offsets[i][j];
					auto e = offsets[i][(j + 1) % offsets[i].size()];
					CGAL_Export_Path_Segment(debug_file, debug_index, name+"_"+std::to_string(i), 1.0, 0.0, 0.0, s, e, 0.006);
				}
			}
			debug_file.clear();
			debug_file.close();
		}

		static void OutputOffsetsFiles(std::vector<string> &offsets_files, std::string path)
		{
			std::ofstream output_file(path);
			output_file << offsets_files.size() << std::endl;
			for (int i = 0; i < offsets_files.size(); i++)
				output_file << offsets_files[i] << std::endl;
			output_file.clear();
			output_file.close();
		}

		static void OutputStrips(const std::string path, const Vector3d2 &offsets, bool name_b=true)
		{
			std::ofstream debug_file(path);
			int debug_index = 1;
			for (int i = 0; i < offsets.size(); i++)
			{
				for (int j = 0; j < offsets[i].size() - 1; j++)
				{
					auto s = offsets[i][j];
					auto e = offsets[i][(j + 1) % offsets[i].size()];
					if (name_b)
						CGAL_Export_Path_Segment(debug_file, debug_index, "name_" + std::to_string(i), 1.0, 0.0, 0.0, s, e, 0.1);
					else
						CGAL_Export_Path_Segment(debug_file, debug_index, "name_" + std::to_string(i) + "_" + std::to_string(j), 1.0, 0.0, 0.0, s, e, 0.1);

				}
			}

			debug_file.clear();
			debug_file.close();
		}

		static void OutputStrip(const std::string path, const Vector3d1 &offsets, bool name_b=true)
		{
			std::ofstream debug_file(path);
			int debug_index = 1;
			for (int j = 0; j < offsets.size() - 1; j++)
			{
				auto s = offsets[j];
				auto e = offsets[(j + 1) % offsets.size()];
				if (name_b)
					CGAL_Export_Path_Segment(debug_file, debug_index, "name", 1.0, 0.0, 0.0, s, e, 0.2);
				else
					CGAL_Export_Path_Segment(debug_file, debug_index, "name_" + IntString(j), 1.0, 0.0, 0.0, s, e, 0.2);
			}

			debug_file.clear();
			debug_file.close();
		}
	};

}

#endif