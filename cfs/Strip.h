#ifndef STRIP_ONCE
#define STRIP_ONCE
#pragma once

#include <stdafx.h>
#include <MathHelper.h>

using namespace hpcg;
using namespace std;

namespace cnc {

	struct Sampling{
	public:
		double par;
		double angle;
		Vector3d v;
		double field;
		bool used;
		Sampling(double d0, double d1, Vector3d i)
		{
			par = d0;
			angle = d1;
			v = i;
			used = false;
		}
	};

	class Strip
	{
	public:
		
		//many kinds of distance computing in 2D
		/***************************************************************************************************/
		static double Distance(Vector2d v0, Vector2d v1)
		{
			return std::sqrt(pow(v0[0] - v1[0],2.0) + pow(v0[1] - v1[1],2.0));
		}

		//many kinds of distance computing in 3D
		/***************************************************************************************************/
		static double Distance(Vector3d v0, Vector3d v1)
		{
			return std::sqrt(pow(v0[0] - v1[0],2.0) + pow(v0[1] - v1[1],2.0) + pow(v0[2] - v1[2],2.0));
		}

		static double Distance(Vector3d v, Vector3d1 &input_points)
		{
			double m_d = MAXDOUBLE;
			for (int i = 0; i < input_points.size() - 1; i++)
			{
				double d = CGAL_3D_Distance_Point_Segment(v, input_points[i],input_points[(i + 1) % input_points.size()]);
				m_d = glm::min(m_d, d);
			}
			return m_d;
		}

		static double Distance(Vector3d v, std::vector< Vector3d1> &input_points)
		{
			double m_d = MAXDOUBLE;
			for (int i = 0; i < input_points.size(); i++)
			{
				m_d = glm::min(m_d, Distance(v, input_points[i]));
			}
			return m_d;
		}
		/***************************************************************************************************/

		//get total distance
		/***************************************************************************************************/
		static double GetTotalLength(Vector3d1 &input_points)
		{
			double length = 0.0;
			if (input_points.size() >= 2)
			{
				for (int i = 0; i < input_points.size() - 1; i++)
				{
					length += CGAL_3D_Distance_Point_Point(input_points[i], input_points[(i + 1) % input_points.size()]);
				}
			}
			return length;
		}

		//There maybe some problems in this function.
		static void Interpolation(Vector3d1 &input_points_0, Vector3d1 &input_points_1, Vector3d1 &result)
		{
			double total_length = GetTotalLength(input_points_0);

			Vector3d1 d_vector;

			for (int i = 0; i < input_points_0.size(); i++)
			{
				d_vector.push_back(input_points_0[i]);

				float d = GetTotalLength(d_vector) / total_length;

				Vector3d p_0 = input_points_0[i];
				Vector3d p_1 = Strip::GetOnePointFromStrip(d, input_points_1);

				Vector3d new_p = (float)(1.0 - d)*p_0 + d*p_1;
				result.push_back(new_p);
			}
			Vector3d1().swap(d_vector);
		}

		static Vector3d GetOnePointFromStrip(double d, Vector3d1 &input_points)
		{
			Vector3d v;

			double total_length = Strip::GetTotalLength(input_points);

			double length = 0.0;

			for (int i = 0; i < input_points.size() - 1; i++)
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

		static void UniformSampling(Vector3d1 &input_points, double delta_length, Vector3d1 &output_points)
		{
			double total_length = Strip::GetTotalLength(input_points);
			double length = 0.0;
			for (int i = 0; i < input_points.size() - 1; i++)
			{
				double l = CGAL_3D_Distance_Point_Point(input_points[i], input_points[i + 1]);
				double d = ((int)(length / delta_length) + 1)*delta_length-length;
				while (true){
					if (d >= 0 && d < l){
						double ll = d / l;
						Vector3d v = (float)(1.0 - ll)*input_points[i] + (float)(ll)*input_points[i + 1];
						output_points.push_back(v);
						d += delta_length;
					}
					else{
						break;
					}
				}
				length += l;
			}
		}


		static void UniformSampling(Vector3d1 &input_points, double delta_length, Vector3d1 &output_points, std::vector<bool> &lables)
		{
			std::vector<bool> new_lables;
			double total_length = Strip::GetTotalLength(input_points);
			double length = 0.0;
			for (int i = 0; i < input_points.size() - 1; i++)
			{
				double l = CGAL_3D_Distance_Point_Point(input_points[i], input_points[i + 1]);
				double d = ((int)(length / delta_length) + 1)*delta_length - length;
				while (true){
					if (d >= 0 && d < l){
						double ll = d / l;
						Vector3d v = (float)(1.0 - ll)*input_points[i] + (float)(ll)*input_points[i + 1];
						output_points.push_back(v);
						new_lables.push_back(lables[i]);
						d += delta_length;
					}
					else{
						break;
					}
				}
				length += l;
			}
			lables = new_lables;
		}

		static bool compare_double(const double &d0, const double &d1)
		{
			return d0 > d1;
		}
		static void AdaptiveSampling(Vector3d1 &input_path, double chord_error)
		{
			/*
			double total_length = Strip::GetTotalLength(input_path);
			double delta = 0.01;
			int sampling_points_nb = total_length / delta;

			//uniform sampling
			Vector3d1 temp;
			for (int i = 0; i < input_path.size(); i++)
				temp.push_back(input_path[i]);
			Vector3d1().swap(input_path);

			for (int i = 0; i <= sampling_points_nb; i++)
			{
				double par = i*1.0 / (double)sampling_points_nb;
				Vector3d v = Strip::GetOnePointFromStrip(par, temp);
				input_path.push_back(v);
			}
			*/

			//selectiong based on chord_error
			Vector3d1 temp = input_path;
			Vector3d1().swap(input_path);

			input_path.push_back(temp[0]);

			int start_index = 0;
			for (int i = 1; i < temp.size()-1; i++)
			{
				int current_index = i;
				for (int j = current_index - 1; j >= start_index + 1; j--)
				{
					double d = CGAL_3D_Distance_Point_Segment(temp[j], temp[start_index], temp[current_index]); 
					if (d > chord_error)
					{
						input_path.push_back(temp[i]);
						start_index = i;
						break;
					}
				}
			}
			input_path.push_back(temp[temp.size()-1]);
		}

		//static Adaptive
		static void AdaptiveSampling(double toolpath_size, Vector3d1 &input_path, int sampling_points_nb, int select_points_nb)
		{
			double total_length = GetTotalLength(input_path);

			//uniform sampling
			std::vector<Sampling> sampling_points;
			if (sampling_points_nb < 3)
			{
				return;
			}

			if (sampling_points_nb < select_points_nb)
			{
				Vector3d1 uniform_sampling;
				UniformSampling(input_path, sampling_points_nb, uniform_sampling);
				Vector3d1().swap(input_path);
				for (int i = 0; i < uniform_sampling.size(); i++){
					input_path.push_back(uniform_sampling[i]);
				}
				Vector3d1().swap(uniform_sampling);
				return;
			}

			//uniform sampling
			std::vector<double> save_smooth;
			Vector3d1 uniform_sampling;
			UniformSampling(input_path, sampling_points_nb, uniform_sampling);
			for (int i = 0; i < uniform_sampling.size(); i++)
			{
				double par = i*1.0 / (double)sampling_points_nb;
				sampling_points.push_back(Sampling(par, 0.0, uniform_sampling[i]));
				save_smooth.push_back(0.0);
			}

			for (int i = 1; i < sampling_points.size() - 1; i++)
			{
				Vector3d v0 = sampling_points[i - 1].v;
				Vector3d v1 = sampling_points[i].v;
				Vector3d v2 = sampling_points[i + 1].v;

				double angle = getAngleBetween(v0-v1,v2-v1);
				sampling_points[i].angle = MM_PI - angle;
			}
			sampling_points[0].angle = sampling_points[1].angle;
			sampling_points[sampling_points.size() - 1].angle = sampling_points[sampling_points.size() - 2].angle;


			//updating angle

			int iter_nb = 2.0*toolpath_size / (total_length / (double)sampling_points_nb);

			for (int iter = 0; iter < iter_nb; iter++)
			{
				for (int i = 0; i < sampling_points.size(); i++)
					save_smooth[i] = sampling_points[i].angle;

				for (int i = 1; i < sampling_points.size() - 1; i++)
				{
					double angle_0 = sampling_points[i - 1].angle;
					double angle_1 = sampling_points[i].angle;
					double angle_2 = sampling_points[i + 1].angle;
					save_smooth[i] = (angle_0 + angle_1 + angle_2) / 3.0;
				}

				for (int i = 0; i < save_smooth.size(); i++)
					sampling_points[i].angle = save_smooth[i];

				std::cout << iter << " " << sampling_points[1].angle << std::endl;
			}

			double min_angle = MAXDOUBLE;
			double max_angle = -MAXDOUBLE;

			for (int i = 1; i < sampling_points.size() - 1; i++)
			{
				min_angle = min(sampling_points[i].angle, min_angle);
				max_angle = max(sampling_points[i].angle, max_angle);
			}

			for (int i = 1; i < sampling_points.size() - 1; i++)
			{
				sampling_points[i].field = (sampling_points[i].angle - min_angle) / (max_angle - min_angle);
			}

			//select
			int select_nb = 0;

			std::cout << "Sampling points" << std::endl;
			int save_duration0 = -1;
			while (select_nb<select_points_nb)
			{
				double random_0 = rand() / double(RAND_MAX);
				int index = (int)(random_0*(double)sampling_points_nb);

				bool b = false;
				if (index > 0 && index < sampling_points.size() - 1)
				{
					double random_1 = rand() / double(RAND_MAX);
					if (!sampling_points[index].used&&random_1 < sampling_points[index].field)
					{
						select_nb++;
						sampling_points[index].used = true;
						b = true;
					}
				}
			}
			sampling_points[0].used = true;
			sampling_points[sampling_points.size()-1].used = true;

			//update points
			Vector3d1().swap(input_path);

			for (int i = 0; i < sampling_points.size(); i++)
			{
				if (sampling_points[i].used)
					input_path.push_back(sampling_points[i].v);
			}
		}

		static void RemovingShortLines(Vector3d1 &input_points, double d)
		{
			Vector3d1 temp = input_points;
			Vector3d1().swap(input_points);

			input_points.push_back(temp[0]);
			for (int i = 1; i < temp.size(); i++)
			{
				if (CGAL_3D_Distance_Point_Point(input_points[input_points.size() - 1], temp[i])>d)
					input_points.push_back(temp[i]);
			}
			Vector3d1().swap(temp);
		}

		static void SmoothingLines(Vector3d1 &input_points, int iter)
		{
			for (int j = 0; j < iter; j++)
			{
				Vector3d1 temp;
				temp.push_back(input_points[0]);
				for (int i = 1; i < input_points.size()-1; i++)
				{
					Vector3d v0 = input_points[i - 1];
					Vector3d v1 = input_points[i];
					Vector3d v2 = input_points[i + 1];
					v1 = (v0 + v1 + v2)/(float)3.0;
					temp.push_back(v1);
				}

				temp.push_back(input_points[input_points.size() - 1]);
				input_points = temp;
				Vector3d1().swap(temp);
			}
		}

		static void SmoothingLines(Vector3d1 &input_points, double shortest_line, double chord_error, double smallest_angle)
		{
			Vector3d1 original_input_points = input_points;

			while (true)
			{
				int nb0 = 0;
				int nb = 0;
				Vector3d1 temp;
				temp.push_back(input_points[0]);
				for (int i = 1; i < input_points.size() - 1; i++)
				{
					Vector3d v0 = input_points[i - 1];
					Vector3d v1 = input_points[i];
					Vector3d v2 = input_points[i + 1];
					double angle = getAngleBetween(v1-v0,v2-v1);
					angle = angle / MM_PI*180.0;
					Vector3d v = (v0 + v1 + v2) / (float)3.0;
					double d = getLength(v - original_input_points[i]);

					if (angle >= smallest_angle&&d <= chord_error)
					{
						nb++;
						v1 = v;
					}

					if (angle >= smallest_angle)
					{
						nb0++;
					}

					temp.push_back(v1);
				}
				temp.push_back(input_points[input_points.size() - 1]);
				input_points = temp;
				Vector3d1().swap(temp);

				std::cout << nb0 << std::endl;

				if (nb == 0) break;
			}
		}

		static double GetTotalAngles(Vector3d1 &input_points)
		{
			double total_angle = 0.0;

			for (int i = 1; i < input_points.size()-1; i++)
			{
				Vector3d v0 = input_points[i - 1];
				Vector3d v1 = input_points[i];
				Vector3d v2 = input_points[i+1];
				double angle = getAngleBetween(v0-v1,v2-v1);
				total_angle += angle;
			}
			return total_angle;
		}

		static void ConnectingSameDirectionLines(Vector3d1 &input_points)
		{
			Vector3d1 temp = input_points;
			Vector3d1().swap(input_points);

			input_points.push_back(temp[0]);
			input_points.push_back(temp[1]);


			for (int i = 2; i < temp.size(); i++)
			{
				Vector3d v0 = input_points[input_points.size() - 2];
				Vector3d v1 = input_points[input_points.size() - 1];
				Vector3d v2 = temp[i];
				double angle = getAngleBetween(v1-v0,v2-v1);
				angle = angle / MM_PI*180.0;

				if (isAlmostZero(angle))
				{
					input_points.erase(input_points.begin()+input_points.size()-1);
				}

				input_points.push_back(temp[i]);
			}
			Vector3d1().swap(temp);
		}
	};
}
#endif