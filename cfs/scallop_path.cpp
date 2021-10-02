#include "stdafx.h"
#include "CFSCNC.h"
#include <Circuit.h>


namespace cnc {


	Vector3d Binary_Search(Vector3d1& vecs, std::vector<int>& face_id_0, std::vector<int>& face_id_1, std::vector<int>& face_id_2,
		Vector3d center, Vector3d CC, Vector3d tangent, double min_angle, double max_angle, double distance)
	{
		Vector3d project = CGAL_3D_Projection_Point_Plane(CC, center, tangent);
		Vector3d base_p = Functs::SetVectorLength(project - center, 2.0) + center;
		double middle_angle = (min_angle + max_angle) / 2.0;
		Vector3d goal;
		while (true)
		{
			goal = Functs::RotationAxis(base_p - center, middle_angle, tangent) + center;
			double d = CGAL_3D_Distance_Point_Triangles(goal, vecs, face_id_0, face_id_1, face_id_2);
			if (abs(d - distance) < 0.001) break;
			if (d > distance)
			{
				max_angle = middle_angle;
				middle_angle = (min_angle + max_angle) / 2.0;
			}
			else
			{
				min_angle = middle_angle;
				middle_angle = (min_angle + max_angle) / 2.0;
			}
		}
		return goal;
	}
	void Construct_Local_Surface(Vector3d1& vecs, std::vector<int>& face_id_0, std::vector<int>& face_id_1, std::vector<int>& face_id_2, std::vector<int> triangle,
		Vector3d1 &local_vecs,std::vector<int> &local_face_id_0, std::vector<int> &local_face_id_1,std::vector<int> &local_face_id_2)
	{
		std::vector<int> m;
		for (int i = 0; i < triangle.size(); i++)
		{
			int id = triangle[i];
			int vertice_id_0 = VectorContainReturnIndex(m, face_id_0[id]);
			int vertice_id_1 = VectorContainReturnIndex(m, face_id_1[id]);
			int vertice_id_2 = VectorContainReturnIndex(m, face_id_2[id]);
			if (vertice_id_0 < 0){
				m.push_back(face_id_0[id]);
				vertice_id_0 = m.size() - 1;
			}

			if (vertice_id_1 < 0){
				m.push_back(face_id_1[id]);
				vertice_id_1 = m.size() - 1;
			}

			if (vertice_id_2 < 0){
				m.push_back(face_id_2[id]);
				vertice_id_2 = m.size() - 1;
			}

			local_face_id_0.push_back(vertice_id_0);
			local_face_id_1.push_back(vertice_id_1);
			local_face_id_2.push_back(vertice_id_2);
		}
		for (int i = 0; i < m.size(); i++)
			local_vecs.push_back(vecs[m[i]]);
	}

	void Trimming(Vector3d1 &first_path, Vector3d1& second_path, std::vector<bool> &debug0)
	{
		std::vector<bool> deleted_points(first_path.size(), false);

		for (int i = 0; i < first_path.size(); i++)
		{
			int index = i;
			int next_index = i;
			for (int j = 0; j < first_path.size(); j++)
			{
				if (!deleted_points[(i + j + 1) % first_path.size()])
				{
					next_index = (i + j + 1) % first_path.size();
					break;
				}
			}
			if (index == next_index) break;

			Vector3d s0 = second_path[index];
			Vector3d s1 = second_path[next_index];
			Vector3d cl0 = first_path[index];
			Vector3d cl1 = first_path[next_index];
			Vector3d n0 = Functs::GetCrossproduct(s1 - s0, cl0 - s0);
			Vector3d n1 = Functs::GetCrossproduct(s1 - cl1, cl0 - cl1);

			if (Math_PI - Functs::GetAngleBetween(n0, n1) > Math_PI / 2.0)
			{
				deleted_points[index] = true;
				deleted_points[next_index] = true;

				for (int j = 0; j < first_path.size(); j++)
				{
					if (!deleted_points[(i - j - 1 + first_path.size()) % first_path.size()])
					{
						next_index = (i - j - 1 + first_path.size()) % first_path.size();
						break;
					}
				}
				if (i != next_index)
					i = next_index;
			}
		}

		debug0 = deleted_points;

		Vector3d1 path;
		for (int i = 0; i < deleted_points.size(); i++)
			if (!deleted_points[i])
				path.push_back(second_path[i]);

		second_path = path;
	}


	void Searching_from_cl_to_scalllop(Vector3d1& vecs, std::vector<int>& face_id_0, std::vector<int>& face_id_1, std::vector<int>& face_id_2, 
		Vector3d1 &one_closed_cl_path, Vector3d1 &cc_path, Vector3d1 &one_closed_scallop_path, Vector3d1 &debug, std::vector<bool> &debug0 )
	{
		Vector3d1().swap(one_closed_scallop_path);

		Vector3d1 tangents;
		Circuit::Tangent_line(one_closed_cl_path, tangents);
		vector<std::vector<int>> triangles;

		CGAL_3D_Mesh_Near_Triangles(vecs, face_id_0, face_id_1, face_id_2, one_closed_cl_path, 2.0, triangles);

		for (int i = 0; i < one_closed_cl_path.size(); i++)
		{
			std::cout << i << " / " << one_closed_cl_path.size() << std::endl;
			int index = i;
			Vector3d1 local_vecs;
			std::vector<int> local_face_id_0;
			std::vector<int> local_face_id_1;
			std::vector<int> local_face_id_2;
			Construct_Local_Surface(vecs, face_id_0, face_id_1, face_id_2, triangles[index],
				local_vecs, local_face_id_0, local_face_id_1, local_face_id_2);

			Vector3d scallop_point = Binary_Search(local_vecs, local_face_id_0, local_face_id_1, local_face_id_2,
				one_closed_cl_path[index], cc_path[index], tangents[index], 0.0, Math_PI / 2.0, 0.2);

			one_closed_scallop_path.push_back(scallop_point);

			Vector3d1().swap(local_vecs);
			std::vector<int>().swap(local_face_id_0);
			std::vector<int>().swap(local_face_id_1);
			std::vector<int>().swap(local_face_id_2);

			debug.push_back(scallop_point);
			debug.push_back(one_closed_cl_path[i]);
		}
		/*****************************************************************************************/

		Trimming(one_closed_cl_path, one_closed_scallop_path, debug0);
		/*************************************************************************************************/
		Vector3d1().swap(tangents);
		vector<std::vector<int>>().swap(triangles);
	}

	void Searching_from_scallop_to_cl(Vector3d1& vecs, std::vector<int>& face_id_0, std::vector<int>& face_id_1, std::vector<int>& face_id_2,
		Vector3d1 &one_closed_scallop_path, Vector3d1 &cc_path, Vector3d1 &one_closed_cl_path)
	{
		Vector3d1().swap(one_closed_cl_path);

		Vector3d1 tangents;
		Circuit::Tangent_line(one_closed_scallop_path, tangents);

		vector<std::vector<int>> triangles;
		CGAL_3D_Mesh_Near_Triangles(vecs, face_id_0, face_id_1, face_id_2, one_closed_scallop_path, 2.0, triangles);

		for (int i = 0; i < one_closed_scallop_path.size(); i++)
		{
			std::cout << i << " / " << one_closed_scallop_path.size() << std::endl;
			int index = i;
			Vector3d1 local_vecs;
			std::vector<int> local_face_id_0;
			std::vector<int> local_face_id_1;
			std::vector<int> local_face_id_2;
			Construct_Local_Surface(vecs, face_id_0, face_id_1, face_id_2, triangles[index],
				local_vecs, local_face_id_0, local_face_id_1, local_face_id_2);

			Vector3d next_cl_point = Binary_Search(local_vecs, local_face_id_0, local_face_id_1, local_face_id_2,
				one_closed_scallop_path[index], cc_path[index], tangents[index], 0.0, Math_PI, 2.0);

			one_closed_cl_path.push_back(next_cl_point);
			
			Vector3d1().swap(local_vecs);
			std::vector<int>().swap(local_face_id_0);
			std::vector<int>().swap(local_face_id_1);
			std::vector<int>().swap(local_face_id_2);
		}
		//Trimming(one_closed_scallop_path, one_closed_cl_path);

		Vector3d1().swap(tangents);
		vector<std::vector<int>>().swap(triangles);
	}

	void CFSCNC::Extract_ISO_Scallop_Contours(std::string obj_path, std::string  full_obj_path, Vector3d2  &cutter_locations)
	{
		Vector3d1 vecs;
		std::vector<int> face_id_0;
		std::vector<int> face_id_1;
		std::vector<int> face_id_2;
		CGAL_3D_Read_Triangle_Mesh(obj_path, vecs, face_id_0, face_id_1, face_id_2);

#if 0

		//find cutter contact points
		/*******************************************************************************/
		Vector3d1 one_closed_cc_path;
		/*******************************************************************************/
		Vector3d2 boundarys;
		CGAL_3D_Triangle_Mesh_Boundary(vecs, face_id_0, face_id_1, face_id_2, boundarys);
		one_closed_cc_path = boundarys[0];

		//cutter_locations.push_back(one_closed_cc_path);

		//find cutter location points
		/*******************************************************************************/
		Vector3d1 one_closed_cl_path;
		/*******************************************************************************/
		Vector3d1 normals;
		CGAL_Normal_Mesh(full_obj_path, boundarys[0], normals);
		Circuit::Laplace_Smoothing(normals, 200);
		for (int i = 0; i < boundarys[0].size(); i++)
			boundarys[0][i] = boundarys[0][i] + Functs::SetVectorLength(normals[i], 2.0);
		one_closed_cl_path = boundarys[0];

		//cutter_locations.push_back(one_closed_cl_path);

		//find scallop curve points
		Vector3d1 one_closed_scallop_path;
		Searching_from_cl_to_scalllop(vecs, face_id_0, face_id_1, face_id_2, one_closed_cl_path, one_closed_cc_path, one_closed_scallop_path);
		cutter_locations.push_back(one_closed_scallop_path);


		for (int iter = 0; iter < 1; iter++)
		{
			Vector3d1 projects = CGAL_Project_Points_Onto_Surface(vecs, face_id_0, face_id_1, face_id_2, one_closed_scallop_path);
			Searching_from_scallop_to_cl(vecs, face_id_0, face_id_1, face_id_2, one_closed_scallop_path, projects, one_closed_cl_path);
			//cutter_locations.push_back(one_closed_cl_path);

			std::ofstream file("Z:\\Documents\\Windows\\SmartSFC\\workspace\\CFS\\save.obj");
			file.precision(8);
			file << one_closed_cl_path.size() << std::endl;

			for (int i = 0; i < one_closed_cl_path.size(); i++)
				file << one_closed_cl_path[i][0] << " " << one_closed_cl_path[i][1] << " " << one_closed_cl_path[i][2] << std::endl;

			file.clear();
			file.close();


			//projects = CGAL_Project_Points_Onto_Surface(vecs, face_id_0, face_id_1, face_id_2, one_closed_cl_path);
			//Searching_from_cl_to_scalllop(vecs, face_id_0, face_id_1, face_id_2, one_closed_cl_path, projects, one_closed_scallop_path);
			//cutter_locations.push_back(one_closed_scallop_path);

		}

#else
		Vector3d1 one_closed_cl_path;
		Vector3d1 one_closed_scallop_path;

		int nb;
		std::ifstream file("Z:\\Documents\\Windows\\SmartSFC\\workspace\\CFS\\save.obj", std::ios::in);

		file >> nb;
		for (int i = 0; i < nb; i++)
		{
			double x, y, z;
			file >> x >> y >> z;
			one_closed_cl_path.push_back(Vector3d(x, y, z));
		}


		file.clear();
		file.close();

		Vector3d1 projects = CGAL_Project_Points_Onto_Surface(vecs, face_id_0, face_id_1, face_id_2, one_closed_cl_path);
		Searching_from_cl_to_scalllop(vecs, face_id_0, face_id_1, face_id_2, one_closed_cl_path, projects, one_closed_scallop_path,debug,debug0);
		cutter_locations.push_back(one_closed_cl_path);
		cutter_locations.push_back(projects);
		cutter_locations.push_back(one_closed_scallop_path);

#endif // 0


	}
}

