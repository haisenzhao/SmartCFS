#include "stdafx.h"
#include "CFSCNC.h"
#include <Circuit.h>

namespace cnc {



	void CFSCNC::Output_Path(std::string path)
	{
		//Vector3d1 normals;
		//CGAL_Normal_Mesh(input_path + "_"+ "_one_region.off", single_path, normals);

		//samplng
		double gap = ComputeGapFromScallop(0.0, 2.0, max_scallop);

		Vector3d1 output_points;
		Strip::UniformSampling(single_path, gap*0.1, output_points, single_path_fixed_label);	
		single_path = output_points;


		//output
		std::ofstream file(path);
		file.precision(8);

		file << max_scallop << std::endl;
		file << drill_radius << std::endl;
		file << toolpath_size << std::endl;
		file << single_path.size() << std::endl;

		//fix fixed_label bug

		if (single_path_fixed_label.size() != 0)
		{
			double d = 0.0;

			int index = 0;
			while (d < gap*2.0)
			{

				d = d + CGAL_3D_Distance_Point_Point(single_path[index], single_path[index+1]);
				single_path_fixed_label[index] = 1;
				index++;
			}
		}


		for (int i = 0; i < single_path.size(); i++)
		{
			if (single_path_fixed_label.size() != 0)
			{
				file << single_path[i][0] << " " << single_path[i][1] << " " << single_path[i][2] << " " 
					 << " " << single_path_fixed_label[i] << std::endl;
			}
			else
				file << single_path[i][0] << " " << single_path[i][1] << " " << single_path[i][2] << std::endl;
		}

		file.clear();
		file.close();
	}

	void CFSCNC::Output_Heat_Loop_Curvature_Factor(std::string in_path, std::string out_path)
	{
		Vector3d1 heat_surface_vertices;
		std::vector<int> heat_face_id_0;
		std::vector<int> heat_face_id_1;
		std::vector<int> heat_face_id_2;
		CGAL_3D_Read_Triangle_Mesh(in_path, heat_surface_vertices, heat_face_id_0, heat_face_id_1, heat_face_id_2);

		std::vector<double> base_vertices_max_cur;
		std::vector<double> base_vertices_min_cur;
		Vector3d1 base_vertices_max_cur_direction;
		Vector3d1 base_vertices_min_cur_direction;

		CGAL_Curvature_Mesh(input_path + "\\" + IntString(cfs_index) + "_full.off",
			heat_surface_vertices, base_vertices_max_cur, base_vertices_min_cur, base_vertices_max_cur_direction, base_vertices_min_cur_direction);

		std::ofstream file(out_path);
		file << heat_surface_vertices.size() << std::endl;

		for (int i = 0; i < heat_surface_vertices.size(); i++)
		{
			file << base_vertices_max_cur[i] << " " << base_vertices_min_cur[i] << " "
				<< base_vertices_max_cur_direction[i][0] << " " << base_vertices_max_cur_direction[i][1] << " " << base_vertices_max_cur_direction[i][2] << " "
				<< base_vertices_min_cur_direction[i][0] << " " << base_vertices_min_cur_direction[i][1] << " " << base_vertices_min_cur_direction[i][2] << std::endl;
		}

		file.clear();
		file.close();

		Vector3d1().swap(heat_surface_vertices);
		std::vector<int>().swap(heat_face_id_0);
		std::vector<int>().swap(heat_face_id_1);
		std::vector<int>().swap(heat_face_id_2);
		std::vector<double>().swap(base_vertices_max_cur);
		std::vector<double>().swap(base_vertices_min_cur);
	}

	void CFSCNC::Output_Boundary(std::string in_path, std::string out_path)
	{
		std::cout << "Output boundary source for heat geodesic computing" << std::endl;
		Vector3d1 heat_surface_vertices;
		std::vector<int> heat_face_id_0;
		std::vector<int> heat_face_id_1;
		std::vector<int> heat_face_id_2;

		CGAL_3D_Read_Triangle_Mesh(in_path, heat_surface_vertices, heat_face_id_0, heat_face_id_1, heat_face_id_2);
		std::vector<bool> heat_vertices_boundary;
		CGAL_3D_Triangle_Mesh_Boundary(heat_surface_vertices, heat_face_id_0, heat_face_id_1, heat_face_id_2, heat_vertices_boundary);

		int nb = 0;
		for (int i = 0; i < heat_vertices_boundary.size(); i++){
			if (heat_vertices_boundary[i]){
				nb++;
			}
		}

		std::ofstream file(out_path);
		file << nb << std::endl;

		for (int i = 0; i < heat_vertices_boundary.size(); i++){
			if (heat_vertices_boundary[i]){
				file << "1 " << i + 1 << " " << std::endl;;
			}
		}

		file.clear();
		file.close();

		Vector3d1().swap(heat_surface_vertices);
		std::vector<int>().swap(heat_face_id_0);
		std::vector<int>().swap(heat_face_id_1);
		std::vector<int>().swap(heat_face_id_2);
		std::vector<bool>().swap(heat_vertices_boundary);
	}



	void CFSCNC::Output_Obj_Cur_Normals(std::string path)
	{
		Vector3d1 surface_vertices;
		std::vector<std::vector<int>> surface_faces;
		Vector3d1 surface_vertices_normal;
		std::vector<double> surface_vectices_min_curs;
		std::vector<double> surface_vectices_max_curs;
		Vector3d1 max_curvature_directions;
		Vector3d1 min_curvature_directions;

		CGAL_3D_Read_Triangle_Mesh(input_path + "\\path\\" + IntString(cfs_index) + "_split.obj", surface_vertices, surface_faces);
		//CGAL_3D_Mesh_Curvature(surface_vertices, surface_faces, surface_vectices_max_curs, surface_vectices_min_curs, max_curvature_directions, min_curvature_directions, surface_vertices_normal);

		CGAL_Curvature_Mesh(input_path + "\\" + IntString(cfs_index) + "_full.off",
			surface_vertices, surface_vectices_max_curs, surface_vectices_min_curs, max_curvature_directions, min_curvature_directions);
		CGAL_Normal_Mesh(input_path + "\\" + IntString(cfs_index) + "_full.off",
			surface_vertices, surface_vertices_normal);

		std::ofstream file(path);
		file.precision(8);

		file << surface_vertices.size() << std::endl;

		file << "normal_x normal_y normal_z max_cur max_cur_x max_cur_y max_cur_z min_cur min_cur_x min_cur_y min_cur_z" << std::endl;;
		for (int i = 0; i < surface_vertices.size(); i++)
		{
			file << surface_vertices_normal[i][0] << " " << surface_vertices_normal[i][1] << " " << surface_vertices_normal[i][2] << " "
				<< surface_vectices_max_curs[i] << " " << max_curvature_directions[i][0] << " " << max_curvature_directions[i][1] << " " << max_curvature_directions[i][2] << " "
				<< surface_vectices_min_curs[i] << " " << min_curvature_directions[i][0] << " " << min_curvature_directions[i][1] << " " << min_curvature_directions[i][2] << std::endl;
		}

		file.clear();
		file.close();

	}

	void CFSCNC::Output_Path_with_Normal(std::string path)
	{
		std::cout << input_path + "\\" + IntString(cfs_index) + ".off" << std::endl;

		Vector3d1 normals;
		CGAL_Normal_Mesh(input_path + "\\" + IntString(cfs_index) + "_full.off", single_final_path, normals);

		int index = -1;
		double max_y = -100000.0;

		for (int i = 0; i < single_final_path.size(); i++)
		{
			if (single_final_path[i][1]>max_y){
				max_y = single_final_path[i][1];
				index = i;
			}
			Vector3d n = single_final_path[i] + SetVectorLength(normals[i],2.0)+Vector3d(0.0,-2.0,0.0);
		}

		index = 0;
		
		//std::ofstream file0("D:\\r.txt");

		std::ofstream file(path);
		file.precision(8);
		file << toolpath_size << std::endl;
		file << single_final_path.size() << std::endl;

		for (int i = index; i < single_final_path.size(); i++)
		{
			file << single_final_path[i][0] << " " << single_final_path[i][1] << " " << single_final_path[i][2] << " " << normals[i][0] << " " << normals[i][1] << " " << normals[i][2] << std::endl;
			//file << single_path[i][0] << " " << single_path[i][1] << " " << single_path[i][2] << std::endl;
			//file0 << single_final_path[i][0] << " " << single_final_path[i][1] << " " << single_final_path[i][2] << std::endl;
		}

		for (int i = 0; i < index; i++)
		{
			file << single_final_path[i][0] << " " << single_final_path[i][1] << " " << single_final_path[i][2] << " " << normals[i][0] << " " << normals[i][1] << " " << normals[i][2] << std::endl;
			//file << single_path[i][0] << " " << single_path[i][1] << " " << single_path[i][2] << std::endl;
		}

		file.clear();
		file.close();

		//file0.clear();
		//file0.close();
	}

	void SmoothTheNormals(Vector3d1 &normal)
	{
		Vector3d1 new_normals;
		new_normals.push_back(normal[0]);

		for (int i = 1; i < normal.size() - 1; i++)
		{
			Vector3d n0 = normal[i - 1];
			Vector3d n1 = normal[i];
			Vector3d n2 = normal[i + 1];

			Vector3d n(0.0, 0.0, 0.0);
			n = n + n0;
			n = n + n1;
			n = n + n2;
			n[0] = n[0] / 3.0;
			n[1] = n[1] / 3.0;
			n[2] = n[2] / 3.0;

			new_normals.push_back(n);
		}

		new_normals.push_back(normal[normal.size() - 1]);

		normal = new_normals;
	}

	void CFSCNC::Output_Path2(std::string path)
	{
		///////////////////////////////////////////////////////////////////////////////
		std::cout << "Path: " << input_path + "0.off" << std::endl;

		Vector3d1().swap(single_final_path_normal);
		Vector3d1().swap(single_final_path_RMDF_normal);

		//single_final_path_normal
		//Vector3d1 normals;
		CGAL_Normal_Mesh(input_path + "0.off", single_final_path, single_final_path_normal);

		//single_final_path_normal = normals;
		for (int i = 0; i < 200; i++)
			SmoothTheNormals(single_final_path_normal);

		//single_final_path_RMDF_normal
		for (int i = 0; i < single_final_path.size() - 1; i++)
		{
			Vector3d n = getCrossproduct(single_final_path_normal[i], single_final_path[i + 1] - single_final_path[i]);
			//	inline Vector3d RotationAxis(Vector3d p, double angle, Vector3d n)
			single_final_path_RMDF_normal.push_back(RotationAxis(single_final_path_normal[i], MM_PI / 10, n));
		}
		single_final_path_RMDF_normal.push_back(single_final_path_normal[single_final_path_normal.size() - 1]);

		for (int i = 0; i < 200; i++)
			SmoothTheNormals(single_final_path_RMDF_normal);

		//cc=>cl
		//single_final_CL_path

		for (int i = 0; i < single_final_path.size(); i++)
		{
			Vector3d surface_normal = single_final_path_normal[i];
			SetVectorLength(surface_normal, 2.0);

			Vector3d drill_orientation = single_final_path_RMDF_normal[i];
			SetVectorLength(drill_orientation, 2.0);

			Vector3d center = single_final_path[i] - surface_normal + drill_orientation;
			single_final_CL_path.push_back(center);
		}

		std::ofstream file(path);
		file.precision(8);
		file << toolpath_size << std::endl;

		file << single_final_path.size() << std::endl;
		for (int i = 0; i < single_final_path.size(); i++)
		{
			file << single_final_path[i][0] << " " << single_final_path[i][1] << " " << single_final_path[i][2] <<
				" " << single_final_path_normal[i][0] << " " << single_final_path_normal[i][1] << " " << single_final_path_normal[i][2] <<
				" " << single_final_path_RMDF_normal[i][0] << " " << single_final_path_RMDF_normal[i][1] << " " << single_final_path_RMDF_normal[i][2] <<
				" " << single_final_CL_path[i][0] << " " << single_final_CL_path[i][1] << " " << single_final_CL_path[i][2] << std::endl;
			//file << single_final_path[i][0] << " " << single_final_path[i][1] << " " << single_final_path[i][2] << std::endl;
		}

		file.clear();
		file.close();

		//getCrossproduct
	}


	




	void CFSCNC::OutputStripNGC(const std::string path, const Vector3d1 &offsets, bool name_b)
	{

		std::ofstream debug_file(path);

		debug_file << "%" << std::endl;
		debug_file << "G20" << std::endl;
		debug_file << "G90 G94 G40 G17" << std::endl;
		debug_file << "G53 G0 Z0" << std::endl;
		debug_file << "G0 A0 B0" << std::endl;
		debug_file << "G0 X0" << std::endl;
		debug_file << "M5" << std::endl;
		debug_file << "M0" << std::endl;
		debug_file << "T4 M6" << std::endl;
		debug_file << "G43" << std::endl;
		debug_file << "S8500 M3" << std::endl;
		debug_file << "F20" << std::endl;


		for (int j = 0; j < offsets.size(); j++)
		{
			auto s = offsets[j];
			debug_file << "G1 X" + std::to_string(s[0]) + " Y" + std::to_string(s[1]) + " Z" + std::to_string(s[2]) << std::endl;
		}

		debug_file << "G49" << std::endl;
		debug_file << "G53 G0 Z0" << std::endl;
		debug_file << "M5" << std::endl;
		debug_file << "M30" << std::endl;
		debug_file << "%" << std::endl;


		debug_file.clear();
		debug_file.close();
	}



	void CFSCNC::Output_Path1(std::string path)
	{
		//output most hightest path
		///////////////////////////////////////////////////////////////////////////////

		int index = -1;
		double max_y = -10000000.0;
		for (int i = 0; i < single_final_path.size(); i++)
		{
			if (single_final_path[i][1]>max_y)
			{
				max_y = single_final_path[i][1];
				index = i;
			}
		}
		Vector3d1 new_points;
		for (int i = index; i < single_final_path.size(); i++)
		{
			new_points.push_back(single_final_path[i]);
		}
		for (int i = 0; i < index; i++)
		{
			new_points.push_back(single_final_path[i]);
		}

		Vector3d1().swap(single_final_path);
		single_final_path = new_points;
		///////////////////////////////////////////////////////////////////////////////

		std::cout << "Path: " << input_path + "0.off" << std::endl;

		Vector3d1 normals;
		CGAL_Normal_Mesh(input_path + "0.off", single_final_path, normals);

		single_final_path_normal = normals;

		for (int i = 0; i < 1; i++)
			SmoothTheNormals(single_final_path_normal);

		std::ofstream file(path);

		file.precision(8);

		file << toolpath_size << std::endl;

		file << single_final_path.size() << std::endl;
		for (int i = 0; i < single_final_path.size(); i++)
		{
			file << single_final_path[i][0] << " " << single_final_path[i][1] << " " << single_final_path[i][2] << " " << single_final_path_normal[i][0] << " " << single_final_path_normal[i][1] << " " << single_final_path_normal[i][2] << std::endl;
			//file << single_final_path[i][0] << " " << single_final_path[i][1] << " " << single_final_path[i][2] << std::endl;
		}

		file.clear();
		file.close();
	}

	void CFSCNC::Load_Zigzag_Path(std::string path)
	{
		//zigzag_final_path
		int nb;
		std::ifstream file(path, std::ios::in);

		file >> nb;
		for (int i = 0; i < nb; i++)
		{
			double x, y, z;
			file >> x >> y >> z;
			zigzag_final_path.push_back(Vector3d(x, y, z));
		}


		file.clear();
		file.close();
	}

	void CFSCNC::Load_RMDF_Path(std::string path)
	{
		std::ifstream file(path, std::ios::in);

		int nb = 0;
		file >> toolpath_size >> nb;

		for (int i = 0; i < nb; i++)
		{
			double x, y, z;
			double n_x, n_y, n_z;
			double rmdf_n_x, rmdf_n_y, rmdf_n_z;
			file >> x >> y >> z;
			file >> n_x >> n_y >> n_z;
			file >> rmdf_n_x >> rmdf_n_y >> rmdf_n_z;

			single_final_path.push_back(Vector3d(x, y, z));
			single_final_path_normal.push_back(Vector3d(n_x, n_y, n_z));
			single_final_path_RMDF_normal.push_back(Vector3d(rmdf_n_x, rmdf_n_y, rmdf_n_z));

		}


		file.clear();
		file.close();

	}

	void CFSCNC::Load_Path(std::string path)
	{
		std::ifstream file(path, std::ios::in);

		int nb = 0;
		file >> max_scallop >> drill_radius;
		file >> toolpath_size >> nb;

		Vector3d1 normals;
		for (int i = 0; i < nb; i++)
		{
			double x, y, z;
			bool b;
			double n_x, n_y, n_z;
			file >> x >> y >> z >>n_x>>n_y>>n_z>> b;
			single_path.push_back(Vector3d(x, y, z));
			normals.push_back(Vector3d(n_x, n_y, n_z));
			single_path_fixed_label.push_back(b);
		}

		file.clear();
		file.close();
	}

	bool CFSCNC::Load_Final_Path(std::string path)
	{
		std::ifstream file(path, std::ios::in);

		if (file)
		{
			int nb = 0;
			file >> max_scallop >> drill_radius;
			file >> nb;
			for (int i = 0; i < nb; i++)
			{
				double x, y, z;
				bool b;
				file >> x >> y >> z>>b;
				single_final_path.push_back(Vector3d(x, y, z));

				single_path_fixed_label.push_back(b);
			}
			return true;
		}
		else
		{
			file.clear();
			file.close();

			std::cout << "Load final path error: ???" << std::endl;

			return false;
		}
	}

	



	void CFSCNC::Extract_ISOContours_From_Heat(std::string obj_path, std::string dist_path, std::vector<double> &dists, Vector3d2 &offsets)
	{
		if (std::ifstream(dist_path, std::ios::in) && std::ifstream(obj_path, std::ios::in))
		{
			Vector3d1 heat_surface_vertices;
			std::vector<int> heat_face_id_0;
			std::vector<int> heat_face_id_1;
			std::vector<int> heat_face_id_2;

			CGAL_3D_Read_Triangle_Mesh(obj_path, heat_surface_vertices, heat_face_id_0, heat_face_id_1, heat_face_id_2);
			std::vector<double> vertices_d;
			std::vector<int> vertices_index;

			std::ifstream file(dist_path);
			for (int i = 0; i < heat_surface_vertices.size(); i++){
				double d;
				int index;
				file >> d>>index;
				vertices_d.push_back(d);
				vertices_index.push_back(index);
			}
			file.clear();
			file.close();

			toolpath_size = 2.0*sqrt(2.0*max_scallop*drill_radius);
			toolpath_size = ComputeGapFromScallop(0.0, 2.0, max_scallop);
			//toolpath_size = 0.537401/2.0;
			
			double distance = toolpath_size / 2.0;
			int index = 0;
			while (CGAL_3D_Mesh_Extract_Isoline(heat_surface_vertices, heat_face_id_0, heat_face_id_1, heat_face_id_2, vertices_d, distance, offsets))
			{
				for (int i = 0; i < offsets.size() - index; i++) dists.push_back(distance);
				index = offsets.size();
				distance = distance + toolpath_size;
			}
		}
	}



	void CFSCNC::LoadOffsetsFiles(std::vector<string> &offsets_files, std::string path)
	{
		std::ifstream file(path, std::ios::in);
		int nb;
		file >> nb;

		for (int i = 0; i < nb; i++)
		{
			std::string offset_path;
			file >> offset_path;
			offsets_files.push_back(offset_path);
		}
		file.clear();
		file.close();
	}

	void CFSCNC::LoadContour(std::string path, Vector3d2 &offsets, std::vector<double> &dists)
	{
		std::ifstream file(path, std::ios::in);

		if (!file){
			std::cout << "input contour error..." << std::endl;
			return;
		}

		file >> max_scallop;
		file >> drill_radius;
		file >> toolpath_size;
		int contour_number;
		file >> contour_number;

		for (int i = 0; i < contour_number; i++){
			int point_number;
			double distance;
			file >> point_number;
			file >> distance;
			Vector3d1 one_path;
			for (int j = 0; j < point_number; j++){
				double x, y, z;
				file >> x >> y >> z;
				one_path.push_back(Vector3d(x, y, z));
			}

			if (one_path.size() > 3 && Circuit::GetTotalLength(one_path)>toolpath_size*4.0){
				dists.push_back(distance);
				offsets.push_back(one_path);
			}
		}
	}

	void CFSCNC::LoadContours(Vector3d3 &offsetses, Vector3d2 &offsets, bool re_running)
	{
		DWORD start_time = GetTickCount();
		//DWORD end_time = GetTickCount();
		//std::cout << "[TIME] LoadContours: " << (end_time - start_time) / 1000.0 << std::endl;

		std::vector<double> dists;
		if (whether_using_heat)
		{
			int xin = 1;

			if (xin)
			{
				if (!std::ifstream(input_path + "\\path\\heat.offsets", std::ios::in))
				{
					string com_isoline = "";
					com_isoline += "D:\\task2\\CNC\\step_2_1_compute_isolines\\x64\\Release\\isoline.exe";
					com_isoline += " ";
					com_isoline += input_path + "\\path\\" + IntString(cfs_index) + "_heat_sub_split.obj";
					com_isoline += " ";
					//com_isoline += DoubleString(max_scallop*3.5);
					com_isoline += DoubleString(max_scallop);
					com_isoline += " ";
					com_isoline += DoubleString(drill_radius);
					com_isoline += " ";
					com_isoline += input_path + "\\path\\heat.dist";
					com_isoline += " ";
					com_isoline += input_path + "\\path\\heat.offsets";
					system(com_isoline.c_str());

					LoadContour(input_path + "\\path\\heat.offsets", offsets, dists);
				}
				else
				{
					LoadContour(input_path + "\\path\\heat.offsets", offsets, dists);
				}
				
				Vector3d2 offsets0 = offsets;
				Vector3d2().swap(offsets);
				std::vector<double> dists0 = dists;
				std::vector<double>().swap(dists);


				for (int i = 0; i < offsets0.size(); i++)
				{
					Strip::SmoothingLines(offsets0[i], 3);
				}

				for (int i = 0; i < offsets0.size(); i++)
				{
					if (offsets0[i].size()>3 && Circuit::GetTotalLength(offsets0[i]) > toolpath_size*0.5)
					{
						offsets.push_back(offsets0[i]);
						dists.push_back(dists0[i]);
					}
				}
		

				Vector3d2().swap(offsets0);
				std::vector<double>().swap(dists0);
			}
			else
			{
				if (!std::ifstream(input_path + "\\path\\heat.offsets", std::ios::in))
				{
					Extract_ISOContours_From_Heat(input_path + "\\path\\" + IntString(cfs_index) + "_heat_sub_split.obj",
						input_path + "\\path\\heat.dist", dists, offsets);

					Vector3d2 offsets0 = offsets;
					Vector3d2().swap(offsets);
					std::vector<double> dists0 = dists;
					std::vector<double>().swap(dists);

					for (int i = 0; i < offsets0.size(); i++)
					{
						if (offsets0[i].size()>3 && Circuit::GetTotalLength(offsets0[i]) > toolpath_size*4.0)
						{
							offsets.push_back(offsets0[i]);
							dists.push_back(dists0[i]);
						}
					}
					Vector3d2().swap(offsets0);
					std::vector<double>().swap(dists0);

					//output
					/////////////////////////////////////////////////////
					std::string path(input_path + "\\path\\heat.offsets");
					std::ofstream output_file(path);
					output_file << max_scallop << std::endl;
					output_file << drill_radius << std::endl;
					output_file << toolpath_size << std::endl;
					output_file << offsets.size() << std::endl;
					for (int i = 0; i < offsets.size(); i++)
					{
						output_file << offsets[i].size() << " " << dists[i] << std::endl;
						for (int j = 0; j < offsets[i].size(); j++)
							output_file << offsets[i][j][0] << " " << offsets[i][j][1] << " " << offsets[i][j][2] << std::endl;
					}

					output_file.clear();
					output_file.close();
					/////////////////////////////////////////////////////
				}
				else{
					LoadContour(input_path + "\\path\\heat.offsets", offsets, dists);
				}
			}

		
		}
		else{

			std::vector<string> offsets_files;
			LoadOffsetsFiles(offsets_files, input_path + "\\path\\" + IntString(cfs_index) + ".offsets_path");

			for (int iter = 0; iter < offsets_files.size(); iter++){
				string path = offsets_files[iter];

				std::cout << "*****************************************" << std::endl;
				std::cout << "Cout: " << path << std::endl;
				//construct the offsets
				LoadContour(offsets_files[iter], offsets, dists);
			}
		}

		//smoothing
		/****************************************************************************************************/
		for (int i = 0; i < offsets.size(); i++)
		{
			//std::cout << i << " " << offsets.size() << std::endl;
			//Circuit::Smoothing(offsets[i],0.01);
		}
		
		/****************************************************************************************************/

		/****************************************************************************************************/
		double min_d = toolpath_size/2.0;
		
		std::vector<bool> dist_labels(dists.size(),true);

		while (true)
		{
			Vector3d2 pathes;

			for (int i = 0; i < dists.size(); i++){

				if (hpcg::areAlmostEqual(dists[i], min_d))
				{
					pathes.push_back(offsets[i]);

					dist_labels[i] = false;
				}
			}

			if (pathes.size()>0)
				offsetses.push_back(pathes);
			bool goon = false;
			for (int i = 0; i < dists.size(); i++)
			{
				if (dist_labels[i])
				{
					goon = true;
				}
			}
			if (!goon)break;

			min_d += toolpath_size;
		}


		/****************************************************************************************************/
		//uniform direction
		//have not considered there are two circles of the first layer

		//if (corrent_direction) UnifiedContourDirection(toolpath_size, offsetses, offsets);

		Vector3d3().swap(offsetses);

		min_d = toolpath_size / 2.0;

		for (int i = 0; i < dists.size(); i++)
		{
			dist_labels[i] = true;
		}

		while (true)
		{
			Vector3d2 pathes;

			for (int i = 0; i < dists.size(); i++)
			{
				if (hpcg::areAlmostEqual(dists[i], min_d))
				{
					pathes.push_back(offsets[i]);
					dist_labels[i] = false;
				}
			}

			if (pathes.size()>0)
				offsetses.push_back(pathes);
			//else
				//break;

			bool goon = false;
			for (int i = 0; i < dists.size(); i++)
			{
				if (dist_labels[i])
				{
					goon = true;
				}
			}
			if (!goon)break;

			min_d += toolpath_size;
		}

		Vector3d2().swap(offsets);

		for (int i = 0; i < offsetses.size(); i++)
		{
			for (int j = 0; j < offsetses[i].size(); j++)
			{
				offsets.push_back(offsetses[i][j]);
			}
		}

		DWORD end_time = GetTickCount();
		std::cout << "[TIME] LoadContours: " << (end_time - start_time) / 1000.0 << std::endl;

		//construct the offsetses
		/*
		Vector3d2 pathes;
		double distance = d[0];
		for (int i = 0; i < d.size(); i++)
		{
			if (hpcg::areAlmostEqual(d[i], distance))
			{
				pathes.push_back(offsets[i]);
			}
			else
			{
				offsetses.push_back(pathes);
				Vector3d2().swap(pathes);
				pathes.push_back(offsets[i]);
				distance = d[i];
			}
		}
		offsetses.push_back(pathes);
		Vector3d2().swap(pathes);
		*/





		
		/*
		Vector3d3().swap(offsetses);
		distance = d[0];
		for (int i = 0; i < d.size(); i++)
		{
			if (hpcg::areAlmostEqual(d[i], distance))
			{
				pathes.push_back(offsets[i]);
			}
			else
			{
				offsetses.push_back(pathes);
				Vector3d2().swap(pathes);
				pathes.push_back(offsets[i]);
				distance = d[i];
			}
		}
		offsetses.push_back(pathes);
		*/
	}



	//generate data
	Vector2d1 GenerateConcavityData(int nb, int period, double max_r, double min_r)
	{
		Vector2d1 cnc_path;

		for (int i = 0; i < nb; i++)
		{
			double angle = (double)i / (double)nb*MM_PI*2.0;

			double r = (cos(angle*period) + 1.0) / 2.0*(max_r - min_r) + min_r;

			double x, y;
			x = r*cos(angle);
			y = r*sin(angle);
			cnc_path.push_back(Vector2d(x, y));
		}
		return cnc_path;
	}

	void CFSCNC::GenerateAdaptiveCurvatureMesh()
	{
		std::string path = "D:\\CNCProduction\\new_setup\\Concavity\\shape\\14.obj";

		//get boundary
		////////////////////////////////////////////////////////
		Vector3d1 surface_vertices;
		std::vector<int> face_id_0;
		std::vector<int> face_id_1;
		std::vector<int> face_id_2;
		std::vector<int> vertices_index;

		CGAL_3D_Read_Triangle_Mesh(path, surface_vertices, face_id_0, face_id_1, face_id_2);
		std::vector<bool> vertices_boundary;
		Vector3d2 boundarys;
		CGAL_3D_Triangle_Mesh_Boundary(surface_vertices, face_id_0, face_id_1, face_id_2, vertices_boundary);
		CGAL_3D_Triangle_Mesh_Boundary(surface_vertices, face_id_0, face_id_1, face_id_2, boundarys);

		//compute neighboring
		////////////////////////////////////////////////////////
		std::vector<std::vector<int>> neighs;
		CGAL_3D_Triangle_Mesh_Vecs_Neighbors(surface_vertices, face_id_0, face_id_1, face_id_2, neighs);

		//hull outside points
		////////////////////////////////////////////////////////
		Vector2d1 hull_outside;
		for (int i = 0; i < boundarys[0].size(); i++)
			hull_outside.push_back(Vector2d(boundarys[0][i][0], boundarys[0][i][2]));

		//sampling
		////////////////////////////////////////////////////////
		int control_points_nb = 5;
		std::vector<double> hull_heights(hull_outside.size(), 0.0);
		Vector2d1 hull_inside(1,Vector2d(0.0,0.0));
		//CGAL_2D_Polygon_Dart_Sampling(hull_outside, 8.0, hull_inside);

		double max_height = 9.0;
		for (int i = 0; i < hull_inside.size(); i++)
		{
			//hull_heights.push_back(max_height*rand() / double(RAND_MAX));
			hull_heights.push_back(max_height);
		}
		std::vector<double> hull_heights0(hull_outside.size() + hull_inside.size(), 0.0);

		//hull_heights0
		CGAL_2D_Polygon_Triangulation(hull_outside, hull_inside, hull_heights, "off", "Z:\\Documents\\Windows\\SmartSFC\\workspace\\CFS\\haisen.off");

		Vector3d1 inside_vecs;
		for (int i = 0; i < surface_vertices.size(); i++)
		{
			if (!vertices_boundary[i])
				inside_vecs.push_back(surface_vertices[i]);
		}
		Vector3d1 inside_normal_vecs(inside_vecs.size(), Vector3d(0.0, 1.0, 0.0));
		Vector3d1 inters;
		CGAL_3D_Intersection_Ray_Mesh(inside_vecs, inside_normal_vecs, "Z:\\Documents\\Windows\\SmartSFC\\workspace\\CFS\\haisen.off", inters);

		std::vector<double> desired_heights;
		std::vector<double> current_heights;
		int index = 0;
		for (int i = 0; i < surface_vertices.size(); i++)
		{
			if (vertices_boundary[i])
				desired_heights.push_back(0.0);
			else{
				desired_heights.push_back(inters[index][1]);
				index++;
			}
			current_heights.push_back(0.0);
			surface_vertices[i][1] = desired_heights[desired_heights.size() - 1];
		}

		CGAL_Output_Obj("Z:\\Documents\\Windows\\SmartSFC\\workspace\\CFS\\remesh.obj", surface_vertices, face_id_0, face_id_1, face_id_2);
		CGAL_Mesh_Laplace_Smooth("Z:\\Documents\\Windows\\SmartSFC\\workspace\\CFS\\remesh.obj", "Z:\\Documents\\Windows\\SmartSFC\\workspace\\CFS\\haisen.obj", 20);
		return;

	}



	void Laplace_Mesh_Vertex_Values(Vector3d1 &vecs, std::vector<int> &face_id_0, std::vector<int> &face_id_1, std::vector<int> &face_id_2, std::vector<double> &values, int iterations)
	{
		std::vector<std::vector<int>> neighs;
		CGAL_3D_Triangle_Mesh_Vecs_Neighbors(vecs, face_id_0, face_id_1, face_id_2, neighs);

		std::vector<double> iter_values(vecs.size(), 0.0);

		for (int i = 0; i < iterations; i++)
		{
			for (int j = 0; j < vecs.size(); j++)
			{
				for (int k = 0; k < neighs[j].size(); k++)
				{
					iter_values[j] += values[neighs[j][k]];
				}
				iter_values[j] = iter_values[j] / (double)neighs[j].size();

				//iter_values[j] = 0.5*iter_values[j] + 0.5*values[i];
			}

			std::vector<double>().swap(values);
			values = iter_values;
			for (int j = 0; j < vecs.size(); j++)
				iter_values[j] = 0.0;
		}
	}


	std::vector<double> CFSCNC::Heat_to_Gradient_0(std::string path, std::string obj_path, std::string full_obj_path, std::string heat_path, std::string scale_path)
	{
		/*******************************************/
		Vector3d1 vecs;
		std::vector<bool> vecs_boundary;
		std::vector<int> face_id_0;
		std::vector<int> face_id_1;
		std::vector<int> face_id_2;
		/*******************************************/
		CGAL_3D_Read_Triangle_Mesh(obj_path, vecs, face_id_0, face_id_1, face_id_2);

		CGAL_3D_Triangle_Mesh_Boundary(vecs, face_id_0, face_id_1, face_id_2, vecs_boundary);

		/*******************************************/
		std::vector<double> dists;
		/*******************************************/
		std::ifstream file(heat_path);
		for (int i = 0; i < vecs.size(); i++)
		{
			double d;
			int index;
			file >> d;
			file >> index;
			dists.push_back(d);
		}
		file.clear();
		file.close();

		/*******************************************/
		Vector3d1 vecs_gradients(vecs.size(), Vector3d(0.0, 0.0, 0.0));
		std::vector<double> vecs_gradients_value(vecs.size(), 0.0);

		std::vector<double> areas(vecs.size(), 0.0);
		/*******************************************/

		/*******************************************/
		Vector3d1 face_gradients;
		/*******************************************/
		for (int i = 0; i < face_id_0.size(); i++)
		{
			int index_0 = face_id_0[i];
			int index_1 = face_id_1[i];
			int index_2 = face_id_2[i];

			Vector3d v_0 = vecs[index_0];
			Vector3d v_1 = vecs[index_1];
			Vector3d v_2 = vecs[index_2];

			Vector3d   n = Math::GetCrossproduct(v_0 - v_2, v_1 - v_2);
			Math::SetVectorLength(n, 1.0);

			double area = CGAL_3D_One_Triangle_Area(v_0, v_1, v_2);

			Vector3d gradient(0.0, 0.0, 0.0);

			double d_0 = (float)dists[index_0];
			double d_1 = (float)dists[index_1];
			double d_2 = (float)dists[index_2];

			gradient += (float)dists[index_0] * Math::GetCrossproduct(n, v_2 - v_1);
			gradient += (float)dists[index_1] * Math::GetCrossproduct(n, v_0 - v_2);
			gradient += (float)dists[index_2] * Math::GetCrossproduct(n, v_1 - v_0);

			gradient = gradient / (float)(2.0 * area);


			face_gradients.push_back(gradient);

			vecs_gradients[index_0] += (float)area*gradient;
			areas[index_0] += area;

			vecs_gradients[index_1] += (float)area*gradient;
			areas[index_1] += area;

			vecs_gradients[index_2] += (float)area*gradient;
			areas[index_2] += area;

			//Vector3d middle = (v_0 + v_1 + v_2) / (float)3.0;
			//Math::SetVectorLength(gradient,0.3);
			//CGAL_Export_Path_Segment(export_fie,export_int,"name_"+Math::IntString(i),1.0,0.0,0.0,middle,middle+gradient,0.02);
		}

		std::vector<double> base_vertices_max_cur;
		std::vector<double> base_vertices_min_cur;
		Vector3d1 base_vertices_max_cur_direction;
		Vector3d1 base_vertices_min_cur_direction;

		CGAL_Curvature_Mesh(full_obj_path,
			vecs, base_vertices_max_cur, base_vertices_min_cur, base_vertices_max_cur_direction, base_vertices_min_cur_direction);

		/*********************/
		/*
		std::ifstream input_file("Z:\\Documents\\Windows\\SmartSFC\\workspace\\CFS\\0_heat_factor.source");

		int nb = 0;
		input_file >> nb;

		for (int i = 0; i < nb; i++)
		{
		double max_cur;
		double min_cur;
		double max_x;
		double max_y;
		double max_z;
		double min_x;
		double min_y;
		double min_z;
		input_file >> max_cur >> min_cur >> max_x >> max_y >> max_z >> min_x >> min_y >> min_z;

		base_vertices_max_cur.push_back(max_cur);
		base_vertices_min_cur.push_back(min_cur);
		base_vertices_max_cur_direction.push_back(Vector3d(max_x, max_y, max_z));
		base_vertices_min_cur_direction.push_back(Vector3d(min_x, min_y, min_z));
		}

		input_file.clear();
		input_file.close();*/
		/*********************/

		int convex_nb = 0;
		int concave_nb = 0;

		for (int i = 0; i < vecs.size(); i++)
		{
			vecs_gradients[i] = vecs_gradients[i] / (float)areas[i];

			//Math::SetVectorLength(vecs_gradients[i], 0.3);
			//CGAL_Export_Path_Segment(export_fie, export_int, "name_" + Math::IntString(i), 1.0, 0.0, 0.0, vecs[i], vecs[i] + vecs_gradients[i], 0.02);

			double min_cur = base_vertices_min_cur[i];
			double max_cur = base_vertices_max_cur[i];
			Vector3d min_cur_dir = base_vertices_min_cur_direction[i];
			Vector3d max_cur_dir = base_vertices_max_cur_direction[i];

			double max_angle = Math::GetAngleBetween(max_cur_dir, vecs_gradients[i]);
			double min_angle = Math::GetAngleBetween(min_cur_dir, vecs_gradients[i]);

			Math::SetVectorLength(min_cur_dir, 0.3);
			Math::SetVectorLength(max_cur_dir, 0.3);

			if (max_angle > MM_PI / 2.0) max_angle = MM_PI - max_angle;
			if (min_angle >MM_PI / 2.0) min_angle = MM_PI - min_angle;

			vecs_gradients_value[i] = (1.0 - max_angle / (MM_PI / 2.0))*max_cur + (1.0 - min_angle / (MM_PI / 2.0))*min_cur;

			if (vecs_gradients_value[i] < 0)
			{
				concave_nb++;
			}
			else
			{
				convex_nb++;
			}
		}

		std::cout << "Convex Nb: " << convex_nb << std::endl;
		std::cout << "Concave Nb: " << concave_nb << std::endl;

		std::ofstream ofile(scale_path);
		ofile << Math::IntString(vecs.size()) << std::endl;

		double flat_width = ComputeGapFromScallop(0.0, 2.0, max_scallop);
		Vector3d1 colors;
		std::vector<double> scales;
		for (int i = 0; i < vecs.size(); i++)
		{
			double w = ComputeGapFromScallop(vecs_gradients_value[i], 2.0, max_scallop);

			if (vecs_boundary[i])
				w = flat_width;

			if (w > flat_width*3.0) w = flat_width*2.0;
			if (w < 0.25) w = 0.25;

			w = w*0.95;

			double factor = w * w / flat_width / flat_width;

			//factor = std::pow(factor,2.2);

			//ofile << factor << std::endl;
			scales.push_back(factor);
		}

		Laplace_Mesh_Vertex_Values(vecs, face_id_0, face_id_1, face_id_2, scales, 5);

		for (int i = 0; i < vecs.size(); i++)
		{
			ofile << scales[i] << std::endl;

			double r, g, b;
			double d = (scales[i] - 0.28) / (1.0 - 0.28);
			Math::ColorMapping(d, r, g, b);
			//colors.push_back(Vector3d(r, g, b));

			if (scales[i]>1.0)
				colors.push_back(Vector3d(1.0, 0.0, 0.0));
			else
				colors.push_back(Vector3d(0.0, 1.0, 0.0));

		}

		ofile.clear();
		ofile.close();

		//smooth
		//handle the boundary 

		file.clear();
		file.close();

		CGAL_Output_Obj(path + "\\path\\density_map.obj", vecs, colors, face_id_0, face_id_1, face_id_2);

		return scales;
	}

}
