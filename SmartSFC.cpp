// SmartSFC.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "iostream"
#include "CFSCNC.h"
#include "RoughMachining.h"
#include "AM.h"
#include "MshLoader.h"

using namespace PyMesh;
using namespace cnc;


void mold()
{
	std::string input_tet_file("E:\\Dropbox\\Mold\\TetWild\\fig19\\input_.msh");
	std::string input_surf_file("E:\\Dropbox\\Mold\\TetWild\\fig19\\input__sf.obj");

	Vector3d1 surf_vecs;
	Vector1i2 surf_faces(3, Vector1i1());

	CGAL_3D_Read_Triangle_Mesh(input_surf_file, surf_vecs, surf_faces[0], surf_faces[1], surf_faces[2]);
	MshLoader ml(input_tet_file);

	std::cerr << ml.m_nodes.cols() << std::endl;
	std::cerr << ml.m_nodes.rows() << std::endl;
	std::cerr << ml.m_elements.cols() << std::endl;
	std::cerr << ml.m_elements.rows() << std::endl;

	int tet_number = ml.m_elements.rows() / 4;
	int vec_number = ml.m_nodes.rows() / 3;
	for (int i = 0; i < tet_number; i++)
	{
		auto index_0 = ml.m_elements[i * 4];
		auto index_1 = ml.m_elements[i * 4 + 1];
		auto index_2 = ml.m_elements[i * 4 + 2];
		auto index_3 = ml.m_elements[i * 4 + 3];
		for (int j = 0; j <= 3; j++)
		{
			auto index = ml.m_elements[i * 4 + j];
			if (!(index >= 0 && index < vec_number))
			{
				//Functs::MAssert("if (!(index >= 0 && index < vec_number))");
			}
		}
	}

	Vector3d1 obj_vectors(vec_number, Vector3d());
	for (int i = 0; i < vec_number; i++)
	{
		obj_vectors[i][0] = ml.m_nodes[i * 3];
		obj_vectors[i][1] = ml.m_nodes[i * 3 + 1];
		obj_vectors[i][2] = ml.m_nodes[i * 3 + 2];
	}

	//Functs::OutputObj3d("E:\\Dropbox\\Task1\\postdoc\\IST\\Desktop\\temp.obj", obj_vectors);


	//std::ofstream export_file("E:\\temp.obj");
	//int export_index = 1;
	//output_point(export_file, export_index,"",0.5,0.5,0.5, Vector3d(0,0,0),0.5);
	//export_file.close();
	system("pause");
};


int main(int argc, char* argv[])
{
	mold();
	return 0;

	//stock_direction / stock_3.stl(25.413, 46.112, 41.005) = > (20.403, 34.292, 21.889)
	//	stock_direction / stock_4.stl(32.281, 45.978, 42.832) = > (13.438, 34.914, 35.563)
	//	stock_direction / stock_5.stl(7.633, 46.132, 28.473) = > (27.152, 34.479, 32.147)
	//	stock_direction / stock_6.stl(10.264, 45.727, 18.948) = > (16.056, 34.964, 38.464)


	//Vector3d1 vec;
	//vec.emplace_back(25.413, 46.112, 41.005);
	//vec.emplace_back(20.403, 34.292, 21.889);
	//vec.emplace_back(32.281, 45.978, 42.832);
	//vec.emplace_back(13.438, 34.914, 35.563);
	//vec.emplace_back(7.633, 46.132, 28.473);
	//vec.emplace_back(27.152, 34.479, 32.147);
	//vec.emplace_back(10.264, 45.727, 18.948);
	//vec.emplace_back(16.056, 34.964, 38.464);

	//for (int i = 0; i < vec.size(); i = i + 2)
	//{
	//	auto a = vec[i + 1] - vec[i];
	//	std::cerr << a[0] << " " << a[1] << " " << a[2] << std::endl;
	//}

	std::cout << "try project " << std::endl;

	if (argc == 2)
	{
		std::string path = argv[1];

		std::cerr << path << std::endl;

		DWORD start_time = GetTickCount();
		RoughMachining rm(path, true, true);
		DWORD end_time = GetTickCount();
		std::cout << "[TIME] tool path generation times " << (end_time - start_time) / 1000.0 << std::endl;
	}

	if (argc == 4)
	{
		std::string path = argv[1];
		double toolpath_size = stod(std::string(argv[2]));
		double layer_thickness = stod(std::string(argv[3]));

		DWORD start_time = GetTickCount();
		RoughMachining rm(path, true, true, toolpath_size, layer_thickness);
		DWORD end_time = GetTickCount();
		std::cout << "[TIME] tool path generation times " << (end_time - start_time) / 1000.0 << std::endl;
	}


	if (argc == 6)
	{
		std::string path = argv[1];
		double toolpath_size = stod(std::string(argv[2]));
		double layer_thickness = stod(std::string(argv[3]));
		bool re_running = stod(std::string(argv[4]));
		bool output_debug = stod(std::string(argv[5]));

		DWORD start_time = GetTickCount();
		RoughMachining rm(path, re_running, output_debug, toolpath_size, layer_thickness);
		DWORD end_time = GetTickCount();
		std::cout << "[TIME] tool path generation times " << (end_time - start_time) / 1000.0 << std::endl;
	}


	if (argc == 7)
	{
		std::string path = argv[1];
		double toolpath_size = stod(std::string(argv[2]));
		double layer_thickness = stod(std::string(argv[3]));
		bool re_running = stod(std::string(argv[4]));
		bool output_debug = stod(std::string(argv[5]));
		bool segment_or_tool_path = stod(std::string(argv[6]));

		DWORD start_time = GetTickCount();
		//RoughMachining rm(path, re_running, output_debug, toolpath_size, layer_thickness, segment_or_tool_path);
		AMVolume am(path, re_running, output_debug, toolpath_size, layer_thickness, segment_or_tool_path);
		DWORD end_time = GetTickCount();
		std::cout << "[TIME] tool path generation times " << (end_time - start_time) / 1000.0 << std::endl;

		//if (segment_or_tool_path) system("pause");
	}

	if (argc == 8)
	{
		std::string path = argv[1];
		double toolpath_size = stod(std::string(argv[2]));
		double layer_thickness = stod(std::string(argv[3]));
		bool re_running = stod(std::string(argv[4]));
		bool output_debug = stod(std::string(argv[5]));
		bool segment_or_tool_path = stod(std::string(argv[6]));
		int index = stod(std::string(argv[7]));

		DWORD start_time = GetTickCount();
		RoughMachining rm(path, re_running, output_debug, toolpath_size, layer_thickness, segment_or_tool_path, index);
		DWORD end_time = GetTickCount();
		std::cout << "[TIME] tool path generation times " << (end_time - start_time) / 1000.0 << std::endl;

		//if (segment_or_tool_path) system("pause");
	}

	//cfscnc.OffsetsBasedCFSLinking(path, true);
	//system("pause");
	return 0;
}

