// SmartSFC.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "iostream"
#include "CFSCNC.h"
#include "RoughMachining.h"
#include "AM.h"

using namespace cnc;

int main(int argc, char* argv[])
{

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

