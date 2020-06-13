// SmartSFC.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "iostream"
#include "CFSCNC.h"
#include "RoughMachining.h"

using namespace cnc;

int main(int argc, char* argv[])
{
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
		RoughMachining rm(path, re_running, output_debug, toolpath_size, layer_thickness, segment_or_tool_path);
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
	
	return 0;
}

