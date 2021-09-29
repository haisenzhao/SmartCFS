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
using namespace Math;

//microstructure accessibility analysis
void MoldAccessibilityAnalysis(const string& input_tet_file, const string& input_surf_file)
{
#pragma region Lambda-Functions
	//load tet file
	auto Load_ML = [](const string& input_tet_file)
	{
		//initialization
		MshLoader ml(input_tet_file);
		int tet_number = ml.m_elements.rows() / 4;
		int vec_number = ml.m_nodes.rows() / 3;
		//output some information
		std::cerr << "Tet file path:  " << input_tet_file << std::endl;
		std::cerr << "Tet elements: " << ml.m_elements.rows() << std::endl;
		std::cerr << "Tet nodes:      " << ml.m_nodes.rows() << std::endl;
		//check validation of element index
		std::cerr << "Checking the element index..." << std::endl;
		for (int i = 0; i < tet_number; i++)
		{
			if (i % (tet_number / 100) == 0)std::cerr << (double)i / (double)tet_number << "% ";
			for (int j = 0; j <= 3; j++)
			{
				auto index = ml.m_elements[i * 4 + j];
				if (!(index >= 0 && index < vec_number))
					MAssert("if (!(index >= 0 && index < vec_number))");
			}
		}
		std::cerr << std::endl;
		//return ml
		return ml;
	};

	//get tet center
	auto Get_Tet_Center = [](const MshLoader& tet, const int& tet_index)
	{
		//initialization
		Vector3d center(0.0, 0.0, 0.0);
		//read the coordinate of each vectice
		for (int j = 0; j < 4; j++)
		{
			center[0] += tet.m_nodes[tet.m_elements[tet_index * 4 + j]];
			center[1] += tet.m_nodes[tet.m_elements[tet_index * 4 + j] + 1];
			center[2] += tet.m_nodes[tet.m_elements[tet_index * 4 + j] + 2];
		}
		center[0] = center[0] / 4.0;
		center[1] = center[1] / 4.0;
		center[2] = center[2] / 4.0;
		return center;
	};

	auto Ray_Intersection = []()
	{

	};

#pragma endregion

	//read tet 
	MshLoader tet = Load_ML(input_tet_file);
	int tet_number = tet.m_elements.rows() / 4;
	int vec_number = tet.m_nodes.rows() / 3;

	//build polyhedron_tree
	Polyhedron_3 polyhedron;
	Construct_Polyhedron(polyhedron, input_surf_file);
	Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);
	tree.accelerate_distance_queries();
	
	//Ray_3 ray(Point_3(0.0, 0.0, 0.0), Vector_3(1.0, 0.0, 0.0));
	//auto inter = tree.do_intersect(ray);
	//std::cerr << "Inter: " << inter << std::endl;

	//accessibility analysis for each tet element
	for (int i = 0; i < tet_number; i++)
	{
		Vector3d center = Get_Tet_Center(tet, i);

	}
};

void mold()
{
	//ray intersection
	std::string input_tet_file("E:\\Dropbox\\Mold\\TetWild\\fig19\\input_.msh");
	std::string input_surf_file("E:\\Dropbox\\Mold\\TetWild\\fig19\\input__sf.obj");

	MoldAccessibilityAnalysis(input_tet_file, input_surf_file);
};


int main(int argc, char* argv[])
{
	mold();
	system("pause");
	return 0;
}

