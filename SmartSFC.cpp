// SmartSFC.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "iostream"
#include "CFSCNC.h"
#include "RoughMachining.h"
#include "AM.h"
#include "MshLoader.h"
#include "pgl_functs.hpp"


using namespace PyMesh;
using namespace cnc;
//using namespace PGL;


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
					PGL::Functs::MAssert("if (!(index >= 0 && index < vec_number))");
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

	//ray intersection
	auto Ray_Intersection = [](const Tree&tree, const Vector3d& origin, const Vector3d& target)
	{
		Ray_3 ray(Point_3(origin[0], origin[1], origin[2]), Vector_3(target[0], target[1], target[1]));
		return  tree.do_intersect(ray);
	};

	//uniform ray directions
	//control distance between points???
	//control percentage ????
	auto Random_Directions = [](const int& dns)
	{
		//if (dns <= 0) MAssert("if (dns <= 0)");//check validation of dns

		double gaussion_sphere_radius = 1.0;
		double idea_distance = 4 * gaussion_sphere_radius / sqrt(dns);

		int dis_iters = 100;
		Vector3d1 directions;
		for (int i = 0; i < dns; i++)
		{
			for (int j = 0; j < dis_iters; j++)
			{
				double alpha_angle = rand() / double(RAND_MAX) *2.0* PGL::Math_PI;
				double alpha_beta = rand() / double(RAND_MAX) *2.0* PGL::Math_PI;
				auto direction_0 = PGL::Functs::RotationAxis(Vector3d(gaussion_sphere_radius, 0.0, 0.0), alpha_angle, Vector3d(0.0, 1.0, 0.0));
				auto direction_axis = PGL::Functs::GetCrossproduct(direction_0, Vector3d(0.0, 1.0, 0.0));
				auto direction_1 = PGL::Functs::RotationAxis(direction_0, alpha_beta, direction_axis);



				directions.push_back(direction_1);
			}

		}
		return directions;
		//Vector3d RotationAxis(Vector3d p, double angle, Vector3d n)
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
	
	//accessibility analysis for each tet element
	for (int i = 0; i < tet_number; i++)
	{
		Vector3d center = Get_Tet_Center(tet, i);
	}

	//test direction generations
	auto directions = Random_Directions(100);

	std::ofstream file("E:\\Dropbox\\Task1\\postdoc\\IST\\Desktop\\temp.obj");
	int export_index = 1;
	for (int i = 0; i < directions.size(); i++)
	{
		CGAL_Export_Path_Point(file, export_index, "", 0.5, 0.5, 0.5, directions[i],0.05);
	}

	file.clear();
	file.close();

	//void CGAL_Export_Path_Point(std::ofstream &export_file_output, int &export_index,
	//	std::string s_name, double r, double g, double b, Vector3d point, double radius);

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

