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

void MoldAccessibilityAnalysis(string input_tet_file, string input_surf_file)
{
	auto Load_ML = [](string input_tet_file)
	{
		MshLoader ml(input_tet_file);
		int tet_number = ml.m_elements.rows() / 4;
		int vec_number = ml.m_nodes.rows() / 3;
		for (int i = 0; i < tet_number; i++)
		{
			for (int j = 0; j <= 3; j++)
			{
				auto index = ml.m_elements[i * 4 + j];
				if (!(index >= 0 && index < vec_number))
				{
					MAssert("if (!(index >= 0 && index < vec_number))");
				}
			}
		}
		return ml;
	};

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
		for (int j = 0; j <= 3; j++)
		{
			tet.m_nodes[tet.m_elements[i * 4 + j]];
			tet.m_nodes[tet.m_elements[i * 4 + j] + 1];
			tet.m_nodes[tet.m_elements[i * 4 + j] + 2];
		}
	}
};

void mold()
{
	//ray intersection
	std::string input_tet_file("E:\\Dropbox\\Mold\\TetWild\\fig19\\input_.msh");
	std::string input_surf_file("E:\\Dropbox\\Mold\\TetWild\\fig19\\input__sf.obj");

	MoldAccessibilityAnalysis(input_tet_file, input_surf_file);

	//ray mesh intersection
	Polyhedron_3 polyhedron;
	Construct_Polyhedron(polyhedron, input_surf_file);
	Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);
	tree.accelerate_distance_queries();
	Ray_3 ray(Point_3(0.0,0.0,0.0), Vector_3(1.0,0.0,0.0));
	auto inter = tree.do_intersect(ray);
	std::cerr << "Inter: " << inter << std::endl;


	//Vector3d1 surf_vecs;
	//Vector1i2 surf_faces(3, Vector1i1());
	//CGAL_3D_Read_Triangle_Mesh(input_surf_file, surf_vecs, surf_faces[0], surf_faces[1], surf_faces[2]);

	//read tet file
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
}

