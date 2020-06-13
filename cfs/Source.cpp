#include "stdafx.h"
#include "CFSCNC.h"
#include <Circuit.h>

namespace cnc {


	//void CFSCNC::CGAL_Cut_Surface()
	//{
	//	std::string full_path = "D:\\CNCProduction\\new_setup\\bunny\\bunny\\0_404_full.obj";
	//	std::string path = "D:\\CNCProduction\\new_setup\\bunny\\bunny\\0_404.obj";

	//	Vector3d1 vecs;
	//	std::vector<int> face_id_0;
	//	std::vector<int> face_id_1;
	//	std::vector<int> face_id_2;
	//	CGAL_3D_Read_Triangle_Mesh(path, vecs, face_id_0, face_id_1, face_id_2);
	//	Vector3d2 boundarys;
	//	CGAL_3D_Triangle_Mesh_Boundary(vecs, face_id_0, face_id_1, face_id_2, boundarys);
	//	Vector3d1 boundary = boundarys[0];

	//	Vector3d1 full_vecs;
	//	std::vector<int> full_face_id_0;
	//	std::vector<int> full_face_id_1;
	//	std::vector<int> full_face_id_2;
	//	CGAL_3D_Read_Triangle_Mesh(full_path, full_vecs, full_face_id_0, full_face_id_1, full_face_id_2);

	//	Vector3d1 projects = CGAL_Project_Points_Onto_Surface(full_vecs, full_face_id_0, full_face_id_1, full_face_id_2, boundary);

	//}


}

