#include "stdafx.h"
#include "cgalpackage.h"


//IO mesh
/***************************************************************************************************/
int get_first_integer(const char *v){
	int ival;
	std::string s(v);
	std::replace(s.begin(), s.end(), '/', ' ');
	sscanf(s.c_str(), "%d", &ival);
	return ival;
}

void CGAL_Load_Obj(std::string path, std::vector<double> &coords, std::vector<int> &tris)
{
	double x, y, z;
	char line[1024], v0[1024], v1[1024], v2[1024];

	// open the file, return if open fails
	FILE *fp = fopen(path.c_str(), "r");
	if (!fp) return;

	while (fgets(line, 1024, fp)){
		if (line[0] == 'v'){
			sscanf(line, "%*s%lf%lf%lf", &x, &y, &z);
			coords.push_back(x);
			coords.push_back(y);
			coords.push_back(z);
		}
		else if (line[0] == 'f'){
			sscanf(line, "%*s%s%s%s", v0, v1, v2);
			tris.push_back(get_first_integer(v0) - 1);
			tris.push_back(get_first_integer(v1) - 1);
			tris.push_back(get_first_integer(v2) - 1);
		}
	}
	fclose(fp);
}

void CGAL_Output_Obj(std::string path, Vector3d1 &vecs)
{
	if (vecs.size() < 3)
	{
		std::cout << "CGAL_Output_Obj error: vecs.size() < 3 " << std::endl;
		return;
	}

	std::ofstream file(path);
	for (int i = 0; i<vecs.size(); i++)
	{
		file << "v " << vecs[i][0] << " " << vecs[i][1] << " " << vecs[i][2]  << std::endl;
	}

	file.close();
}

void CGAL_Output_Obj(std::string path, Vector3d1 &vecs, Vector3d1 &colors, std::vector<std::vector<int>> &face_ids)
{
	std::vector<int> face_id_0;
	std::vector<int> face_id_1;
	std::vector<int> face_id_2;

	for (int i = 0; i < face_ids.size(); i++)
	{
		face_id_0.push_back(face_ids[i][0]);
		face_id_1.push_back(face_ids[i][1]);
		face_id_2.push_back(face_ids[i][2]);
	}
	CGAL_Output_Obj(path, vecs, colors, face_id_0, face_id_1, face_id_2);
}

void CGAL_Output_Obj(std::string path, Vector3d1 &vecs, Vector3d1 &colors, std::vector<int> &face_id_0, std::vector<int> &face_id_1, std::vector<int> &face_id_2)
{

	if (vecs.size() < 3 || colors.size()<3 || face_id_0.size() < 1 || face_id_1.size() < 1 || face_id_2.size() < 1)
	{
		std::cout << "CGAL_Output_Obj error: vecs.size() < 3 || face_id_0.size() < 1 || face_id_1.size() < 1 || face_id_2.size() < 1" << std::endl;
		return;
	}

	for (int i = 0; i < face_id_0.size(); i++)
	{
		int index_0 = face_id_0[i];
		int index_1 = face_id_1[i];
		int index_2 = face_id_2[i];

		if (index_0 < 0 || index_0 >= vecs.size() || index_1 < 0 || index_1 >= vecs.size() || index_2 < 0 || index_2 >= vecs.size())
		{
			std::cout << "CGAL_Output_Obj error: index_0 < 0 || index_0 >= vecs.size() || index_1 < 0 || index_1 >= vecs.size() || index_2 < 0 || index_2 >= vecs.size()" << std::endl;
			return;
		}
	}

	std::ofstream file(path);
	for (int i = 0; i<vecs.size(); i++)
	{
		file << "v " << vecs[i][0] << " " << vecs[i][1] << " " << vecs[i][2] << " " << colors[i][0] << " " << colors[i][1] << " " << colors[i][2] << std::endl;
	}

	for (int i = 0; i < face_id_0.size(); i++)
	{
		int index_0 = face_id_0[i];
		int index_1 = face_id_1[i];
		int index_2 = face_id_2[i];

		if (index_0 != index_1&&index_0 != index_2&&index_1 != index_2)
			file << "f " << index_0 + 1 << " " << index_1 + 1 << " " << index_2 + 1 << std::endl;
	}
	file.close();
}

void CGAL_Output_Obj(std::string path, Vector3d1 &vecs, std::vector<int> &face_id_0, std::vector<int> &face_id_1, std::vector<int> &face_id_2)
{

	if (vecs.size() < 3 || face_id_0.size() < 1 || face_id_1.size() < 1 || face_id_2.size() < 1)
	{
		std::cout << "CGAL_Output_Obj error: vecs.size() < 3 || face_id_0.size() < 1 || face_id_1.size() < 1 || face_id_2.size() < 1" << std::endl;
		return;
	}

	for (int i = 0; i < face_id_0.size(); i++)
	{
		int index_0 = face_id_0[i];
		int index_1 = face_id_1[i];
		int index_2 = face_id_2[i];

		if (index_0 < 0 || index_0 >= vecs.size() || index_1 < 0 || index_1 >= vecs.size() || index_2 < 0 || index_2 >= vecs.size())
		{
			std::cout << "CGAL_Output_Obj error: index_0 < 0 || index_0 >= vecs.size() || index_1 < 0 || index_1 >= vecs.size() || index_2 < 0 || index_2 >= vecs.size()" << std::endl;
			return;
		}
	}

	std::ofstream file(path);
	for (int i = 0; i<vecs.size(); i++)
	{
		Vector3d v = vecs[i];
		//file << "v " << vecs[i][0] << " " << vecs[i][1] << " " << vecs[i][2] << std::endl;
		file << "v " << v[0] << " " << v[1] << " " << v[2] << std::endl;
	}

	for (int i = 0; i < face_id_0.size(); i++)
	{
		int index_0 = face_id_0[i];
		int index_1 = face_id_1[i];
		int index_2 = face_id_2[i];

		if (index_0 != index_1&&index_0 != index_2&&index_1 != index_2)
			file << "f " << index_0 + 1 << " " << index_1 + 1 << " " << index_2 + 1 << std::endl;
	}
	file.close();
}

void CGAL_Output_Obj(std::string path, Vector3d1 &vecs, std::vector<std::vector<int>> &face_ids)
{
	std::vector<int> face_id_0;
	std::vector<int> face_id_1;
	std::vector<int> face_id_2;

	for (int i = 0; i < face_ids.size(); i++)
	{
		face_id_0.push_back(face_ids[i][0]);
		face_id_1.push_back(face_ids[i][1]);
		face_id_2.push_back(face_ids[i][2]);
	}

	CGAL_Output_Obj(path, vecs, face_id_0, face_id_1, face_id_2);
}
void CGAL_Output_Obj(std::string path, Vector3d1 &vecs, std::vector<std::vector<int>> &face_ids, std::vector<int> &triangles_lables, int index)
{
	std::vector<int> face_id_0;
	std::vector<int> face_id_1;
	std::vector<int> face_id_2;

	std::vector<int> lables(vecs.size(), -1);
	for (int i = 0; i < face_ids.size(); i++)
	{
		face_id_0.push_back(face_ids[i][0]);
		face_id_1.push_back(face_ids[i][1]);
		face_id_2.push_back(face_ids[i][2]);
		if (triangles_lables[i] == index)
		{
			lables[face_ids[i][0]] = 0;
			lables[face_ids[i][1]] = 0;
			lables[face_ids[i][2]] = 0;
		}
	}
	Vector3d1 new_vecs;
	std::vector<int> new_face_id_0;
	std::vector<int> new_face_id_1;
	std::vector<int> new_face_id_2;

	int vertices_nb = 0;
	for (int i = 0; i < vecs.size(); i++)
	{
		if (lables[i] == 0)
		{
			Vector3d v = vecs[i];
			new_vecs.push_back(v);
			lables[i] = vertices_nb;
			vertices_nb++;
		}
	}

	for (int i = 0; i < face_id_0.size(); i++)
	{
		if (triangles_lables[i] == index)
		{
			new_face_id_0.push_back(lables[face_id_0[i]]);
			new_face_id_1.push_back(lables[face_id_1[i]]);
			new_face_id_2.push_back(lables[face_id_2[i]]);
		}
	}

	CGAL_Output_Obj(path, new_vecs, new_face_id_0, new_face_id_1, new_face_id_2);
}


void CGAL_Output_Off(std::string path, Vector3d1 &vecs, std::vector<int> &face_id_0, std::vector<int> &face_id_1, std::vector<int> &face_id_2)
{
	if (vecs.size() < 3 || face_id_0.size() < 1 || face_id_1.size() < 1 || face_id_2.size() < 1)
	{
		std::cout << "CGAL_Output_Off error: vecs.size() < 3 || face_id_0.size() < 1 || face_id_1.size() < 1 || face_id_2.size() < 1" << std::endl;
		return;
	}

	for (int i = 0; i < face_id_0.size(); i++)
	{
		int index_0 = face_id_0[i];
		int index_1 = face_id_1[i];
		int index_2 = face_id_2[i];

		if (index_0 < 0 || index_0 >= vecs.size() || index_1 < 0 || index_1 >= vecs.size() || index_2 < 0 || index_2 >= vecs.size())
		{
			std::cout << "CGAL_Output_Off error: index_0 < 0 || index_0 >= vecs.size() || index_1 < 0 || index_1 >= vecs.size() || index_2 < 0 || index_2 >= vecs.size()" << std::endl;
			return;
		}
	}
	std::ofstream file(path);
	file << "OFF" << std::endl;
	file << vecs.size() << " " << face_id_0.size() << " 0" << std::endl;
	for (int i = 0; i < vecs.size(); i++)
		file << vecs[i][0] << " " << vecs[i][1] << " " << vecs[i][2] << std::endl;
	for (int i = 0; i < face_id_0.size(); i++)
		file << "3 " << face_id_0[i] << " " << face_id_1[i] << " " << face_id_2[i] << " " << std::endl;
	file.close();
}


void CGAL_Rotation_Obj(std::string path, double angle, Vector3d axis)
{
	Vector3d1 vecs;
	std::vector<int> face_id_0;
	std::vector<int> face_id_1;
	std::vector<int> face_id_2;

	CGAL_3D_Read_Triangle_Mesh(path, vecs, face_id_0, face_id_1, face_id_2);
	for (int i = 0; i < vecs.size(); i++)
	{
		Vector3d v = Functs::RotationAxis(vecs[i], angle, axis);
		vecs[i] = v;
	}
	CGAL_Output_Obj(path, vecs, face_id_0, face_id_1, face_id_2);
}

void CGAL_Export_Path_Segment(std::ofstream &export_file_output, int &export_index,
	std::string s_name, double r, double g, double b, Vector3d start, Vector3d end, double radius)
{
	Vector3d normal = end - start;
	Vector3d base_1 = CGAL_3D_Plane_Base_1(start, normal);
	double length_base_1 = glm::length(base_1);

	base_1[0] = base_1[0] / length_base_1 * radius;
	base_1[1] = base_1[1] / length_base_1 * radius;
	base_1[2] = base_1[2] / length_base_1 * radius;

	Vector3d1 vecs;

	for (int i = 0; i < 4; i++)
	{
		double angle = i * 2 * Math_PI / 4;
		Vector3d v = Functs::RotationAxis(normal + base_1, angle, normal);
		vecs.push_back(v + start);
	}
	for (int i = 0; i < 4; i++)
	{
		vecs.push_back(vecs[i] - normal);
	}

	std::vector<std::vector<int>> faces;

	int face_index_0[4] = { 0, 1, 2, 3 };
	int face_index_1[4] = { 5, 1, 0, 4 };
	int face_index_2[4] = { 4, 0, 3, 7 };
	int face_index_3[4] = { 5, 4, 7, 6 };
	int face_index_4[4] = { 7, 3, 2, 6 };
	int face_index_5[4] = { 6, 2, 1, 5 };

	faces.push_back(std::vector<int>(face_index_0, face_index_0 + 4));
	faces.push_back(std::vector<int>(face_index_1, face_index_1 + 4));
	faces.push_back(std::vector<int>(face_index_2, face_index_2 + 4));
	faces.push_back(std::vector<int>(face_index_3, face_index_3 + 4));
	faces.push_back(std::vector<int>(face_index_4, face_index_4 + 4));
	faces.push_back(std::vector<int>(face_index_5, face_index_5 + 4));

	export_file_output << "g " + s_name << std::endl;

	for (int i = 0; i < vecs.size(); i++)
	{
		//export_file_output << "v " << vecs[i][0] << " " << vecs[i][1] << " " << vecs[i][2] << " " << r << " " << g << " " << b << std::endl;
		export_file_output << "v " << vecs[i][0] << " " << vecs[i][1] << " " << vecs[i][2]  << std::endl;
	}

	for (int i = 0; i < faces.size(); i++)
	{
		export_file_output << "f ";

		for (int j = 0; j < faces[i].size(); j++)
		{
			export_file_output << faces[i][j] + export_index << " ";
		}
		export_file_output << "" << std::endl;
	}

	export_index += 8;
}

void CGAL_Export_Path_Point(std::ofstream &export_file_output, int &export_index,
	std::string s_name, double r, double g, double b, Vector3d point, double radius)
{
	Vector3d1 vecs;
	vecs.push_back(Vector3d(0.5, 0.5, 0.5));
	vecs.push_back(Vector3d(-0.5, 0.5, 0.5));
	vecs.push_back(Vector3d(-0.5, 0.5, -0.5));
	vecs.push_back(Vector3d(0.5, 0.5, -0.5));

	vecs.push_back(Vector3d(0.5, -0.5, 0.5));
	vecs.push_back(Vector3d(-0.5, -0.5, 0.5));
	vecs.push_back(Vector3d(-0.5, -0.5, -0.5));
	vecs.push_back(Vector3d(0.5, -0.5, -0.5));

	std::vector<std::vector<int>> faces;

	int face_index_0[4] = { 0, 1, 2, 3 };
	int face_index_1[4] = { 5, 1, 0, 4 };
	int face_index_2[4] = { 4, 0, 3, 7 };
	int face_index_3[4] = { 5, 4, 7, 6 };
	int face_index_4[4] = { 7, 3, 2, 6 };
	int face_index_5[4] = { 6, 2, 1, 5 };

	faces.push_back(std::vector<int>(face_index_0, face_index_0 + 4));
	faces.push_back(std::vector<int>(face_index_1, face_index_1 + 4));
	faces.push_back(std::vector<int>(face_index_2, face_index_2 + 4));
	faces.push_back(std::vector<int>(face_index_3, face_index_3 + 4));
	faces.push_back(std::vector<int>(face_index_4, face_index_4 + 4));
	faces.push_back(std::vector<int>(face_index_5, face_index_5 + 4));

	export_file_output << "g " + s_name << std::endl;

	for (int i = 0; i < vecs.size(); i++)
	{
		vecs[i][0] = vecs[i][0] * radius;
		vecs[i][1] = vecs[i][1] * radius;
		vecs[i][2] = vecs[i][2] * radius;

		vecs[i][0] += point[0];
		vecs[i][1] += point[1];
		vecs[i][2] += point[2];

		//export_file_output << "v " << vecs[i][0] << " " << vecs[i][1] << " " << vecs[i][2] << " " << r << " " << g << " " << b << std::endl;
		export_file_output << "v " << vecs[i][0] << " " << vecs[i][1] << " " << vecs[i][2] << std::endl;
	}
	for (int i = 0; i < faces.size(); i++)
	{
		export_file_output << "f ";

		for (int j = faces[i].size() - 1; j >= 0; j--)
		{
			export_file_output << faces[i][j] + export_index << " ";
		}
		export_file_output << "" << std::endl;
	}
	export_index += 8;
}

