#include "stdafx.h"
#include "CFSCNC.h"
#include "Circuit.h"
#include "cgalpackage.h"

using namespace hpcg;

namespace cnc {


	void CFSCNC::HeatGradient(std::string path, bool re_running)
	{
	
	}


	void CFSCNC::OffsetsGeneration(std::string path, bool re_running)
	{
		input_path = path;

		if (!std::ifstream(path + "\\path\\" + IntString(cfs_index) + "_split.obj", std::ios::in) || re_running)
		{
			string com_split = "";
			com_split += "D:\\CNCProduction\\Release\\geodesic_isoline_usage_shiqing\\Split.exe";
			com_split += " ";
			com_split += path + "\\" + IntString(cfs_index) + ".obj";
			com_split += " ";
			com_split += "0.2";
			com_split += " ";
			com_split += path + "\\path\\" + IntString(cfs_index) + "_split.obj";
			system(com_split.c_str());
		}

		if (!std::ifstream(path + "\\path\\" + IntString(cfs_index) + ".offsets_path", std::ios::in) || re_running)
		{
			if (whether_using_heat)
			{
				if (!std::ifstream(path + "\\path\\" + IntString(cfs_index) + "_heat_sub.obj", std::ios::in) || re_running)
					CGAL_Mesh_Loop_Subdivision_Own_Version(path + "\\" + IntString(cfs_index) + ".obj", sub_div, path + "\\path\\" + IntString(cfs_index) + "_heat_sub.obj");

				if (!std::ifstream(path + "\\path\\" + IntString(cfs_index) + "_heat_sub_split.obj", std::ios::in) || re_running)
				{
					string com_split = "";
					com_split += "D:\\CNCProduction\\Release\\geodesic_isoline_usage_shiqing\\Split.exe";
					com_split += " ";
					com_split += path + "\\path\\" + IntString(cfs_index) + "_heat_sub.obj";
					com_split += " ";
					com_split += "0.4";
					com_split += " ";
					com_split += path + "\\path\\" + IntString(cfs_index) + "_heat_sub_split.obj";
					system(com_split.c_str());
				}

				Output_Boundary(path + "\\path\\" + IntString(cfs_index) + "_heat_sub_split.obj", path + "\\path\\" + IntString(cfs_index) + "_heat.source");
				//Output_Heat_Loop_Curvature_Factor(path + "\\path\\" + IntString(cfs_index) + "_heat_sub_split.obj", path + "\\path\\" + IntString(cfs_index) + "_heat_factor.source");
			}
			else
			{
				std::vector<string> offsets_files;

				for (int i = 0; i < objs[cfs_index].size(); i++)
				{
					std::cout << "Compute offsets: " << IntString(cfs_index) << " " << IntString(objs[cfs_index][i]) << std::endl;
					if (!std::ifstream(path + "\\path\\" + IntString(cfs_index) + "_" + IntString(objs[cfs_index][i]) + "_sub.obj", std::ios::in) || re_running)
					{
						string com_sub = "";
						com_sub += "D:\\CNCProduction\\Release\\geodesic_isoline_usage_shiqing\\Subdivide.exe";
						com_sub += " ";
						com_sub += path + "\\" + IntString(cfs_index) + "_" + IntString(objs[cfs_index][i]) + ".obj";
						com_sub += " ";
						com_sub += DoubleString(0.1);
						com_sub += " ";
						com_sub += path + "\\path\\" + IntString(cfs_index) + "_" + IntString(objs[cfs_index][i]) + "_sub.obj";
						system(com_sub.c_str());
					}

					if (!std::ifstream(path + "\\path\\" + IntString(cfs_index) + "_" + IntString(objs[cfs_index][i]) + "_sub_split.obj", std::ios::in) || re_running)
					{
						string com_split = "";
						com_split += "D:\\CNCProduction\\Release\\geodesic_isoline_usage_shiqing\\Split.exe";
						com_split += " ";
						com_split += path + "\\path\\" + IntString(cfs_index) + "_" + IntString(objs[cfs_index][i]) + "_sub.obj";
						com_split += " ";
						com_split += "0.1";
						com_split += " ";
						com_split += path + "\\path\\" + IntString(cfs_index) + "_" + IntString(objs[cfs_index][i]) + "_sub_split.obj";
						system(com_split.c_str());
					}

					if (!std::ifstream(path + "\\path\\" + IntString(cfs_index) + "_" + IntString(objs[cfs_index][i]) + ".offsets", std::ios::in) || re_running)
					{
						string com_isoline = "";
						com_isoline += "D:\\CNCProduction\\Release\\geodesic_isoline_usage_shiqing\\isoline.exe";
						com_isoline += " ";
						com_isoline += path + "\\path\\" + IntString(cfs_index) + "_" + IntString(objs[cfs_index][i]) + "_sub_split.obj";
						//com_isoline += path +"\\"+ IntString(cfs_index) + "_" + IntString(i) + ".obj";
						com_isoline += " ";
						com_isoline += DoubleString(max_scallop);
						com_isoline += " ";
						com_isoline += DoubleString(drill_radius);
						com_isoline += " ";
						com_isoline += path + "\\path\\" + IntString(cfs_index) + "_" + IntString(objs[cfs_index][i]) + ".offsets";
						system(com_isoline.c_str());
					}
					/***********************************************************************************************/
					offsets_files.push_back(path + "\\path\\" + IntString(cfs_index) + "_" + IntString(objs[cfs_index][i]) + ".offsets");
				}

				Circuit::OutputOffsetsFiles(offsets_files, path + "\\path\\" + IntString(cfs_index) + ".offsets_path");
			}
		}
	}

	double  CFSCNC::ComputeScale(double cur)
	{
		double flat_width = ComputeGapFromScallop(0.0, 2.0, max_scallop);

		double w = ComputeGapFromScallop(cur, 2.0, max_scallop);

		/*	if (w > flat_width*2.0) w = flat_width*2.0;
		if (w < 0.3) w = 0.3;*/
		return w * w / flat_width / flat_width;
	}

	//Get all cutting parts of node id
	//node_id: node id
	Vector3d2 CFSCNC::GetNodeCuttingParts(int node_id)
	{
		Vector3d2 cutting_parts;
		for (int i = 0; i < contour_nodes[node_id].edges_id.size(); i++)
		{
			int edge_id = contour_nodes[node_id].edges_id[i];
			int edge_index = contour_nodes[node_id].edges_index[i];

			if (contour_edges[edge_id].cutting_parts_v.size()>0)
				cutting_parts.push_back(contour_edges[edge_id].cutting_parts_v[edge_index]);
		}
		return 	cutting_parts;
	}

	bool CFSCNC::GetEdgeid(int index_0, int index_1, int &edge_id)
	{
		for (int i = 0; i < contour_edges.size(); i++)
		{
			if (contour_edges[i].node_id[0] == index_0 && contour_edges[i].node_id[1] == index_1)
			{
				edge_id = i;
				return true;
			}
			if (contour_edges[i].node_id[0] == index_1 && contour_edges[i].node_id[1] == index_0)
			{
				edge_id = i;
				return true;
			}
		}
		return false;
	}

	bool CFSCNC::GetEdgeid(int index_0, int index_1, int &edge_id, int &source_index, int &target_index)
	{
		for (int i = 0; i < contour_edges.size(); i++)
		{
			if (contour_edges[i].node_id[0] == index_0 && contour_edges[i].node_id[1] == index_1)
			{
				edge_id = i;
				source_index = 0;
				target_index = 1;
				return true;
			}
			if (contour_edges[i].node_id[0] == index_1 && contour_edges[i].node_id[1] == index_0)
			{
				edge_id = i;
				source_index = 1;
				target_index = 0;
				return true;
			}
		}
		return false;
	}

}