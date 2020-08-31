#ifndef CFSCNC_ONCE
#define CFSCNC_ONCE
#pragma once
#include "MathHelper.h"
#include <iostream>
#include <fstream>
#include <cassert>
#include <vector>
#include <string>
#include <algorithm>

namespace cnc {


	struct Part
	{
	public:
		//int node_id;
		Vector3d1 points;

		std::vector<int> connecting_edges_id;
		std::vector<int> connecting_edges_index;

		std::vector<int> related_part_id;
		std::vector<bool> related_part_end_bool;
	};

	struct BaseNode
	{
		Vector3d2 boundary;
		int layer;
		int index;
		bool used;
		std::vector<int> neighbors;
		double plane_d;
		BaseNode()
		{
			used = false;
		};

		void PushNeighbors(int node)
		{
			if (std::find(neighbors.begin(), neighbors.end(), node) == neighbors.end())
			{
				neighbors.emplace_back(node);
			}
		}
	};

	struct Node
	{
	public:
		int id;
		int layer_id;
		bool has_connected;
		string connecting_label;
		//distinguished by differenct connecting method
		//"spiral": smooth fermat spiral
		//"bruteforce": bruteforce method

		int degree;
		//tree degree
		//"leaf": degree =1
		//"trunk": degree >2 
		//other: degree=2

		Vector3d1 points;
		Vector3d1 surface_normal;
		std::vector<int> edges_id;//related to edge id
		std::vector<int> edges_index;//related to edge index
		//each element in "edges_id" and "edges_index" is one-to-one correspondence

		std::vector<std::vector<double>> cutting_parts;
		//cutting parts along the node points
		//Each element has two double, referring to one cutting "short" part
		std::vector<int> cutting_parts_order;

		int GetIndex(int order)
		{
			int index = -1;
			for (int i = 0; i < cutting_parts_order.size(); i++)
			{
				if (cutting_parts_order[i] == order)
				{
					index = i;
					break;
				}
			}
			return index;
		}

		void GetCuttingPartsOrders()
		{
			if (cutting_parts.size() <= 0)
				return;

			for (int i = 0; i < edges_id.size(); i++)
			{
				int order = 0;
				for (int j = 0; j < edges_id.size(); j++)
				{
					if (i != j && ((cutting_parts[i][0] + cutting_parts[i][1])>(cutting_parts[j][0] + cutting_parts[j][1])))
					{
						order++;
					}
				}
				cutting_parts_order.push_back(order);
			}
		}
	};

	struct Edge
	{
	public:
		int edge_id;
		std::vector<int> node_id;//two node id of current edge 
		//std::vector<std::vector<double>> sharing_parts;
		Vector3d3 sharing_parts_v;//the possible connecting part between the two nodes' contour

		//std::vector<std::vector<double>> cutting_parts;
		Vector3d2 cutting_parts_v;//????
		Vector3d2 connecting_segments;
		//connecting segments between two nodes of this edge
		//Once time, two segments will be pushed into "connecting_segments"
	};

	class CFSCNC
	{
	public:

		bool corrent_direction;
		double toolpath_size;
		double max_scallop;
		double drill_radius;
		double chord_error;
		double shortest_segments_length;

		int sub_div;

		std::string input_path;

		Vector3d1 debug;
		Vector3d1 debug1;
		std::vector<bool> debug0;

		std::string cls_path;
		std::string scallop_base_mesh;
		std::string scallop_simulated_mesh;
		std::string scallop_map_path;
		std::string measure_path_output;

		std::string nc_path;
		Vector3d1 sharp_turns;

		std::vector<Node> contour_nodes;
		std::vector<Edge> contour_edges;

		Vector3d2 offsets;
		Vector3d2 single_paths;
		Vector3d1 single_path;
		std::vector<bool> single_path_fixed_label;

		//

		bool whether_using_heat;

		Vector3d1 zigzag_final_path;

		Vector3d1 single_final_path;
		Vector3d1 single_final_path_0;
		std::vector<double> single_final_path_v;

		Vector3d1 single_final_CL_path;
		Vector3d1 single_final_path_normal;
		Vector3d1 single_final_path_RMDF_normal;

		/**********************************************************************/
		//Vector3d1 surface_vertices;
		//Vector3d1 surface_vertices_normals;
		/**********************************************************************/

		std::vector<string> offsets_name;
		std::vector<string> single_path_name;
		std::vector<string> single_final_path_name;
		std::vector<string> single_final_path_point_name;

		Vector3d2 single_path_parts;

		Vector3d1 anchors;

		std::vector<Part> parts;
		std::vector<int> merge_edge;
		std::vector<int> removed_merge_edge;
		std::vector<bool> used;

		bool fermat_success;

		int rerouting_points_nb;

		std::vector<std::vector<int>> objs;
		int cfs_index;

		
		CFSCNC();
		void ClearAll();

		//rough machining


		void OffsetsBasedCFSLinking(const std::string &path, const Vector3d2 &boundaries, const bool &re_running, const bool magic_number=true);
		void OffsetsBasedCFSLinking_Framework(std::string path, Vector3d2 &boundaries, bool re_running);
		void OffsetsGeneration(std::string path, bool re_running);
		void HeatGradient(std::string path, bool re_running);

		void ShowPathResults(std::string path);
		void GenerateFermatSpirals(std::vector<int> &mst, int start_nodel_id);
		bool GenerateFermatSpirals_HybridMethod(std::vector<int> &mst);
		/***************************************************************************************************************** /
		/*Generate fermat spirals
		//MinimalSpanningTree() must be called before calling this function.
		//The solution is based on a depth-first searching.
		//The whole process starts from a leaf node.
		//mst: input tree
		//start_node_id: the first node to be searched
		*****************************************************************************************************************/

		void GenerateSmoothFermatSpirals(std::vector<int> &mst, int start_nodel_id);

		void InputDataStructure(Vector3d3 &offsetses, 
			Vector3d3 &offset_graph_sharing_parts, std::vector<int> &mst, std::vector<int> &offset_graph);

		void FermatSpiral(std::vector<int> sequence);
		//For a local spirallable iso-contours, generate a Fermat spiral with a brute-force pattern
		//sequence: the input sequence node ids of the spirallable iso-contours

		void SmoothFermatSpiral(std::vector<int> sequence);

		bool SmoothFermatSpiral_2(std::vector<int> sequence);

		bool GetEdgeid(int index_0, int index_1, int &edge_id, int &source_index, int &target_index);
		bool GetEdgeid(int index_0, int index_1, int &edge_id);
		bool ConnectTwoContours(int index_0, int index_1);

		bool CheckCutValid(int edge_index, int index, double source, double target);
		bool CheckCutValid(int edge_index, int index, Vector3d1 part);

		void ShowZigzag(std::string  path);
		void ShowZigzag1(std::string  path);
		void ShowUGCLS(std::string  path);
		double ComputeGapFromScallop(double surface_curvature, double R_cutter, double scallop);
		void ShowScallopComputation();
		void ShowZigzag3(std::string  path);
		void ShowSimulationPath();
		void ShowNCPath(std::string  path);
		double InterferenceOneRayDistance(Vector3d v, Vector3d n);
		double ClosedPoint(Vector3d v);



		void ScallopHeightSimulation(Vector3d2 &cc_pathes, std::string base_mesh, std::string output_path);
		void ScallopHeightSimulation(Vector3d1 &cc_path, std::string base_mesh, std::string output_path);

		bool CheckValidDirection(Vector3d source_p_0, Vector3d2 delta_ps_0, Vector3d source_p_1, Vector3d2 delta_ps_1);

		Vector3d2 GetNodeCuttingParts(int node_id);
		std::vector<std::vector<double>> GetNodeCuttingPartsPar(int node_id);

		void GetOneSinglePath();
		bool GetOneSinglePath_0();
		bool ConnectAllParts(int start_part, int start_part_end);


		void Output_Boundary(std::string in_path, std::string out_path);
		void Output_Heat_Loop_Curvature_Factor(std::string in_path, std::string out_path);
		void Extract_ISOContours_From_Heat(std::string obj_path, std::string dist_path, std::vector<double> &dists, Vector3d2 &offsets);
		void Extract_ISO_Scallop_Contours(std::string obj_path, std::string  full_obj_path, Vector3d2  &cutter_locations);
		void Output_Path(std::string path);
		void Output_Path_with_Normal(std::string path);
		void Output_Obj_Cur_Normals(std::string path);
		void Output_Path1(std::string path);
		void Output_Path2(std::string path);


		void OutputStripNGC(const std::string path, const Vector3d1 &offsets, bool name_b = true);
		void OutputStrip(const std::string path, const Vector3d1 &offsets, bool name_b = true);



		void Load_Path(std::string path);
		void Load_RMDF_Path(std::string path);
		void Load_Zigzag_Path(std::string path);


		void LoadOffsetsFiles(std::vector<string> &offsets_files, std::string path);

		bool Load_Final_Path(std::string path);
		void LoadContour(std::string path, Vector3d2 &offsets, std::vector<double> &dists);
		void LoadContours(Vector3d3 &offsetses, Vector3d2 &offsets, bool re_running);
		void UnifiedContourDirection(double toolpath_size, Vector3d3 &offsetses, Vector3d2 &offsets);
	
		double GeodesicDistance(Vector3d source, Vector3d target);

		void GenerateAdaptiveCurvatureMesh();

		std::vector<double> Heat_to_Gradient_0(std::string path, std::string obj_path, std::string full_obj_path, std::string heat_path, std::string scale_path);

		double  CFSCNC::ComputeScale(double cur);

		void LoadBoundary(std::string path, Vector3d2 &input_boundaries);

		void BuildOffsets(const Vector3d2 &input_boundaries, Vector3d3 &offsetses, Vector3d2 &offsets);

	};

}
#endif