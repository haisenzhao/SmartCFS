#ifndef ROUGHCNC_ONCE
#define ROUGHCNC_ONCE
#pragma once
#include "MathHelper.h"
#include <iostream>
#include <fstream>
#include <cassert>
#include <vector>
#include <string>
#include <algorithm>

namespace cnc {

	class RoughMachining
	{
	public:

		bool output_debug_offset;
		bool output_rerun;

		double area_threshold;
		double length_threshold;
		double removed_length_threshold;
		double layer_connecting_threshold;

		int volume_index;

		CFSCNC cfs;
		RoughMachining() :output_debug_offset(true), output_rerun(false)
		{
		};

		RoughMachining(std::string path, bool re_running = true,
			bool output_debug = false, double set_toolpath_size = -1.0, 
			double set_layer_thickness = -1.0, 
			bool segment_or_tool_path = true, int volume_index_ = -1);

		Vector3d1 GenerateToolpathForOneVolume(const int &segment_index, const std::string &path,
			const Vector3d3 &rough_boundaries, const std::vector<double> &plane_d, const bool &re_running, const glm::dmat4 &rm);
		void BuildRelationShip(const std::string &path, Vector3d3 &layer_boundaries, const std::vector<double> &plane_d, const glm::dmat4 &rm, const double toolpath_size,
			std::vector<BaseNode> &base_nodes, std::vector<std::vector<int>> &layer_nodes, std::vector<int> &mst);
		void SelectSequeueFromTree(const std::string &path, const Vector3d3 &layer_boundaries, std::vector<BaseNode> &base_nodes,
			const std::vector<std::vector<int>> layer_nodes, std::vector<std::vector<int>> &sequences);
		void SegmentMeshReconstruct(const std::string &path, const std::vector<double> &plane_d, const Vector3d &plane_n, const double &layer_thichness,
			const std::vector<BaseNode> &base_nodes, const std::vector<std::vector<int>> &sequences, const glm::dmat4 &rm);
		void SegmentMeshReconstruct2(const std::string &path, const std::vector<double> &plane_d, const Vector3d &plane_n, const double &layer_thichness,
			const std::vector<BaseNode> &base_nodes, const std::vector<std::vector<int>> &sequences, const glm::dmat4 &rm);

		void OutputCLS(const std::string folder_name, const std::string output_path, const Vector3d1 &rough_machining_path, const Vector3d &plane_n);
		void OutputFinalNGC(const std::string output_path, const Vector3d1 &rough_machining_path, const Vector3d &plane_n);

		void MeshSlicer(const std::string &path, const std::string off_path, const Vector3d &plane_n, const std::vector<double> &plane_d, const glm::dmat4 &rm,
			const double toolpath_size, Vector3d3 &rough_boundaries);
	};

}
#endif