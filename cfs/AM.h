#ifndef AM_ONCE
#define AM_ONCE
#pragma once
#include "MathHelper.h"
#include <iostream>
#include <fstream>
#include <cassert>
#include <vector>
#include <string>
#include <algorithm>

namespace cnc {

	class AMVolume
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
		AMVolume(std::string path, bool re_running = true,
			bool output_debug = false, double set_toolpath_size = -1.0,
			double set_layer_thickness = -1.0, bool segment_or_tool_path = true);
		
		void MeshSlicer(const std::string &path, const std::string off_path, const Vector3d &plane_n, const double &layer_thichness,
			const std::vector<double> &plane_d, const glm::dmat4 &rm,
			const double toolpath_size, Vector3d3 &rough_boundaries);

		void BuildRelationShip(const std::string &path, Vector3d3 &layer_boundaries, const std::vector<double> &plane_d, const glm::dmat4 &rm, const double toolpath_size,
			std::vector<BaseNode> &base_nodes, std::vector<std::vector<int>> &layer_nodes, std::vector<int> &mst);
	};

}
#endif