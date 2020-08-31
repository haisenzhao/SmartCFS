#include "stdafx.h"
#include "CFSCNC.h"
#include <Circuit.h>
#include "Tree.h"
#include "RoughMachining.h"

namespace cnc {

	Vector3d1 Helix(const Vector3d &start, const Vector3d &end, int helix_nb=48)
	{
		return Vector3d1();
		double radius = 0.2;
		auto normal = end - start;
		double total_depth = Math::GetLength(normal);

		Vector3d base_1 = CGAL_3D_Plane_Base_1(start, normal);
		double length_base_1 = glm::length(base_1);

		base_1[0] = base_1[0] / length_base_1 * radius;
		base_1[1] = base_1[1] / length_base_1 * radius;
		base_1[2] = base_1[2] / length_base_1 * radius;

		Vector3d1 vecs(1, base_1 + start);

		for (int i = 0; i < helix_nb; i++)
		{
			double depth = (double)i / (double)47 * total_depth;
			auto nnn = normal;
			auto aaa = Math::SetVectorLength(nnn, depth);

			double angle = i *  Math::Math_PI / 4;
			Vector3d v = Math::RotationAxis(base_1 + aaa, angle, normal);
			vecs.push_back(v + start);
		}

		return vecs;
	};

	//load setting
	void LoadSetting(const std::string &path, Vector3d &plane_n, Vector3d &origin_n,
		std::vector<double> &plane_d, double &toolpath_size, double &layer_thichness, 
		glm::dmat4 &rm, std::string &obj_file, std::string &off_file, double set_toolpath_size, double set_layer_thickness)
	{
		std::ifstream file(path + "cutter_direction.txt", std::ios::in);
		double minimal_d;
		double maximal_d;

		file >> toolpath_size;
		file >> plane_n.x >> plane_n.y >> plane_n.z;
		file >> layer_thichness;
		file >> obj_file;

		if (set_toolpath_size > 0.0&&set_layer_thickness > 0.0)
		{
			toolpath_size = set_toolpath_size;
			layer_thichness = set_layer_thickness;
		}

		plane_n = -plane_n;
		origin_n = plane_n;
		rm = RotationMatrix(Vector3d(0.0, 0.0, 1.0), plane_n);
		plane_n = Vector3d(0.0, 0.0, 1.0);

		Vector3d1 vecs;
		std::vector<int> face_id_0, face_id_1, face_id_2;
		CGAL_3D_Read_Triangle_Mesh(path + obj_file, vecs, face_id_0, face_id_1, face_id_2);

		vecs = PosApplyM(vecs, rm);

		CGAL_Output_Off(path + "transform.off", vecs, face_id_0, face_id_1, face_id_2);
		CGAL_Output_Obj(path + "transform.obj", vecs, face_id_0, face_id_1, face_id_2);

		off_file = "transform.off";

		Vector3d min_corner, max_corner;
		CGAL_3D_Mesh_Boundingbox(vecs, min_corner, max_corner);
		minimal_d = min_corner[2];
		maximal_d = max_corner[2];

		for (double d = minimal_d; d < maximal_d; d = d + layer_thichness)
			plane_d.emplace_back(d);

		//plane_d.emplace_back(plane_d.back() + layer_thichness);
		//plane_d.emplace_back(plane_d.back() + layer_thichness);

		file.clear();
		file.close();
	};
	
	RoughMachining::RoughMachining(std::string path, bool re_running, bool output_debug,
		double set_toolpath_size, double set_layer_thickness, bool segment_or_tool_path, int volume_index_)
	{

		volume_index = volume_index_;
		//threshold
		area_threshold = 3.0;
		length_threshold = 10.0;
		removed_length_threshold = 5.0;
		layer_connecting_threshold = 0.8;


		output_debug_offset = output_debug;
		output_rerun = re_running;

		auto RoughPathGeneration = [&](const std::vector<double> &plane_d, const Vector3d &plane_n, const double &layer_thichness,
			const std::vector<BaseNode> &base_nodes, const std::vector<std::vector<int>> &sequences, const glm::dmat4 &rm)
		{
			double extracting_back_distance = 20.0;

			DWORD start_time = GetTickCount();

			Vector3d1 rough_machining_path;
			for (int i = 0; i < sequences.size(); i++)
			{
				//if (i != 1)continue;
				std::cerr << "RoughMachining Iteration: " << i << std::endl;

				Vector3d3 rough_boundaries_0;
				std::vector<double> plane_d_0;

				for (int seq_index = 0; seq_index < sequences[i].size();seq_index++)
				{
					auto seq = sequences[i][seq_index];

					rough_boundaries_0.emplace_back(base_nodes[seq].boundary);
					plane_d_0.emplace_back(base_nodes[seq].plane_d);
					
					if (seq_index == sequences[i].size() - 1)
					{
						auto boundary_1 = base_nodes[seq].boundary;
						for (auto& a : boundary_1)for (auto&b : a)b[2] = b[2] + layer_thichness;
						rough_boundaries_0.emplace_back(boundary_1);
						plane_d_0.emplace_back(base_nodes[seq].plane_d + layer_thichness);

						auto boundary_2 = base_nodes[seq].boundary;
						for (auto& a : boundary_2)for (auto&b : a)b[2] = b[2] + 2*layer_thichness;
						rough_boundaries_0.emplace_back(boundary_2);
						plane_d_0.emplace_back(base_nodes[seq].plane_d + 2*layer_thichness);

					}
				}

				auto segment_path = GenerateToolpathForOneVolume(i, path, rough_boundaries_0, plane_d_0, re_running,rm);

				if (i == 0)
				{
					//helix
					rough_machining_path.emplace_back(segment_path.front());
					rough_machining_path.back()[2] = plane_d_0[0] + 1.0*(-plane_n.z);
					auto helix = Helix(rough_machining_path.back(), segment_path.front());
					for (auto h : helix)
						rough_machining_path.emplace_back(h);

					for (auto p : segment_path)
						rough_machining_path.emplace_back(p);
					rough_machining_path.emplace_back(segment_path.back());
					rough_machining_path.back()[2] = plane_d[0] + extracting_back_distance*(-plane_n.z);
				}
				else
				{
					rough_machining_path.emplace_back(segment_path.front());
					rough_machining_path.back()[2] = plane_d[0] + extracting_back_distance*(-plane_n.z);

					//helix
					rough_machining_path.emplace_back(segment_path.front());
					rough_machining_path.back()[2] = plane_d_0[0] + 1.0*(-plane_n.z);
					auto helix = Helix(rough_machining_path.back(), segment_path.front());
					for (auto h : helix)
						rough_machining_path.emplace_back(h);

					for (auto p : segment_path)
						rough_machining_path.emplace_back(p);
					rough_machining_path.emplace_back(segment_path.back());
					rough_machining_path.back()[2] = plane_d[0] + extracting_back_distance*(-plane_n.z);
				}
			}

			DWORD end_time = GetTickCount();
			std::cout << "[TIME] ABCD TIME: " << (end_time - start_time) / 1000.0 << std::endl;


			return rough_machining_path;
		};

		//clear folder
		std::cerr << "#Clear working folder..." << std::endl;
		Math::ClearFolder(path + "boundary\\");
		Math::ClearFolder(path + "path\\");
		Math::ClearFolder(path + "debug\\");
		//if (re_running)Math::ClearFolder(path + "output\\");

		//load setting
		std::cerr << "# Load setting..." << std::endl;
		std::vector<double> plane_d;
		Vector3d plane_n;
		double layer_thichness;
		std::string obj_file, off_file;
		glm::dmat4 rm;
		Vector3d origin_n;
	
		LoadSetting(path, plane_n, origin_n, plane_d, cfs.toolpath_size, layer_thichness, rm, obj_file, off_file, set_toolpath_size, set_layer_thickness);
		auto strings = Split(path,"\\");
		std::string folder_name = strings[strings.size()-2];

		//mesh slicer
		std::cerr << "# Mesh slicer..." << std::endl;
		Vector3d3 rough_boundaries;
		MeshSlicer(path, path + off_file, plane_n, layer_thichness, plane_d, rm, cfs.toolpath_size, rough_boundaries);
		int rough_b_nb = 0;
		for (auto & boud : rough_boundaries) rough_b_nb += boud.size();
		if (rough_b_nb == 0)
		{
			std::cerr << "if (rough_b_nb == 0)" << std::endl;
			system("pause");
		}

		if (output_debug_offset)
		{
			std::cerr << "#Output offsets..." << std::endl;
			for (int i = 0; i < rough_boundaries.size(); i++)
			{
				auto a = PosApplyM(rough_boundaries[i], glm::inverse(rm));
				Circuit::OutputOffsets(path + "boundary\\boudary_original_" + std::to_string(i) + ".obj", a, "boudary_original_" + std::to_string(i));
			}
		}


		//build relationship
		std::cerr << "#build laryer relationship..." << std::endl;
		std::vector<BaseNode> base_nodes;
		std::vector<std::vector<int>> sequences;
		std::vector<std::vector<int>> layer_nodes;
		std::vector<int> mst;
		BuildRelationShip(path, rough_boundaries, plane_d, rm, cfs.toolpath_size, base_nodes, layer_nodes, mst);

		if (base_nodes.empty())
		{
			std::cerr << "if (base_nodes.empty())" << std::endl;
			system("pause");
		}

		SelectSequeueFromTree(path, rough_boundaries, base_nodes,layer_nodes, sequences);

		for (int i = 0; i < sequences.size(); i++)
		{
			std::cerr << i << " : ";
			for (int j = 0; j < sequences[i].size(); j++) std::cerr << sequences[i][j] << " ";
			std::cerr << "\n";
		}

		if (false)
		if (output_debug_offset)
		{
			for (int seg_index = 0; seg_index < sequences.size(); seg_index++)
			{
				for (int i = 0; i < sequences[seg_index].size(); i++)
				{
					auto  node_index = sequences[seg_index][i];
					auto name = "layer_boundary_" + std::to_string(seg_index) + "_" + std::to_string(i) + "_" + IntString(node_index);
					auto a = PosApplyM(base_nodes[node_index].boundary, glm::inverse(rm));
					Circuit::OutputOffsets(path + "boundary\\" + name + ".obj", a, name);
				}
			}
		}
		
		if (segment_or_tool_path)
		{
			//output reconstruction mesh
			//std::cerr << "#output reconstruction mesh..." << std::endl;
			//SegmentMeshReconstruct2(path, plane_d, plane_n, layer_thichness, base_nodes, sequences, rm);
			SegmentMeshReconstruct(path, plane_d, plane_n, layer_thichness, base_nodes, sequences,rm);
			
			//return;
		}
		else
		{
			std::cerr << "#Generate rough machining tool path..." << std::endl;

			Vector3d1 rough_machining_path;

			if (re_running || !LoadVectors(path + "output\\rough_machining_path.txt", rough_machining_path))
			{
				//sequence order
				//TODO

				rough_machining_path = RoughPathGeneration(plane_d, plane_n, layer_thichness, base_nodes, sequences, rm);
				OutputVectors(path + "output\\rough_machining_path.txt", rough_machining_path);
			}

			std::cerr << "#Output strips..." << std::endl;
			rough_machining_path = PosApplyM(rough_machining_path, glm::inverse(rm));

			if (output_debug) Circuit::OutputStrip(path + "output\\rough_path.obj", rough_machining_path);
			std::cerr << "#Output g-code..." << std::endl;

			OutputFinalNGC(path + "output\\first.ngc", rough_machining_path, origin_n);

			//std::string folder_name
			OutputCLS(folder_name, path + folder_name + ".cls", rough_machining_path, origin_n);

			double length = Strip::GetTotalLength(rough_machining_path);
			std::cerr << "Total length: " << length << std::endl;

		}

		std::cerr << "#over" << std::endl;
	}

	void OutputNodeTriangles(std::string file_path, Vector3d2 &cut_offsets, const glm::dmat4 &rm)
	{
		std::ofstream file(file_path);

		auto points_2d = Math::Vector3d2d(cut_offsets);
		auto triangles = CGAL_2D_Polygon_Triangulation(points_2d);

		for (int i = 0; i<cut_offsets.size(); i++)
		{
			for (int j = 0; j < cut_offsets[i].size(); j++)
			{
				auto v = PosApplyM(cut_offsets[i][j], glm::inverse(rm));
				file << "v " << v[0] << " " << v[1] << " " << v[2] << std::endl;
			}
		}
		int nb = 1;
		for (int i = 0; i < triangles.size(); i++)
		{
			int index_0 = triangles[i][0] + nb;
			int index_1 = triangles[i][1] + nb;
			int index_2 = triangles[i][2] + nb;
			file << "f " << index_0 << " " << index_1 << " " << index_2 << std::endl;
		}


		file.close();
	}

	void OutputPolygon(std::string file_path, Vector3d2 &cut_offsets, const glm::dmat4 &rm, const double offset=-0.4)
	{
		std::ofstream file(file_path);
		int nb = 1;
		for (int iter = 0; iter < cut_offsets.size(); iter++)
		{
			auto points_2d = Math::Vector3d2d(cut_offsets[iter]);
			std::vector<std::vector<Vector2d>> offsets_xys;
			
			CGAL_2D_Polygon_One_Offsets(points_2d, offset, offsets_xys);

			auto points_3d = Math::Vector2d3d(offsets_xys[0], cut_offsets[iter].front()[2]);
			auto triangles = CGAL_2D_Polygon_Triangulation(Math::Vector3d2d(points_3d));

			for (int i = 0; i<points_3d.size(); i++)
			{
				auto v = PosApplyM(points_3d[i], glm::inverse(rm));
				file << "v " <<v[0] << " " << v[1] << " " << v[2] << std::endl;
				//file << "v " << points_3d[i][0] << " " << points_3d[i][1] << " " << points_3d[i][2] << std::endl;
			}
			for (int i = 0; i < triangles.size(); i++)
			{
				int index_0 = triangles[i][0] + nb;
				int index_1 = triangles[i][1] + nb;
				int index_2 = triangles[i][2] + nb;
				file << "f " << index_0 << " " << index_1 << " " << index_2 << std::endl;
			}
			nb += points_3d.size();
		}
		file.close();
	};

	void RoughMachining::MeshSlicer(const std::string &path, const std::string off_path,
		const Vector3d &plane_n, const double &layer_thichness, const std::vector<double> &plane_d, const glm::dmat4 &rm, const double toolpath_size, Vector3d3 &rough_boundaries)
	{
		Vector3d3 removed_offsets;
		
		if (output_rerun || !LoadVectors(path + "output\\mesh_slicer.txt", rough_boundaries) || !LoadVectors(path + "output\\mesh_slicer_removed_offsets.txt", removed_offsets))
		{
			Vector3d plane_p(0.0, 0.0, 0.0);

			Vector3d3 original_offsetses;
			Vector3d2 original_offsets;
			CGAL_Slicer_Mesh(off_path, Vector3d(0.0, 0.0, 1.0), plane_d, original_offsetses, original_offsets);

			if (original_offsetses.empty() || original_offsets.empty())
			{
				std::cerr << "if (original_offsetses.empty() || original_offsets.empty())" << std::endl;
				system("pause");
			}

			//extend the final layer
			//original_offsetses[original_offsetses.size() - 2] = original_offsetses[original_offsetses.size() - 3];
			//original_offsetses[original_offsetses.size() - 1] = original_offsetses[original_offsetses.size() - 2];
			//for (auto&offsets : original_offsetses[original_offsetses.size() - 2]) for (auto&offset : offsets)offset[2] = plane_d[plane_d.size() - 2];
			//for (auto&offsets : original_offsetses[original_offsetses.size() - 1]) for (auto&offset : offsets)offset[2] = plane_d[plane_d.size() - 1];

			if (output_debug_offset)
			{
				std::cerr << "#Output offsets..." << std::endl;
				for (int i = 0; i < original_offsetses.size(); i++)
				{
					auto a = PosApplyM(original_offsetses[i], glm::inverse(rm));
					Circuit::OutputOffsets(path + "boundary\\boudary_original_" + std::to_string(i) + ".obj", a, "boudary_original_" + std::to_string(i));
				}
			}

			for (int i = 0; i < original_offsetses.size(); i++)
			{
				Vector3d2 layer_offsets;
				Vector3d2 layer_removed_offsets;
				for (int j = 0; j < original_offsetses[i].size(); j++)
				{
					if (original_offsetses[i][j].size() < 3)continue;

					Vector2d1 points_2d;
					CGAL_3D_Plane_3D_to_2D_Points(plane_p, Vector3d(0.0, 0.0, 1.0), original_offsetses[i][j], points_2d);
					auto clear_points_2d = Circuit::ClearCircuit(points_2d);
				
					if (clear_points_2d.size() < 3)continue;

					auto length = Circuit::GetTotalLength(Math::Vector2d3d(clear_points_2d));
					//auto area = CGAL_2D_Polygon_Area(clear_points_2d);
					if (length > length_threshold)
					{
						if (!CGAL_2D_Polygon_Simple(clear_points_2d))
						{
							auto area = CGAL_2D_Polygon_Area(clear_points_2d);


							auto v = PosApplyM(Math::Vector2d3d(clear_points_2d, plane_d[i]), glm::inverse(rm));
							Circuit::OutputOffsets1(path + "debug\\debug.obj", v);

							v = PosApplyM(Math::Vector2d3d(points_2d, plane_d[i]), glm::inverse(rm));
							Circuit::OutputOffsets1(path + "debug\\debug1.obj", v);

							std::cerr << "if (!CGAL_2D_Polygon_Simple(clear_points_2d))" << std::endl;
							system("pause");
							int d = 0;
						}
						else
						{
							auto area = CGAL_2D_Polygon_Area(clear_points_2d);

							if (area > area_threshold)
							{
								layer_offsets.emplace_back(Math::Vector2d3d(clear_points_2d, plane_d[i]));
							}
							else
							{
								if (length>removed_length_threshold) layer_removed_offsets.emplace_back(Math::Vector2d3d(clear_points_2d, plane_d[i]));
							}
				
						}
					}
					else
					{
						if (length>removed_length_threshold) layer_removed_offsets.emplace_back(Math::Vector2d3d(clear_points_2d, plane_d[i]));
					}
				}
				rough_boundaries.emplace_back(layer_offsets);
				removed_offsets.emplace_back(layer_removed_offsets);
			}

			//Circuit::OutputOffsets(path + "boundary\\boudary_removed.obj", removed_offsets,"boudary_removed");
			OutputVectors(path + "output\\mesh_slicer.txt", rough_boundaries);
			OutputVectors(path + "output\\mesh_slicer_removed_offsets.txt", removed_offsets);
		}

		if (output_rerun || !LoadExisting(path + "output\\initial_cut_contours.obj"))
		{
			//remove small mesh
			std::vector<BaseNode> base_nodes;
			std::vector<std::vector<int>> layer_nodes;
			std::vector<int> mst;
			BuildRelationShip(path, removed_offsets, plane_d, rm, toolpath_size,   base_nodes, layer_nodes, mst);

			Vector3d2 cut_offsets;
			std::vector<std::vector<int>>  cc = Tree::ConnectedComponents(base_nodes.size(), mst);
			for (int i = 0; i < cc.size(); i++)
			{
				int out_node = -1;
				int min_layer = 10000000;
				for (int j = 0; j < cc[i].size(); j++)
				{
					if (base_nodes[cc[i][j]].layer < min_layer)
					{
						min_layer = base_nodes[cc[i][j]].layer;
						out_node = cc[i][j];
					}
				}
				for (auto b : base_nodes[out_node].boundary)cut_offsets.emplace_back(b);
			}

			Circuit::OutputOffsets(path + "boundary\\cut_offsets.obj", cut_offsets, "cut_offsets");

			//Vector3d2(1,cut_offsets[8])


			OutputPolygon(path + "output\\initial_cut_contours.obj", cut_offsets,rm);
			//CGAL_Cut_Surface_by_Multi_Boundaries(Vector3d2(1, cut_offsets[18]), base_nodes.front().boundary.front().front(), path + "volume.off", path + "output\\initial_cut.obj");
			//CGAL_Cut_Surface_by_Multi_Boundaries(cut_offsets, base_nodes.front().boundary.front().front(), path + "volume.off", path + "output\\initial_cut.obj");

			//std::cerr << "Generate initial cut manually..." << std::endl;
			//system("pause");
			Math::ClearFolder(path + "boundary\\");
		}

	

	};


	void RoughMachining::SelectSequeueFromTree(const std::string &path, const Vector3d3 &layer_boundaries, std::vector<BaseNode> &base_nodes,
		const std::vector<std::vector<int>> layer_nodes, std::vector<std::vector<int>> &sequences)
	{
		//decomposition
		for (int layer = 0; layer < layer_boundaries.size(); layer++)
		{
			auto current_nodes = layer_nodes[layer];
			for (auto node : current_nodes)
			{
				if (!base_nodes[node].used)
				{
					base_nodes[node].used = true;
					std::vector<int> sequence(1, node);
					while (true)
					{
						int size = sequence.size();
						for (auto candidate_node : base_nodes[sequence.back()].neighbors)
						{
							if (!base_nodes[candidate_node].used &&
								base_nodes[candidate_node].layer>base_nodes[sequence.back()].layer)
							{
								base_nodes[candidate_node].used = true;
								sequence.emplace_back(candidate_node);
								break;
							}
						}
						if (size == sequence.size())break;
					}
					sequences.emplace_back(sequence);
				}
			}
		};
	}

	void RoughMachining::BuildRelationShip(const std::string &path, Vector3d3 &layer_boundaries, const std::vector<double> &plane_d, const glm::dmat4 &rm, 
		const double toolpath_size,
		std::vector<BaseNode> &base_nodes, std::vector<std::vector<int>> &layer_nodes, std::vector<int> &mst)
	{
		auto OutputBaseNodesSequences = [](std::string path, const std::vector<BaseNode> &base_nodes, std::vector<std::vector<int>> &layer_nodes, std::vector<int> &mst)
		{
			std::ofstream file(path);

			////sequences
			//file << sequences.size() << std::endl;
			//for (int i = 0; i < sequences.size(); i++)
			//{
			//	file << sequences[i].size() << " ";
			//	for (int j = 0; j < sequences[i].size(); j++)
			//		file << sequences[i][j] << " ";
			//	file << "\n";
			//}

			//mst
			file << mst.size() << std::endl;
			for (int i = 0; i < mst.size(); i++)
				file << mst[i] << " ";
			file << "\n";

			//layer nodes
			file << layer_nodes.size() << std::endl;
			for (int i = 0; i < layer_nodes.size(); i++)
			{
				file << layer_nodes[i].size() << " ";
				for (int j = 0; j < layer_nodes[i].size(); j++)
					file << layer_nodes[i][j] << " ";
				file << "\n";
			}

			//BaseNodes
			file << base_nodes.size() << std::endl;
			for (const auto &node : base_nodes)
			{
				file << node.layer << " " << node.index << " " << node.used << " " << node.plane_d << std::endl;
				file << node.neighbors.size() << " ";
				for (const auto &neighbor : node.neighbors)file << neighbor << " "; file << "\n";

				file << node.boundary.size() << std::endl;
				for (int i = 0; i < node.boundary.size(); i++)
				{
					file << node.boundary[i].size() << std::endl;
					for (int j = 0; j < node.boundary[i].size(); j++)
						file << node.boundary[i][j][0] << " " << node.boundary[i][j][1] << " " << node.boundary[i][j][2] << " ";
					file << "\n";
				}
			}

			file.clear();
			file.close();
		};

		auto LoadBaseNodesSequences = [](std::string path, std::vector<BaseNode> &base_nodes, std::vector<std::vector<int>> &layer_nodes, std::vector<int> &mst)
		{
			std::ifstream file(path, std::ios::in);

			if (!file) return false;

			//sequences
	/*		int seq_nb_0;
			file >> seq_nb_0;
			for (int i = 0; i < seq_nb_0; i++)
			{
				int seq_nb_1;
				file >> seq_nb_1;
				std::vector<int> sequence(seq_nb_1, 0);
				for (int j = 0; j < seq_nb_1; j++) file >> sequence[j];
				sequences.emplace_back(sequence);
			}*/

			//mst
			int mst_nb;
			file >> mst_nb;
			for (int i = 0; i < mst_nb; i++)
			{
				mst.emplace_back(0);
				file >> mst.back();
			}

			//layer nodes
			int seq_nb_0;
			file >> seq_nb_0;
			for (int i = 0; i < seq_nb_0; i++)
			{
				int seq_nb_1;
				file >> seq_nb_1;
				std::vector<int> sequence(seq_nb_1, 0);
				for (int j = 0; j < seq_nb_1; j++) file >> sequence[j];
				layer_nodes.emplace_back(sequence);
			}

			//BaseNodes
			int node_nb;
			file >> node_nb;
			for (int i = 0; i < node_nb; i++)
			{
				BaseNode node;
				file >> node.layer >> node.index >> node.used >> node.plane_d;

				int neighbor_nb;
				file >> neighbor_nb;
				node.neighbors = std::vector<int>(neighbor_nb, -1);
				for (int j = 0; j < neighbor_nb; j++)file >> node.neighbors[j];

				int boundary_nb_0 = 0;
				file >> boundary_nb_0;
				for (int j = 0; j < boundary_nb_0; j++)
				{
					int boundary_nb_1 = 0;
					file >> boundary_nb_1;
					Vector3d1 vec(boundary_nb_1, Vector3d(0.0, 0.0, 0.0));
					for (int k = 0; k < boundary_nb_1; k++)
						file >> vec[k][0] >> vec[k][1] >> vec[k][2];
					node.boundary.emplace_back(vec);
				}
				base_nodes.emplace_back(node);
			}

			file.clear();
			file.close();



			return true;
		};

		struct BaseEdge
		{
			int index;
			int node_0;
			int node_1;
		};

		if (output_rerun || !LoadBaseNodesSequences(path + "output\\nodes_sequences.txt", base_nodes, layer_nodes,mst))
		{
			layer_nodes = std::vector<std::vector<int>>(layer_boundaries.size(), std::vector<int>());

			std::vector<BaseEdge> base_edges;

			for (int layer_index = 0; layer_index < layer_boundaries.size(); layer_index++)
			{
				auto &layer_boundary = layer_boundaries[layer_index];

				if (layer_boundary.size() == 1)
				{
					BaseNode n;
					n.boundary = layer_boundary;
					n.index = base_nodes.size();
					n.layer = layer_index;
					n.plane_d = plane_d[layer_index];
					layer_nodes[layer_index].emplace_back(n.index);
					base_nodes.emplace_back(n);
				}
				else
				{
					std::vector<bool> bbb(layer_boundary.size(), true);

					while (true)
					{
						bool goon = false;
						for (auto b : bbb)
						{
							if (b) goon = true;
						}
						if (!goon) break;

						std::vector<std::vector<int>> children(layer_boundary.size(), std::vector<int>());
						std::vector<std::vector<int>> parent(layer_boundary.size(), std::vector<int>());
						for (int j = 0; j < children.size(); j++)
						{
							for (int k = j + 1; k < children.size(); k++)
							{
								if (bbb[j] && bbb[k])
								{
									auto outside_py = Math::Vector3d2d(layer_boundary[j]);
									auto inside_py = Math::Vector3d2d(layer_boundary[k]);

									if (CGAL_2D_Detect_Polygon_Inside(outside_py, inside_py))
									{
										children[j].emplace_back(k);
										parent[k].emplace_back(j);
									}

									if (CGAL_2D_Detect_Polygon_Inside(inside_py, outside_py))
									{
										children[k].emplace_back(j);
										parent[j].emplace_back(k);
									}
								}
							}
						}

						while (true)
						{
							int index = -1;
							for (int j = 0; j < parent.size(); j++)
							{
								if (parent[j].empty() && bbb[j])
								{
									index = j;
									break;
								}
							}
							if (index < 0)break;

							bbb[index] = false;

							Vector3d2 dsd;
							dsd.emplace_back(layer_boundary[index]);
							for (int j = 0; j < children[index].size(); j++)
							{
								if (parent[children[index][j]].size() == 1)
								{
									dsd.emplace_back(layer_boundary[children[index][j]]);
									bbb[children[index][j]] = false;
								}
							}

							BaseNode n;
							n.boundary = dsd;
							n.index = base_nodes.size();
							n.layer = layer_index;
							n.plane_d = plane_d[layer_index];
							layer_nodes[layer_index].emplace_back(n.index);
							base_nodes.emplace_back(n);
						}
					}
				}
			}

			//directions
			for (auto &base_node : base_nodes)
			{
				for (int i = 0; i < base_node.boundary.size(); i++)
				{
					bool b = CGAL_2D_Polygon_Is_Clockwise_Oriented(Math::Vector3d2d(base_node.boundary[i]));
					if ((i == 0 && b) || (i != 0 && !b))
						std::reverse(base_node.boundary[i].begin(), base_node.boundary[i].end());
				}
			}

			//debug
			if (output_debug_offset)
			{
				for (int i = 0; i < base_nodes.size(); i++)
				{
					auto name = "base_node_" + IntString(i) + "_" + std::to_string(base_nodes[i].layer);
					auto a = PosApplyM(base_nodes[i].boundary, glm::inverse(rm));
					Circuit::OutputOffsets(path + "boundary\\" + name + ".obj", a, name);
					name = "base_node_tri_" + IntString(i) + "_" + std::to_string(base_nodes[i].layer);
					OutputNodeTriangles(path + "boundary\\" + name + ".obj", base_nodes[i].boundary, rm);
				}
			}

			
			for (int layer = 0; layer < layer_boundaries.size() - 1; layer++)
			{
				auto current_nodes = layer_nodes[layer];
				auto next_nodes = layer_nodes[layer + 1];

				for (auto next_node : next_nodes)
				{
					bool goon = true;
					double min_dis = 10000000000.0;
					int min_cur_node = -1;
					for (auto cur_node : current_nodes)
					{
						auto outside_boundary = Math::Vector3d2d(base_nodes[cur_node].boundary);
						auto inside_boundary = Math::Vector3d2d(base_nodes[next_node].boundary);

						double dis = CGAL_2D_Distance_Polygon_Polygon(outside_boundary, inside_boundary);

						if (dis < min_dis)
						{
							min_dis = dis;
							min_cur_node = cur_node;
						}

						std::vector<std::pair<int,int>> inside_outer;
						int all_inside_nb = 0;
						for (int i = 0; i < inside_boundary.size(); i++)
						{
							for (int j = 0; j < inside_boundary[i].size(); j++)
							if (!CGAL_2D_Detect_Polygon_Inside(outside_boundary, inside_boundary[i][j]))
							{
								inside_outer.emplace_back(i, j);
							}
							all_inside_nb += inside_boundary[i].size();
						}
			
						double rate = (double)inside_outer.size() / (double)all_inside_nb;

						if (rate <1.0)
						{
							//edge
							BaseEdge be;
							be.index = base_edges.size();
							be.node_0 = cur_node;
							be.node_1 = next_node;
							base_edges.emplace_back(be);
							mst.emplace_back(cur_node);
							mst.emplace_back(next_node);

							base_nodes[cur_node].PushNeighbors(next_node);
							base_nodes[next_node].PushNeighbors(cur_node);

							goon = true;
							continue;
						}
					}

				
					if (goon&&min_dis< layer_connecting_threshold &&min_cur_node >= 0)
					{
						BaseEdge be;
						be.index = base_edges.size();
						be.node_0 = min_cur_node;
						be.node_1 = next_node;
						base_edges.emplace_back(be);
						mst.emplace_back(min_cur_node);
						mst.emplace_back(next_node);

						base_nodes[min_cur_node].PushNeighbors(next_node);
						base_nodes[next_node].PushNeighbors(min_cur_node);
					}
				}
			};

			//debug
			Output_tree(base_nodes.size(), mst, "Z:\\Documents\\Windows\\SmartSFC\\workspace\\CFS\\layer.gml");

			OutputBaseNodesSequences(path + "output\\nodes_sequences.txt", base_nodes, layer_nodes,mst);
		}

	};


	void RoughMachining::SegmentMeshReconstruct(const std::string &path, const std::vector<double> &plane_d, const Vector3d &plane_n, const double &layer_thichness,
		const std::vector<BaseNode> &base_nodes, const std::vector<std::vector<int>> &sequences, const glm::dmat4 &rm)
	{
		auto delta_plane_d = plane_n[2] / abs(plane_n[2])*layer_thichness;

		for (int seg_index = 0; seg_index < sequences.size(); seg_index++)
		{
			//volume_index

			std::ofstream o_file;

			if (volume_index<0)
				o_file.open(path + "boundary\\carvable_volume_" + std::to_string(seg_index + 1) + ".obj");
			else
				o_file.open(path + "boundary\\carvable_volume_" + std::to_string(volume_index) 
				+ "_" + std::to_string(seg_index + 1) + ".obj");

			int o_index = 1;

			for (int i = 0; i < sequences[seg_index].size(); i++)
			{
				auto node_index = sequences[seg_index][i];
				auto cur_points = base_nodes[node_index].boundary;
				auto cur_points_2d = Math::Vector3d2d(cur_points);
				auto pre_points = Math::Vector2d3d(cur_points_2d, base_nodes[node_index].plane_d - delta_plane_d);

				if (false)
				{
					auto pre_points_2d = cur_points_2d;
					auto outside_points = base_nodes[sequences[seg_index][i - 1]].boundary;
					auto outside_points_2d = Math::Vector3d2d(outside_points);

					for (int j = 0; j < pre_points_2d.size(); j++)
					for (int k = 0; k < pre_points_2d[j].size(); k++)
						pre_points_2d[j][k] = CGAL_2D_Nearest_Point_Polygon(pre_points_2d[j][k], outside_points_2d);
					pre_points = Math::Vector2d3d(pre_points_2d, base_nodes[node_index].plane_d - delta_plane_d);
				}

				auto triangles = CGAL_2D_Polygon_Triangulation(cur_points_2d);

				//upper face
				int nb = 0;
				std::vector<std::vector<int>> cur_vecs;
				for (int j = 0; j < cur_points.size(); j++)
				{
					cur_vecs.emplace_back(std::vector<int>());
					for (int k = 0; k < cur_points[j].size(); k++)
					{
						cur_vecs.back().emplace_back(nb + o_index);
						auto v = PosApplyM(cur_points[j][k],glm::inverse(rm));
						o_file << "v " << v[0] << " " << v[1] << " " << v[2] << std::endl;
						nb++;
					}
				}

				for (int j = 0; j < triangles.size(); j++)
					o_file << "f " << triangles[j][2] + o_index << " " << triangles[j][1] + o_index << " " << triangles[j][0] + o_index << std::endl;
				o_index += nb;

				//lower face
				nb = 0;
				std::vector<std::vector<int>> pre_vecs;
				for (int j = 0; j < pre_points.size(); j++)
				{
					pre_vecs.emplace_back(std::vector<int>());
					for (int k = 0; k < pre_points[j].size(); k++)
					{
						pre_vecs.back().emplace_back(nb  + o_index);
						auto v = PosApplyM(pre_points[j][k], glm::inverse(rm));
						o_file << "v " << v[0] << " " << v[1] << " " << v[2] << std::endl;
						nb++;
					}
				}

				for (int j = 0; j < triangles.size(); j++)
					o_file << "f " << triangles[j][0] + o_index << " " << triangles[j][1] + o_index << " " << triangles[j][2] + o_index << std::endl;
				o_index += nb;


				//side surfaces (points)
				//for (int j = 0; j < cur_points.size(); j++)
				//{
				//	for (int k = 0; k < cur_points[j].size(); k++)
				//	{
				//		auto p0 = cur_points[j][k];
				//		auto p1 = cur_points[j][(k + 1) % cur_points[j].size()];
				//		auto p2 = pre_points[j][k];
				//		auto p3 = pre_points[j][(k + 1) % cur_points[j].size()];

				//		o_file << "v " << p0[0] << " " << p0[1] << " " << p0[2] << std::endl;
				//		o_file << "v " << p1[0] << " " << p1[1] << " " << p1[2] << std::endl;
				//		o_file << "v " << p2[0] << " " << p2[1] << " " << p2[2] << std::endl;
				//		o_file << "v " << p3[0] << " " << p3[1] << " " << p3[2] << std::endl;
				//	}
				//}

				//side surfaces (faces)
				for (int j = 0; j < cur_points.size(); j++)
				{
					for (int k = 0; k < cur_points[j].size(); k++)
					{
						auto a = cur_vecs[j][k];
						auto b = cur_vecs[j][(k + 1) % cur_points[j].size()];
						auto c = pre_vecs[j][k];
						auto d = pre_vecs[j][(k + 1) % cur_points[j].size()];

						o_file << "f " << d << " " << b << " " << a << std::endl;
						o_file << "f " << c << " " << d << " " << a << std::endl;
					}
				}
			}

			o_file.clear();
			o_file.close();
		}
	};


	void RoughMachining::SegmentMeshReconstruct2(const std::string &path, const std::vector<double> &plane_d, const Vector3d &plane_n, const double &layer_thichness,
		const std::vector<BaseNode> &base_nodes, const std::vector<std::vector<int>> &sequences, const glm::dmat4 &rm)
	{
		auto SearchNeighbor = [](const std::vector<BaseNode> &base_nodes, const std::vector<int> &sequence)
		{
			std::vector<int> neighbors;
			for (int i = 0; i < sequence.size(); i++)
			{
				const BaseNode &node = base_nodes[sequence[i]];
				for (int j = 0; j < node.neighbors.size(); j++)
				if (std::find(sequence.begin(), sequence.end(), node.neighbors[j]) == sequence.end())
				{
					neighbors.emplace_back(node.neighbors[j]);
				}
			}
			return neighbors;
		};

		auto GenerateCutComponent = [&](const std::vector<BaseNode> &base_nodes, const BaseNode & h_node, const BaseNode & l_node, 
			Vector3d2 &seg_true, Vector3d2 &seg_false)
		{
			//get bro nodes
			std::vector<int> bro_nodes;
			for (int j = 0; j < h_node.neighbors.size(); j++)
			{
				const auto &node = base_nodes[h_node.neighbors[j]];
				if (node.layer == l_node.layer && node.index != l_node.index) bro_nodes.emplace_back(node.index);
			}

			if (bro_nodes.empty()) return false;

			//get bro boundary
			Vector2d2 bro_boundary;
			for (auto & bro_node : bro_nodes)
			{
				const auto &node = base_nodes[bro_node];
				for (const auto&boundary : node.boundary)bro_boundary.emplace_back(Math::Vector3d2d(boundary));
			}

			//cur boundary
			Vector2d2 cur_boundary = Math::Vector3d2d(l_node.boundary);

			//pare boundary
			Vector2d2 pare_boundary = Math::Vector3d2d(h_node.boundary);
			std::vector<std::vector<bool>> indicator(pare_boundary.size(),std::vector<bool>());
			Vector3d1 points_0, points_1;
			std::vector<std::string> name_0, name_1;
			for (int i = 0; i < pare_boundary.size(); i++)
			{
				for (int j = 0; j < pare_boundary[i].size(); j++)
				{
					//CGAL_2D_Distance_Point_Polygon(Vector2d p, Vector2d2 pys)
					double curd = CGAL_2D_Distance_Point_Polygon(pare_boundary[i][j], cur_boundary);
					double brod = CGAL_2D_Distance_Point_Polygon(pare_boundary[i][j], bro_boundary);
					indicator[i].push_back(curd <= brod);

					if (indicator[i].back())
					{
						points_0.emplace_back(h_node.boundary[i][j]);
						name_0.emplace_back("p_"+std::to_string(i)+"_"+std::to_string(j));
					}
					else
					{
						points_1.emplace_back(h_node.boundary[i][j]);
						name_1.emplace_back("p_" + std::to_string(i) + "_" + std::to_string(j));
					}
				}
			}

			Circuit::OutputPoints(path + "debug\\debug_0.obj", points_0, name_0);
			Circuit::OutputPoints(path + "debug\\debug_1.obj", points_1, name_1);

			//connecting
			std::vector<std::vector<int>> sequences;
			Vector2d2 seges_true;
			Vector2d2 seges_false;
			for (int i = 0; i < indicator.size(); i++)
			{
				Vector2d1 strip;
				std::vector<bool> used(indicator[i].size(),false);

				while (true)
				{
					int start_index = -1;
					for (int j = 0; j < used.size(); j++)
					{
						if (!used[j])
						{
							start_index = j;
							break;
						}
					}
					if (start_index < 0)break;

					//front
					Vector2d1 segs;
					bool b = indicator[i][start_index];
					for (int j = start_index; j < start_index + used.size(); j++)
					{
						if (indicator[i][j] == b&&!used[j])
						{
							segs.emplace_back(pare_boundary[i][j]);
							used[j] = true;
						}
						else
						{
							if (indicator[i][j] != b) segs.emplace_back((pare_boundary[i][j]+segs.back())/(float)2.0);
							break;
						}
					}

	

					//back
					for (int j = (start_index - 1); (j + used.size()) % used.size() != start_index; j--)
					{
						auto jj = (j + used.size()) % used.size();
						if (indicator[i][jj] == b&&!used[jj])
						{
							segs.insert(segs.begin() + 0, pare_boundary[i][jj]);
							used[jj] = true;
						}
						else
						{
							if (indicator[i][jj] != b) segs.insert(segs.begin(),(pare_boundary[i][jj] + segs.front()) / (float)2.0);
							break;
						}
					}

					if (b) seges_true.emplace_back(segs); else seges_false.emplace_back(segs);
				}
			}

			auto ConnectingSegs = [](Vector2d2 &seges, const Vector2d2 &pare_boundary)
			{
				auto CheckInside = [](const Vector2d2 &pare_boundary, const Vector2d &s, const Vector2d &e)
				{
					//bool CGAL_2D_Detect_Polygon_Inside(Vector2d2 outside_pys, Vector2d pppp)
					for (double d = 0.2; d < 1.0; d += 0.2)
					{
						if (!CGAL_2D_Detect_Polygon_Inside(pare_boundary, (float)d*s + (float)(1.0 - d)*e))
							return false;
					}
					return true;
				};

				Vector2d2 output_segs;
				std::vector<bool> used(seges.size(),false);

				while (true)
				{
				
					int start_index = -1;
					for (int i = 0; i < used.size(); i++)
					{
						if (!used[i])
						{
							start_index = i;
							break;
						}
					}
					if (start_index < 0)break;

					Vector2d1 seg = seges[start_index];
					used[start_index] = true;

					while (true)
					{
						int jj = -1;
						int kk = -1;
						double min_d = 100000000000000.0;
						for (int j = 0; j < used.size(); j++)
						{
							if (!used[j])
							{
								double d_front = CGAL_2D_Distance_Point_Point(seg.back(), seges[j].front());
								double d_back = CGAL_2D_Distance_Point_Point(seg.back(), seges[j].back());

								if (CheckInside(pare_boundary, seg.back(), seges[j].front()))
								{
									if (d_front < min_d)
									{
										jj = j;
										kk = 0;
										min_d = d_front;
									}
								}

								if (CheckInside(pare_boundary, seg.back(), seges[j].back()))
								{
									if (d_back < min_d)
									{
										jj = j;
										kk = 1;
										min_d = d_back;
									}
								}
							}
						}

						if (jj < 0 || kk<0)break;
						if (kk == 1) std::reverse(seges[jj].begin(), seges[jj].end());
						for (auto p : seges[jj]) seg.emplace_back(p);
						used[jj] = true;
					}
					output_segs.emplace_back(seg);
				}


			


				return output_segs;
			};



			seg_true = Math::Vector2d3d(ConnectingSegs(seges_true, pare_boundary), h_node.plane_d);
			seg_false = Math::Vector2d3d(ConnectingSegs(seges_false, pare_boundary), h_node.plane_d);

			Circuit::OutputOffsets(path + "debug\\seg_true.obj", seg_true, "seg_true");
			Circuit::OutputOffsets(path + "debug\\seg_false.obj", seg_false, "seg_false");

			int d = 0;

			return true;
		};

		auto SearchNeighbor_1 = [&](const std::vector<BaseNode> &base_nodes, 
			const std::vector<int> &sequence, Vector3d2 &h_multi_boundary, Vector3d2 &l_multi_boundary)
		{
			for (int i = 0; i < sequence.size(); i++)
			{
				const BaseNode &cur_node = base_nodes[sequence[i]];
				if (i == 0)
				{
					for (int j = 0; j < cur_node.neighbors.size(); j++)
					{
						const BaseNode &h_node = base_nodes[cur_node.neighbors[j]];
						if (std::find(sequence.begin(), sequence.end(), h_node.index) == sequence.end())
						{
							if (h_node.layer < cur_node.layer)
							{
								auto mm_center = CGAL_3D_Mesh_Center(h_node.boundary);
								Vector3d v0 = mm_center + Vector3d(-100.0, -100.0, 0.0);
								Vector3d v1 = mm_center + Vector3d( 100.0, -100.0, 0.0);
								Vector3d v2 = mm_center + Vector3d( 100.0, 100.0, 0.0);
								Vector3d v3 = mm_center + Vector3d(-100.0, 100.0, 0.0);
								h_multi_boundary.emplace_back(Vector3d1{ v0, v1, v2, v3 });

								Vector3d2 seg_true, seg_false;
								if (GenerateCutComponent(base_nodes, h_node, cur_node, seg_true, seg_false))
								{
									for (auto seg : seg_false)l_multi_boundary.emplace_back(seg);
								}
								//l_multi_boundary.emplace_back(seg_false);
							}
						}
					}
				}
				
				int lower_node;
				for (int j = 0; j < cur_node.neighbors.size(); j++)
				{
					const BaseNode &next_node = base_nodes[cur_node.neighbors[j]];
					if (std::find(sequence.begin(), sequence.end(), next_node.index) != sequence.end() && next_node.layer>cur_node.layer)
					{
						lower_node = next_node.index;
					}
				}

				if (i < sequence.size() - 1)
				{
					auto & l_node = base_nodes[sequence[i+1]];
					Vector3d2 seg_true, seg_false;
					if (GenerateCutComponent(base_nodes, cur_node, l_node, seg_true, seg_false))
					{
						for (auto seg : seg_false)l_multi_boundary.emplace_back(seg);
					}
				}
	
			}
		};

		for (int i = 0; i < sequences.size(); i++)
		{
			std::cerr << i << " : ";
			for (const auto & id : sequences[i])  std::cerr << id << " ";
			std::cerr << "\n";

			auto neighbors = SearchNeighbor(base_nodes, sequences[i]);
			std::cerr << i << " : ";
			for (const auto & neighbor : neighbors)  std::cerr << neighbor << " ";
			std::cerr << "\n";

			Vector3d2 h_multi_boundary, l_multi_boundary;
			SearchNeighbor_1(base_nodes, sequences[i], h_multi_boundary, l_multi_boundary);

			if (!h_multi_boundary.empty())
				OutputPolygon(path + "output\\carvable_volume_h_"+std::to_string(i)+".obj", h_multi_boundary,rm);
			OutputPolygon(path + "output\\carvable_volume_l_" + std::to_string(i) + ".obj", l_multi_boundary,rm);
			auto inside_node = base_nodes[sequences[i][(int)((double)sequences[i].size() / 2.0)]];
			//Circuit::OutputOffsets(path+"\\boundary\\multi_boundary.obj",multi_boundary,"multi_boundary");
			//CGAL_Cut_Surface_by_Multi_Boundaries(multi_boundary, inside_node.boundary.front().front(), path + "output\\initial_cut.off", path + "output\\volume_" + std::to_string(i) + ".obj");

		}
	};


	void RoughMachining::OutputCLS(const std::string folder_name, const std::string output_path, const Vector3d1 &rough_machining_path, const Vector3d &plane_n)
	{
		std::ofstream file(output_path);

		file << "TOOL PATH/CAVITY_MILL,TOOL,MILL" << std::endl;
		file << "TLDATA/MILL,3.1750,0.0000,50.0000, 0.0000,0.0000" << std::endl;
		file << "SYS/0.0000,0.0000,0.0000,1.0000000,0.0000000,0.0000000,0.0000000,1.0000000,0.0000000" << std::endl;
		file << "$$ centerline data" << std::endl;
		file << "PAINT/PATH" << std::endl;
		file << "PAINT/SPEED,10" << std::endl;
		file << "PAINT/COLOR,186" << std::endl;



		auto cut_normal = -plane_n;
		Math::SetVectorLength(cut_normal,1.0);

		auto s = rough_machining_path.front();

		file << "RAPID" << std::endl;
		file << "GOTO/" << std::to_string(s[0]) << "," << std::to_string(s[1]) << "," << std::to_string(s[2])<<","<<
			std::to_string(cut_normal[0]) << "," << std::to_string(cut_normal[1]) << "," << std::to_string(cut_normal[2])<< std::endl;
		file << "PAINT/COLOR,211" << std::endl;

		file << "PAINT/COLOR,42" << std::endl;
		file << "FEDRAT/MMPM,5000.0000" << std::endl;

		for (int i = 1; i < rough_machining_path.size(); i++)
		{
			auto v = rough_machining_path[i];
			file << "GOTO/" << std::to_string(v[0]) << "," << std::to_string(v[1]) << "," << std::to_string(v[2]) << std::endl;
		}

		//RAPID
		//	GOTO / 25.000, 60.0000, 25.0, 0.0000000, 1.0000000, 0.0000000
		//	PAINT / COLOR, 211
		//	RAPID
		//	GOTO / 25.000, 50.0000, 25.0
		//	PAINT / COLOR, 42
		//	FEDRAT / MMPM, 250.0000
		//	GOTO / 50.4000, 50.0, 25.0
		//	RAPID
		//	GOTO / 25.000, 60.0000, 25.0

		file << "PAINT/SPEED,10" << std::endl;
		file << "PAINT/TOOL,NOMORE" << std::endl;
		file << "END-OF-PATH" << std::endl;

		file.clear();
		file.close();
	}

	void RoughMachining::OutputFinalNGC(const std::string output_path, const Vector3d1 &rough_machining_path, const Vector3d &plane_n)
	{
		//glm::dmat4 rm = RotationMatrix(Vector3d(0.0, 0.0, -1.0), plane_n, Vector3d(0.0,1.0,0.0));
		
		auto GetRMAB = [](const Vector3d &n, glm::dmat4 &m,double &A,double &B)
		{
			Vector3d x,y,z;
			z = n;
			double angle = Math::GetAngleBetween(z, Vector3d(0.0, 1.0, 0.0));
			angle = RadiantoAngle(angle);
			if (Math::IsAlmostZero(angle))
			{
				z = Vector3d(0.0, 1.0, 0.0);
				x = Vector3d(1.0, 0.0, 0.0);
				y = Vector3d(0.0, 0.0, -1.0);
			}
			else
			{
				Vector3d axis = GetCrossproduct(z, Vector3d(0.0, 1.0, 0.0));
				y = RotationAxis(z, MM_PI / 2.0, axis);
				x = GetCrossproduct(y, z);
			}
			Math::ClearVector3d(x);
			Math::ClearVector3d(y);
			Math::ClearVector3d(z);
			x = x / (float)Math::GetLength(x);
			y = y / (float)Math::GetLength(y);
			z = z / (float)Math::GetLength(z);

			//m
			m = RotationMatrixXYZ(x, y, z);

			//AB
			if (Math::IsAlmostZero(angle))
			{
				A = 90.0;
				B = 0.0;
			}
			else
			{
				A = 90.0 - angle;
				B = Math::GetAngleBetween(Vector3d(z[0],0.0,z[2]), Vector3d(0.0, 0.0, 1.0));
				B = RadiantoAngle(B);
				if (z[0] > 0.0)B = 360.0 - B;
			}
		};

		glm::dmat4 rm;
		double A, B;
		GetRMAB(-plane_n, rm, A, B);

		Vector3d1 vecs;
		for (auto p : rough_machining_path)
		{
			//mm => inch (expert)
			Vector3d inch_p(p[0] / 25.4, p[1] / 25.4, p[2] / 25.4);
			//expert => fusion 3d global 
			Vector3d fusion_g(inch_p[0] - 0.815, inch_p[1] + 1.563, inch_p[2] - 1.0);

			// -0.999: -25.3746
			// 1.563:39.7002
			// -0.82: -20.828

			//=>setup coordinate
			//-center
			fusion_g[1] = fusion_g[1] - 1.859;
			//rotation
			//Vector3d rotate_v = RotationAxis(fusion_g, MM_PI / 2.0, Vector3d(0.0, 1.0, 0.0));
			Vector3d rotate_v = fusion_g;
			//=>tool path rotation
			//rotate_v = RotationAxis(rotate_v, MM_PI / 2.0*(plane_n.z), Vector3d(0.0, 1.0, 0.0));
			rotate_v = PosApplyM(rotate_v, rm);
			vecs.emplace_back(rotate_v);
		}

		//output.open( filename.c_str(), ios::out | ios::app )

		bool head = false;
		if (_access("Z:\\Documents\\Windows\\SmartSFC\\workspace\\path.ngc", 0) == -1)
		{
			head = true;
		}

		std::ofstream full_fill("Z:\\Documents\\Windows\\SmartSFC\\workspace\\path.ngc",  ios::out | ios::app);


		std::ofstream file(output_path);
		file << "%" << std::endl;
		file << "(AXIS, stop)" << std::endl;
		file << "(EXPERT_BUNNY)" << std::endl;
		file << "(V2 FIRST OUTPUT)" << std::endl;
		file << "G20" << std::endl;
		file << "G90 G94 G40 G17 G91.1" << std::endl;
		file << "G53 G0 Z0." << std::endl;
		file << "(POCKET4 3)" << std::endl;
		file << "G49" << std::endl;
		file << "M5" << std::endl;
		file << "G53 G0 X2.5 Y2.5" << std::endl;
		file << "M0" << std::endl;
		file << "T1 M6" << std::endl;
		file << "(0)" << std::endl;
		file << "S12000 M3" << std::endl;
		file << "G54 G0" << std::endl;

		if (head)
		{
			full_fill << "%" << std::endl;
			full_fill << "(AXIS, stop)" << std::endl;
			full_fill << "(EXPERT_BUNNY)" << std::endl;
			full_fill << "(V2 FIRST OUTPUT)" << std::endl;
			full_fill << "G20" << std::endl;
			full_fill << "G90 G94 G40 G17 G91.1" << std::endl;
			full_fill << "G53 G0 Z0." << std::endl;
			full_fill << "(POCKET4 3)" << std::endl;
			full_fill << "G49" << std::endl;
			full_fill << "M5" << std::endl;
			full_fill << "G53 G0 X2.5 Y2.5" << std::endl;
			full_fill << "M0" << std::endl;
			full_fill << "T1 M6" << std::endl;
			full_fill << "(0)" << std::endl;
			full_fill << "S12000 M3" << std::endl;
			full_fill << "G54 G0" << std::endl;
		}

	//	if (plane_n.z == -1)
		//	file << "A0.B-90." << std::endl;

		//if (plane_n.z == 1)
		//	file << "A0.B90." << std::endl;

		//////////////////////////////////////////////////////////
		file << "A"+std::to_string(A)+"B"+std::to_string(B) << std::endl;
		full_fill << "A" + std::to_string(A) + "B" + std::to_string(B) << std::endl;

		auto first_point = vecs.front();
		file << "G1 X" + std::to_string(first_point[0]) + " Y" + std::to_string(first_point[1]) + " F5000." << std::endl;
		file << "G43 Z" + std::to_string(vecs.front()[2]) + " H1" << std::endl;

		full_fill << "G1 X" + std::to_string(first_point[0]) + " Y" + std::to_string(first_point[1]) + " F5000." << std::endl;
		full_fill << "G43 Z" + std::to_string(vecs.front()[2]) + " H1" << std::endl;

		//expert CO => 
		for (auto s : vecs)
		{
			file << "G1 X" + std::to_string(s[0]) + " Y" + std::to_string(s[1]) + " Z" + std::to_string(s[2]) << std::endl;
			full_fill << "G1 X" + std::to_string(s[0]) + " Y" + std::to_string(s[1]) + " Z" + std::to_string(s[2]) << std::endl;
		}

		file << "G53 G0 Z0" << std::endl;
		full_fill << "G53 G0 Z0" << std::endl;

		////////////////////////////////////////////////////////////

		file << "G49" << std::endl;
		file << "G53 G0 X2.5 Y2.5" << std::endl;
		file << "A0.B0." << std::endl;
		file << "M30" << std::endl;
		file << "(AXIS, stop)" << std::endl;
		file << "%" << std::endl;


		file.clear();
		file.close();

		full_fill.clear();
		full_fill.close();

	};


	Vector3d1 RoughMachining::GenerateToolpathForOneVolume(const int &segment_index, const std::string &path,
		const Vector3d3 &rough_boundaries, const std::vector<double> &plane_d, const bool &re_running, const glm::dmat4 &rm)
	{
		auto StartEndSelection = [](const Vector3d3 &rough_boundaries, const std::vector<double> &plane_d)
		{
			//3d 2d
			Vector2d3 rough_boundaries_2d;
			for (int i = 0; i < rough_boundaries.size(); i++)
				rough_boundaries_2d.emplace_back(Math::Vector3d2d(rough_boundaries[i]));

			//final boundary
			auto final_boundary = rough_boundaries_2d.back();

			Vector2d inner_vec;
			bool success_b = CGAL_2D_Polygon_Inside_Point(final_boundary, inner_vec);

			if (!success_b)
				std::cerr << "error" << std::endl;

			Vector3d1 start_ends;
			for (int i = 0; i < plane_d.size(); i++)
				start_ends.emplace_back(Math::Vector2d3d(inner_vec, plane_d[i]));

			return start_ends;
		};

		//return Vector3d1();

		//start exist points
		auto start_ends = StartEndSelection(rough_boundaries, plane_d);

		///generate tool path
		Vector3d2 layer_pathes;
		for (int i = 0; i < rough_boundaries.size(); i++)
		{
			//if (plane_d[i] != 31)continue;

			//i = 13;
			std::cerr << "Generate CFS for the layer: " << Math::IntString(i) << " plane_d:" << std::to_string(plane_d[i]) << std::endl;

			cfs.ClearAll();
			cfs.cfs_index = i;
			//OffsetsBasedCFSLinking_Framework(path, rough_boundaries[i], re_running);
			cfs.OffsetsBasedCFSLinking(path, rough_boundaries[i], re_running);

			int location_index; double par;
			Circuit::FindNearestPoint(start_ends[i], cfs.single_path, location_index, par);

			if (cfs.single_path.empty())
			{
				std::cerr << "if (cfs.single_path.empty())" << std::endl;
				system("pause");
				cfs.ClearAll();
				cfs.cfs_index = i;
				//OffsetsBasedCFSLinking_Framework(path, rough_boundaries[i], re_running);
				cfs.OffsetsBasedCFSLinking(path, rough_boundaries[i], re_running);
			}

			Vector3d1 new_single_path;
			for (int j = location_index; j < cfs.single_path.size(); j++)
				new_single_path.emplace_back(cfs.single_path[j]);
			for (int j = 0; j <= location_index; j++)
				new_single_path.emplace_back(cfs.single_path[j]);

			auto layer_path = Math::Vector2d3d(Math::Vector3d2d(new_single_path), plane_d[i]);

			Strip::SmoothingLines(layer_path, 3);
			Strip::ConnectingSameDirectionLines(layer_path);


			auto v = PosApplyM(layer_path, glm::inverse(rm));
			auto pre_name = path + "path\\layer_spiral_" + std::to_string(segment_index) + "_" + std::to_string(i);

			cfs.OutputStripNGC(pre_name + ".ngc", v);
			Circuit::OutputStrip(pre_name + ".obj", v);


			Circuit::OutputOffsets(path + "path\\layer_offsets_" + std::to_string(segment_index) + "_" + std::to_string(i)+".obj", cfs.offsets, "offsets");

			layer_pathes.emplace_back(layer_path);
			//break;
		}

		Vector3d1 rough_machining_path;
		for (int i = 0; i < layer_pathes.size(); i++)
		{
			if (i >= 1)
			{
				auto helix_vecs = Helix(rough_machining_path.back(), layer_pathes[i].front());
				for (auto v : helix_vecs)
					rough_machining_path.emplace_back(v);
			}

			for (int j = 0; j < layer_pathes[i].size(); j++)
				rough_machining_path.emplace_back(layer_pathes[i][j]);

			if (i < layer_pathes.size() - 1)
			{
				auto next_v = layer_pathes[(i + 1)].front();
				next_v[2] = plane_d[i];
				rough_machining_path.emplace_back(next_v);
			}
		}

		//Strip::AdaptiveSampling(rough_machining_path, 0.002);
		return rough_machining_path;
	}


}

//
//if (false)
//{
//	/////////////////////////////////////////////////////////
//	std::ofstream file("D:\\program\\win10\\CycleGrouping_v1.0.0\\CycleGrouping_v1.0.0\\exe_data\\Data\\Hip\\CTR\\haisen1.contour");
//	file.precision(7);
//	file << sequences[i].size() << std::endl;;
//	//for (int j = 0; j < 9; j++)
//	for (int j = 0; j < sequences[i].size(); j++)
//	{
//		auto base_node = base_nodes[sequences[i][j]];
//
//		for (int k = 0; k < base_node.boundary.size(); k++)
//		{
//			//if (true)
//			if (k == 0)
//			{
//				auto &a = base_node.boundary[k];
//				Vector2d1 aaa;
//				CGAL_3D_Plane_3D_to_2D_Points(Vector3d(0.0, 0.0, 0.0), Vector3d(0.0, 0.0, 1.0), a, aaa);
//				if (CGAL_2D_Polygon_Is_Clockwise_Oriented(aaa))
//					std::reverse(a.begin(), a.end());
//			}
//			else
//			{
//				auto &a = base_node.boundary[k];
//				Vector2d1 aaa;
//				CGAL_3D_Plane_3D_to_2D_Points(Vector3d(0.0, 0.0, 0.0), Vector3d(0.0, 0.0, 1.0), a, aaa);
//				if (!CGAL_2D_Polygon_Is_Clockwise_Oriented(aaa))
//					std::reverse(a.begin(), a.end());
//			}
//		}
//
//		double a, b, c, d;
//		CGAL_3D_Plane_ABCD(Vector3d(0.0, 0.0, base_node.plane_d), plane_n, a, b, c, d);
//		//d = -d;
//		file << a << " " << b << " " << c << " " << d << std::endl;
//		int nb = 0;
//		for (auto a : base_node.boundary)for (auto b : a)nb++;
//		file << nb << " " << nb << std::endl;
//		for (auto a : base_node.boundary)for (auto b : a) file << b[0] << " " << b[1] << " " << b[2] << std::endl;
//		nb = 0;
//		for (int k = 0; k < base_node.boundary.size(); k++)
//		{
//			for (int l = 0; l < base_node.boundary[k].size(); l++)
//			{
//				file << l + nb << " " << nb + (l + 1) % base_node.boundary[k].size() << " 0 1" << std::endl;
//			}
//			nb += base_node.boundary[k].size();
//		}
//	}
//
//	file.clear();
//	file.close();
//	/////////////////////////////////////////////////////////
//}