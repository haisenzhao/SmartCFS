#include "stdafx.h"
#include "CFSCNC.h"
#include <Circuit.h>
#include "Tree.h"
#include "AM.h"

namespace cnc {

	//load setting
	void AMLoadSetting(const std::string &path, Vector3d &plane_n, Vector3d &origin_n,
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

	void AMOutputNodeTriangles(std::string file_path, Vector3d2 &cut_offsets, const glm::dmat4 &rm)
	{
		std::ofstream file(file_path);

		auto points_2d = Math::Vector3d2d(cut_offsets);
		auto triangles = CGAL_2D_Polygon_Triangulation(points_2d);

		for (int i = 0; i < cut_offsets.size(); i++)
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
	AMVolume::AMVolume(std::string path, bool re_running, bool output_debug,
		double set_toolpath_size, double set_layer_thickness, bool segment_or_tool_path)
	{
		if (re_running)
		{
			Math::ClearFolder(path + "boundary\\");
			Math::ClearFolder(path + "debug\\");
			Math::ClearFolder(path + "output\\");
			Math::ClearFolder(path + "path\\");
		}


		//threshold
		area_threshold = 3.0;
		length_threshold = 10.0;
		removed_length_threshold = 5.0;
		layer_connecting_threshold = 0.8;
		output_debug_offset = output_debug;
		output_rerun = re_running;


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

		AMLoadSetting(path, plane_n, origin_n, plane_d, cfs.toolpath_size, layer_thichness, rm, obj_file, off_file, set_toolpath_size, set_layer_thickness);
		auto strings = Split(path, "\\");
		std::string folder_name = strings[strings.size() - 2];

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

		for (int layer = 0; layer < layer_nodes.size(); layer++)
		{
			if (layer<=841)continue;

			for (int layer_node = 0; layer_node < layer_nodes[layer].size(); layer_node++)
			{
				auto &node = base_nodes[layer_nodes[layer][layer_node]];

				std::cerr << "Generate CFS for the layer: " << Math::IntString(layer) << "/" << layer_nodes.size() << " plane_d:" << std::to_string(node.plane_d) << std::endl;

				cfs.ClearAll();
				cfs.cfs_index = node.index;
				cfs.OffsetsBasedCFSLinking(path, node.boundary, re_running);

				if (cfs.single_path.empty())
				{
					std::cerr << "if (cfs.single_path.empty())" << std::endl;
					system("pause");
					cfs.ClearAll();
					cfs.cfs_index = node.index;
					//OffsetsBasedCFSLinking_Framework(path, rough_boundaries[i], re_running);
					cfs.OffsetsBasedCFSLinking(path, node.boundary, re_running);
				}

				Strip::SmoothingLines(cfs.single_path, 2);
				Strip::ConnectingSameDirectionLines(cfs.single_path);

				for (auto& s : cfs.single_path)s[2] = node.plane_d;

				auto v = PosApplyM(cfs.single_path, glm::inverse(rm));
				auto pre_name = path + "path\\layer_spiral_" + std::to_string(layer) + "_" + std::to_string(layer_node);

				cfs.OutputStrip(pre_name + ".ngc", v);
				Circuit::OutputStrip(pre_name + ".obj", v);

				//Circuit::OutputOffsets(path + "path\\layer_offsets_" + std::to_string(layer) + "_" + std::to_string(layer_node) + ".obj", cfs.offsets, "offsets");

			
			}
		}

		//227
		//if (CGAL_2D_Detect_Polygon_Inside(outside_py, inside_py))
		//{
		//	children[j].emplace_back(k);
		//	parent[k].emplace_back(j);
		//}
		                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
		//for (auto &base_node : base_nodes)
		//{
		//	for (int i = 0; i < base_node.boundary.size(); i++)
		//	{
		//		bool b = CGAL_2D_Polygon_Is_Clockwise_Oriented(Math::Vector3d2d(base_node.boundary[i]));
		//		if ((i == 0 && b) || (i != 0 && !b))
		//			std::reverse(base_node.boundary[i].begin(), base_node.boundary[i].end());
		//	}
		//}

	}


		void AMVolume::BuildRelationShip(const std::string &path, Vector3d3 &layer_boundaries, const std::vector<double> &plane_d, const glm::dmat4 &rm,
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

			if (output_rerun || !LoadBaseNodesSequences(path + "output\\nodes_sequences.txt", base_nodes, layer_nodes, mst))
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
						AMOutputNodeTriangles(path + "boundary\\" + name + ".obj", base_nodes[i].boundary, rm);
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

							std::vector<std::pair<int, int>> inside_outer;
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
				Output_tree(base_nodes.size(), mst, path+"output\\layer.gml");

				OutputBaseNodesSequences(path + "output\\nodes_sequences.txt", base_nodes, layer_nodes, mst);
			}

		};



	void AMVolume::MeshSlicer(const std::string &path, const std::string off_path,
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
								if (length > removed_length_threshold) layer_removed_offsets.emplace_back(Math::Vector2d3d(clear_points_2d, plane_d[i]));
							}

						}
					}
					else
					{
						if (length > removed_length_threshold) layer_removed_offsets.emplace_back(Math::Vector2d3d(clear_points_2d, plane_d[i]));
					}
				}
				rough_boundaries.emplace_back(layer_offsets);
				removed_offsets.emplace_back(layer_removed_offsets);
			}

			//Circuit::OutputOffsets(path + "boundary\\boudary_removed.obj", removed_offsets,"boudary_removed");
			OutputVectors(path + "output\\mesh_slicer.txt", rough_boundaries);
			OutputVectors(path + "output\\mesh_slicer_removed_offsets.txt", removed_offsets);
		}

	};



}