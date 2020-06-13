#ifndef TREE_ONCE
#define TREE_ONCE
#pragma once

#include <stdafx.h>
#include <MathHelper.h>

#include <Strip.h>
#include <Circuit.h>

#include "kdtree.h"

using namespace hpcg;
using namespace std;

namespace cnc {

	class Tree
	{

	public:

		static int GetOffsetIndex(std::vector<int> &index_int, int i, int j)
		{
			int index = -1;

			for (int l = 0; l < index_int.size() && index<0; l = l + 3)
			{
				if (index_int[l] == i && index_int[l + 1] == j)
				{
					index = index_int[l + 2];
				}
			}
			return index;
		}

		static int GetOffsetGraph(std::vector<int> &offset_graph, int upper_index, int lower_index)
		{
			for (int i = 0; i < offset_graph.size(); i = i + 2)
			{
				if (offset_graph[i] == upper_index&&offset_graph[i + 1] == lower_index)
				{
					return i;
				}
			}

			return -1;
		}

		static void UnifiedContourDirection(double toolpath_size, 
			Vector3d3 &offsetses, 
			Vector3d2 &offsets)
		{
			//build offset graph
			std::vector<int> index_int;

			int index = -1;
			for (int i = 0; i < offsetses.size(); i++)
			{
				for (int j = 0; j < offsetses[i].size(); j++)
				{
					index++;
					index_int.push_back(i);
					index_int.push_back(j);
					index_int.push_back(index);
				}
			}

			//first laryer
			if (offsetses[0].size()>1)
			{
				for (int i = 0; i < offsetses[0].size(); i++)
				{
					int upper_index = GetOffsetIndex(index_int, 0, i);
					for (int j = i + 1; j < offsetses[0].size(); j++)
					{
						int lower_index = GetOffsetIndex(index_int, 0, j);
						Vector3d2 sharing_parts_0;
						Vector3d2 sharing_parts_1;
						bool b0 = GetSharingParts(toolpath_size, offsets[upper_index], offsets[lower_index], sharing_parts_0);
						bool b1 = GetSharingParts(toolpath_size, offsets[lower_index], offsets[upper_index], sharing_parts_1);
						if (b0&&b1)
						{
							if (Circuit::CheckSameDirection(offsets[upper_index], offsets[lower_index], toolpath_size))
							{
								std::reverse(offsets[lower_index].begin(), offsets[lower_index].end());
							}
						}
					}
				}
			}

			/*
			for (int i = 0; i < offsetses.size(); i++)
			{
				for (int j = 0; j < offsetses[i].size() - 1; j++)
				{
					int upper_index = GetOffsetIndex(index_int, i, j);

					for (int k = j + 1; k < offsetses[i].size(); k++)
					{
						int lower_index = GetOffsetIndex(index_int, i, k);

						Vector3d2 sharing_parts_0;
						Vector3d2 sharing_parts_1;

						bool b0 = GetSharingParts(toolpath_size, offsets[upper_index], offsets[lower_index], sharing_parts_0);
						bool b1 = GetSharingParts(toolpath_size, offsets[lower_index], offsets[upper_index], sharing_parts_1);

						if (b0&&b1)
						{
							offset_graph.push_back(upper_index);
							offset_graph.push_back(lower_index);
							offset_graph_sharing_parts.push_back(sharing_parts_0);
							offset_graph_sharing_parts.push_back(sharing_parts_1);
						}
					}
				}
			}
			*/


			//build the offset graph up to down
			for (int i = 0; i < offsetses.size() - 1; i++)
			{
				for (int j = 0; j < offsetses[i].size(); j++)
				{
					int upper_index = GetOffsetIndex(index_int, i, j);

					std::vector<int> lower_indexes;

					for (int k = 0; k < offsetses[i + 1].size(); k++)
					{
						int lower_index = GetOffsetIndex(index_int, i + 1, k);
						if (Circuit::DistanceDouble(offsets[upper_index], offsets[lower_index], toolpath_size))
							lower_indexes.push_back(lower_index);
					}

					if (lower_indexes.size() > 0)
					{
						for (int k = 0; k < lower_indexes.size(); k++)
						{
							if (!Circuit::CheckSameDirection(offsets[upper_index], offsets[lower_indexes[k]], toolpath_size))
							{
								std::reverse(offsets[lower_indexes[k]].begin(), offsets[lower_indexes[k]].end());
							}
							int nearest_int_0;
							int nearest_int_1;
							Circuit::NearestPoint(offsets[upper_index], offsets[lower_indexes[k]], nearest_int_0, nearest_int_1);
							Circuit::MakeIndexAsFirst(offsets[upper_index], nearest_int_0);
							Circuit::MakeIndexAsFirst(offsets[lower_indexes[k]], nearest_int_1);
						}
					}
				}
			}
		}

		static bool GetSharingParts1(double toolpath_size, Vector3d1 &offset_0, Vector3d1 &offset_1,
			Vector3d2 &sharing_parts)
		{
			std::vector<bool> refer_bool(offset_0.size(), false);

			bool have_true = false;
			bool have_false = false;

			struct kdtree *tree = kd_create(3);

			for (int i = 0; i < offset_1.size(); i++)
			{
				void *val = &offset_1[i];
				kd_insert3(tree, offset_1[i][0], offset_1[i][1], offset_1[i][2],val);
			}

			//#pragma omp parallel 
			{
				for (int i = 0; i < offset_0.size(); i++)
				{
					//double d0 = Circuit::Distance(offset_0[i], offset_1);

					double *pos = new double[3];
					pos[0] = offset_0[i][0];
					pos[1] = offset_0[i][1];
					pos[2] = offset_0[i][2];

					struct kdres *r = kd_nearest(tree, pos);
					double position[3];
					*(int*)kd_res_item(r, position);

					double d = Strip::Distance(Vector3d(pos[0], pos[1], pos[2]), Vector3d(position[0], position[1], position[2]));

					//double d = cfs.GeodesicDistance(offset_0[i]);
					/************************************************************/
					/*
					double d = 1000000000.0;
					for (int j = 0; j < offset_1.size(); j++)
					{
					d = min(d, cfs.GeodesicDistance(offset_0[i], offset_1[j]));
					}*/
					/************************************************************/

					//I am so stupid to the current detection method. OMG

					refer_bool[i] = d < toolpath_size + toolpath_size*TOLERANCE;
					//refer_bool[i] = abs(d - toolpath_size) < toolpath_size*TOLERANCE;

					if (d < toolpath_size + toolpath_size*TOLERANCE)
					//if (abs(d - toolpath_size)<toolpath_size*TOLERANCE)
					{
						have_true = true;
					}
					else
					{
						have_false = true;
					}
				}
			}
			
			kd_free(tree);


			if (!have_true)
				return false;

			if (!have_false)
			{
				sharing_parts.push_back(offset_0);
				return true;
			}

			int start_index = -1;
			for (int i = 0; i < refer_bool.size(); i++)
			{
				if (!refer_bool[i] && refer_bool[(i + 1) % refer_bool.size()])
				{
					start_index = i;
					break;
				}
			}

			int index = (start_index + 1) % refer_bool.size();

			Vector3d1 one_sharing_part;
			do
			{
				if (refer_bool[index])
					one_sharing_part.push_back(offset_0[index]);
				else
				{
					if (Strip::GetTotalLength(one_sharing_part) > toolpath_size*2.0)
					{
						sharing_parts.push_back(one_sharing_part);
					}
					Vector3d1().swap(one_sharing_part);
				}

				index = (index + 1) % refer_bool.size();

			} while (index != start_index);


			if (Strip::GetTotalLength(one_sharing_part) > toolpath_size*2.0)
			{
				sharing_parts.push_back(one_sharing_part);
			}
			Vector3d1().swap(one_sharing_part);

			return 	(sharing_parts.size() > 0);
		}


		static bool GetSharingParts(double toolpath_size, Vector3d1 &offset_0, Vector3d1 &offset_1,
			Vector3d2 &sharing_parts)
		{
			std::vector<bool> refer_bool(offset_0.size(),false);

			bool have_true = false;
			bool have_false = false;


			#pragma omp parallel 
			{
				for (int i = 0; i < offset_0.size(); i++)
				{
					double d = Circuit::Distance(offset_0[i], offset_1);

					//double d = cfs.GeodesicDistance(offset_0[i]);
					/************************************************************/
					/*
					double d = 1000000000.0;
					for (int j = 0; j < offset_1.size(); j++)
					{
					d = min(d, cfs.GeodesicDistance(offset_0[i], offset_1[j]));
					}*/

					/************************************************************/

					refer_bool[i] = abs(d - toolpath_size) < toolpath_size*TOLERANCE;

					if (abs(d - toolpath_size)<toolpath_size*TOLERANCE)
					{
						have_true = true;
					}
					else
					{
						have_false = true;
					}
				}
			}



			if (!have_true)
				return false;

			if (!have_false)
			{
				sharing_parts.push_back(offset_0);
				return true;
			}
			
			int start_index = -1;
			for (int i = 0; i < refer_bool.size(); i++)
			{
				if (!refer_bool[i] && refer_bool[(i + 1) % refer_bool.size()])
				{
					start_index = i;
					break;
				}
			}

			int index = (start_index + 1) % refer_bool.size();

			Vector3d1 one_sharing_part;
			do
			{
				if (refer_bool[index])
					one_sharing_part.push_back(offset_0[index]);
				else
				{
					if (Strip::GetTotalLength(one_sharing_part) > toolpath_size*1.1)
					{
						sharing_parts.push_back(one_sharing_part);
					}
					Vector3d1().swap(one_sharing_part);
				}
				
				index = (index + 1) % refer_bool.size();
		
			} while (index!=start_index);


			if (Strip::GetTotalLength(one_sharing_part) > toolpath_size*1.1)
			{
				sharing_parts.push_back(one_sharing_part);
			}
			Vector3d1().swap(one_sharing_part);

			/*************************************************************************/
			if (false)
			if (sharing_parts.size() > 0)
			{
				for (int i = 0; i < sharing_parts.size(); i++)
				{
					Vector3d1  part = sharing_parts[i];

					std::vector<int> index;

					index.push_back(0);

					for (int j = 1; j < part.size() - 1; j++)
					{
						Vector3d p_0 = part[j - 1];
						Vector3d p_1 = part[j];
						Vector3d p_2 = part[j + 1];

						if (!isAlmostZero(getLength(p_1 - p_0)) && !isAlmostZero(getLength(p_2 - p_1)))
						{
							double angle = getAngleBetween(p_1 - p_0, p_2 - p_1);

							if (angle > MM_PI / 2.0)
							{
								index.push_back(j);
							}
						}
					}

					if (index.size() > 1)
					{
						index.push_back(part.size()-1);

						double max_length = -10000.0;
						for (int j = 0; j < index.size()-1; j++)
						{
							int index_0 = index[j];
							int index_1 = index[j+1];
							Vector3d1 one;
							for (int k = index_0; k <= index_1; k++)
							{
								one.push_back(part[k]);
							}
							double length = Strip::GetTotalLength(one);

							if (length > max_length)
							{
								max_length = length;
								sharing_parts[i] = one;
							}
						}
					}
				}
			}

			/*************************************************************************/

			return 	(sharing_parts.size() > 0);
		}

		static bool CheckManuallyGraph(std::vector<int> &manually_graph, int index_0,int index_1)
		{
			if (manually_graph.size() == 0) return true;

			for (int i = 0; i < manually_graph.size(); i = i + 2)
			{
				if ((index_0 == manually_graph[i] && index_1 == manually_graph[i + 1])||
					(index_1 == manually_graph[i] && index_0 == manually_graph[i + 1]))
				{
					return true;
				}
			}
			return false;
		}
	
		static  void BuildOffsetGraph_Nopar(string input_path, double toolpath_size, 
			Vector3d3 &offsetses, 
			Vector3d2 &offsets,
			std::vector<int> &offset_graph,
			Vector3d3 &offset_graph_sharing_parts, 
			std::vector<double> &parts_length)
		{

			////////////////////////////////////////////////////////
			std::vector<int> manually_graph;

			if (std::ifstream(input_path +"\\path\\manually.graph", std::ios::in))
			{
				std::ifstream file(input_path + "\\path\\manually.graph", std::ios::in);
				int nb;
				file >> nb;
				for (int i = 0; i < nb; i++){
					int x, y;
					file >> x;
					file >> y;
					manually_graph.push_back(x);
					manually_graph.push_back(y);
				}
			}
			////////////////////////////////////////////////////////

			DWORD start_time = GetTickCount();
			//DWORD end_time = GetTickCount();
			//std::cout << "[TIME] BuildOffsetGraph_Nopar: " << (end_time - start_time) / 1000.0 << std::endl;

			std::cout << "Build offsets graphs ..." << std::endl;

			//build offset graph
			std::vector<int> index_int;

			int index = -1;
			for (int i = 0; i < offsetses.size(); i++)
			{
				for (int j = 0; j < offsetses[i].size(); j++)
				{
					index++;
					index_int.push_back(i);
					index_int.push_back(j);
					index_int.push_back(index);
				}
			}

			//same laryer

			for (int i = 0; i < offsetses.size(); i++)
			{
				std::cout << i << " " << offsetses.size() << std::endl;

				for (int j = 0; j < offsetses[i].size()-1; j++)
				{
					int upper_index = GetOffsetIndex(index_int, i, j);

					for (int k = j+1; k < offsetses[i].size(); k++)
					{
						int lower_index = GetOffsetIndex(index_int, i, k);

						//if (upper_index == 62 && lower_index == 63)continue;
						//if (upper_index == 63 && lower_index == 62)continue;

						if (!CheckManuallyGraph(manually_graph, upper_index,lower_index))
						{
							continue;
						}


						//if (upper_index == 47 && lower_index == 42)continue;
						//if (upper_index == 42 && lower_index == 47)continue;

						//if (upper_index == 43 && lower_index == 49)continue;
						//if (upper_index == 49 && lower_index == 43)continue;

						//if (upper_index == 27 && lower_index == 22)continue;
						//if (upper_index == 22 && lower_index == 27)continue;

						//if (upper_index == 26 && lower_index == 21)continue;
						//if (upper_index == 21 && lower_index == 26)continue;


						//if (upper_index == 29 && lower_index == 27)continue;
						//if (upper_index == 27 && lower_index == 29)continue;

						//if (upper_index == 30 && lower_index == 24)continue;
						//if (upper_index == 24 && lower_index == 30)continue;

						//if (upper_index == 30 && lower_index == 32)continue;
						//if (upper_index == 32 && lower_index == 30)continue;

						//63 58

						//if (upper_index == 63 && lower_index == 58)continue;
						//if (upper_index == 58 && lower_index == 63)continue;

						Vector3d2 sharing_parts_0;
						Vector3d2 sharing_parts_1;

						bool b0 = GetSharingParts1(toolpath_size, offsets[upper_index], offsets[lower_index], sharing_parts_0);
						bool b1 = GetSharingParts1(toolpath_size, offsets[lower_index], offsets[upper_index], sharing_parts_1);

						if (MyGetUserName() == "debug")
						{
							if ((upper_index == 30 && lower_index == 32) || (lower_index == 30 && upper_index == 32))
							{
								std::ofstream debug_file("Z:\\Documents\\Windows\\SmartSFC\\workspace\\CFS\\sharingpart_" + IntString(upper_index) + "_" + IntString(lower_index) + ".obj");
								int debug_index = 1;

								for (int i = 0; i < sharing_parts_0.size(); i++)
								{
									for (int j = 0; j < sharing_parts_0[i].size() - 1; j++)
									{
										CGAL_Export_Path_Segment(debug_file, debug_index, "debug_0_" + IntString(upper_index) + "_" + IntString(lower_index), 1.0, 0.0, 0.0,
											sharing_parts_0[i][j], sharing_parts_0[i][j + 1], 0.1);
									}
								}
								for (int i = 0; i < sharing_parts_1.size(); i++)
								{
									for (int j = 0; j < sharing_parts_1[i].size() - 1; j++)
									{
										CGAL_Export_Path_Segment(debug_file, debug_index, "debug_1_" + IntString(upper_index) + "_" + IntString(lower_index), 1.0, 0.0, 0.0,
											sharing_parts_1[i][j], sharing_parts_1[i][j + 1], 0.1);
									}
								}
								debug_file.clear();
								debug_file.close();
							}

						}

						if (b0&&b1)
						{
							offset_graph.push_back(upper_index);
							offset_graph.push_back(lower_index);
							offset_graph_sharing_parts.push_back(sharing_parts_0);
							offset_graph_sharing_parts.push_back(sharing_parts_1);
						}

					}

				}
			}


			//upper and lower
			for (int i = 0; i < offsetses.size() - 1; i++)
			{
				std::cout << i << " " << offsetses.size() << std::endl;

				for (int j = 0; j < offsetses[i].size(); j++)
				{
					int upper_index = GetOffsetIndex(index_int, i, j);
					
					for (int k = 0; k < offsetses[i + 1].size(); k++)
					{
						int lower_index = GetOffsetIndex(index_int, i + 1, k);

						/*			if (upper_index == 62 && lower_index == 63)continue;
						if (upper_index == 63 && lower_index == 62)continue;

						if (upper_index == 47 && lower_index == 42)continue;
						if (upper_index == 42 && lower_index == 47)continue;

						if (upper_index == 43 && lower_index == 49)continue;
						if (upper_index == 49 && lower_index == 43)continue;

						if (upper_index == 27 && lower_index == 22)continue;
						if (upper_index == 22 && lower_index == 27)continue;

						if (upper_index == 26 && lower_index == 21)continue;
						if (upper_index == 21 && lower_index == 26)continue;


						if (upper_index == 29 && lower_index == 27)continue;
						if (upper_index == 27 && lower_index == 29)continue;

						if (upper_index == 30 && lower_index == 24)continue;
						if (upper_index == 24 && lower_index == 30)continue;


						if (upper_index == 30 && lower_index == 32)continue;
						if (upper_index == 32 && lower_index == 30)continue;*/


						if (!CheckManuallyGraph(manually_graph, upper_index, lower_index))
						{
							continue;
						}


						//if (upper_index == 63 && lower_index == 58)continue;
						//if (upper_index == 58 && lower_index == 63)continue;

						Vector3d2 sharing_parts_0;
						Vector3d2 sharing_parts_1;


						if ((upper_index == 30 && lower_index == 32) || (lower_index == 30 && upper_index == 32))
						{
							int dasdas = 10;
						}

						bool b0 = GetSharingParts1(toolpath_size, offsets[upper_index], offsets[lower_index], sharing_parts_0);
						bool b1 = GetSharingParts1(toolpath_size, offsets[lower_index], offsets[upper_index], sharing_parts_1);


						if (MyGetUserName() == "debug")
						{
							if ((upper_index == 30 && lower_index == 32) || (lower_index == 30 && upper_index == 32))
							{
								std::ofstream debug_file("Z:\\Documents\\Windows\\SmartSFC\\workspace\\CFS\\sharingpart_" + IntString(upper_index) + "_" + IntString(lower_index) + ".obj");
								int debug_index = 1;

								for (int i = 0; i < sharing_parts_0.size(); i++)
								{
									for (int j = 0; j < sharing_parts_0[i].size() - 1; j++)
									{
										CGAL_Export_Path_Segment(debug_file, debug_index, "debug_0_" + IntString(upper_index) + "_" + IntString(lower_index), 1.0, 0.0, 0.0,
											sharing_parts_0[i][j], sharing_parts_0[i][j + 1], 0.1);
									}
								}
								for (int i = 0; i < sharing_parts_1.size(); i++)
								{
									for (int j = 0; j < sharing_parts_1[i].size() - 1; j++)
									{
										CGAL_Export_Path_Segment(debug_file, debug_index, "debug_1_" + IntString(upper_index) + "_" + IntString(lower_index), 1.0, 0.0, 0.0,
											sharing_parts_1[i][j], sharing_parts_1[i][j + 1], 0.1);
									}
								}
								debug_file.clear();
								debug_file.close();
							}

						}



						if (b0&&b1)
						{
							offset_graph.push_back(upper_index);
							offset_graph.push_back(lower_index);
							offset_graph_sharing_parts.push_back(sharing_parts_0);
							offset_graph_sharing_parts.push_back(sharing_parts_1);
						}

					}
				}
			}

			std::vector<int>().swap(index_int);

			//assign weight length
			std::vector<int> node_degree = ComputeNodeDegree(offsets.size(), offset_graph);
			for (int i = 0; i < offset_graph.size(); i = i + 2)
			{
				if (node_degree[offset_graph[i]] > 2 || node_degree[offset_graph[i + 1]] > 2)
					parts_length.push_back(Circuit::GetOffsetPartsLength(offset_graph_sharing_parts[i]) + Circuit::GetOffsetPartsLength(offset_graph_sharing_parts[i + 1]));
				else
					parts_length.push_back(0.0);
			}
			std::vector<int>().swap(node_degree);

			DWORD end_time = GetTickCount();
			std::cout << "[TIME] BuildOffsetGraph_Nopar: " << (end_time - start_time) / 1000.0 << std::endl;
		}

		static void MinimalSpanningTree(Vector3d2 &offsets, std::vector<int> &offset_graph,
			std::vector<double> &parts_length, std::vector<int> &mst)
		{
			DWORD start_time = GetTickCount();
			//DWORD end_time = GetTickCount();
			//std::cout << "[TIME] MinimalSpanningTree: " << (end_time - start_time) / 1000.0 << std::endl;

			std::vector<int> nodes;
			std::vector<int> edges;
			std::vector<double> costs;

			//initialize nodes/edges/costs

			double maximal_length = 0.0;
			for (int i = 0; i < parts_length.size(); i++)
				maximal_length = std::max(maximal_length, parts_length[i]);

			for (int i = 0; i < offsets.size(); i++)
				nodes.push_back(i);

			for (int i = 0; i < offset_graph.size(); i = i + 2)
			{
				int index_0 = offset_graph[i];
				int index_1 = offset_graph[i + 1];
				edges.push_back(index_0);
				edges.push_back(index_1);
				costs.push_back(maximal_length - parts_length[i / 2]);
				//costs.push_back(parts_length[i / 2]);
			}

			std::vector<std::vector<int>> containers;
			for (int i = 0; i < nodes.size(); i++)
			{
				std::vector<int> container;
				container.push_back(nodes[i]);
				containers.push_back(container);
				std::vector<int>().swap(container);
			}

			std::vector<bool> edges_used;
			for (int i = 0; i < costs.size(); i++)
				edges_used.push_back(false);

			do
			{
				//find a minimal cost edge
				int minimal_cost_edge_index = -1;
				double minimal_cost = MAXDOUBLE;
				#pragma region find_a_minimal_cost_edge

				for (int j = 0; j < costs.size(); j++)
				{
					if (!edges_used[j])
					{
						if (costs[j] < minimal_cost)
						{
							minimal_cost = costs[j];
							minimal_cost_edge_index = j;
						}
					}
				}
				#pragma endregion

				if (minimal_cost_edge_index < 0)
					break;

				//check valid
				int node_index_0 = edges[2 * minimal_cost_edge_index];
				int node_index_1 = edges[2 * minimal_cost_edge_index + 1];


				if (node_index_0 == 16 && node_index_1 == 18)
				{
					int dsad = 0;
				}

				int container_0 = -1;
				int container_0_0 = -1;
				int container_1 = -1;
				int container_1_0 = -1;

				for (int j = 0; j < containers.size() && (container_0 < 0 || container_1 < 0); j++)
				{
					for (int k = 0; k < containers[j].size() && (container_0 < 0 || container_1 < 0); k++)
					{
						if (node_index_0 == containers[j][k])
						{
							container_0 = j;
							container_0_0 = k;
						}
						if (node_index_1 == containers[j][k])
						{
							container_1 = j;
							container_1_0 = k;
						}
					}
				}

				if (!(container_0 >= 0 && container_1 >= 0))
				{
					break;
				}

				if (container_0 == container_1)
				{
					edges_used[minimal_cost_edge_index] = true;
				}
				else
				{
					mst.push_back(node_index_0);
					mst.push_back(node_index_1);
					edges_used[minimal_cost_edge_index] = true;

					for (int i = 0; i < containers[container_1].size(); i++)
					{
						containers[container_0].push_back(containers[container_1][i]);
					}

					containers.erase(containers.begin() + container_1);
				}

			} while (containers.size() != 1);

			std::vector<bool>().swap(edges_used);
			std::vector<std::vector<int>>().swap(containers);

			std::vector<int>().swap(nodes);
			std::vector<int>().swap(edges);
			std::vector<double>().swap(costs);

			DWORD end_time = GetTickCount();
			std::cout << "[TIME] MinimalSpanningTree: " << (end_time - start_time) / 1000.0 << std::endl;
		}


		static std::vector<int> ComputeNodeDegree(int nodes_nb, std::vector<int> &edges)
		{
			std::vector<bool> node_used;
			std::vector<int> node_degree;

			for (int i = 0; i < nodes_nb; i++)
			{
				node_used.push_back(false);
				node_degree.push_back(0);
			}
			for (int i = 0; i < edges.size(); i = i + 2)
			{
				if (!node_used[edges[i]] && !node_used[edges[i + 1]])
				{
					node_degree[edges[i]]++;
					node_degree[edges[i + 1]]++;
				}
			}
			return node_degree;
		}

		//degree of next_node_id is 2
		static int AnotherSideNodeIdofDegree2Node(int node_nb, std::vector<int> edges, int start_node_id, int next_node_id)
		{
			for (int i = 0; i < edges.size(); i = i + 2)
			{
				if (edges[i] == next_node_id&&edges[i + 1] != start_node_id)
					return edges[i + 1];

				if (edges[i+1] == next_node_id&&edges[i] != start_node_id)
					return edges[i];
			}

			assert(false);
			return 0;
		}

		//Get one degree=2 sequence from start_node_id.
		//The next node is next_node_id.
		//The basic searching process is along the direction from start_node_id to next_node_id.
		//contour_nodes: all nodes
		//edges: indicate the edges among these nodes
		//start_node_id: start searching node id
		//next_node_id: next searching node id of this sequence
		static std::vector<int> GetNodeSequence(std::vector<Node> &contour_nodes, std::vector<int> edges, int start_node_id, int next_node_id)
		{
			std::vector<int> sequence;
			sequence.push_back(start_node_id);

			while (contour_nodes[next_node_id].degree == 2)
			{
				int next_id = AnotherSideNodeIdofDegree2Node(contour_nodes.size(), edges, start_node_id, next_node_id);
				sequence.push_back(next_node_id);
				start_node_id = next_node_id;
				next_node_id = next_id;
			}
			sequence.push_back(next_node_id);
			return sequence;
		}

		//Get all degree=2 sequences from start_node_id
		//For each sequence, generate one Fermat spiral
		//contour_nodes: all nodes
		//edges: indicate the edges among these nodes
		//start_node_id: start searching node id
		//node_sequences: output of sequences
		static void GetNodeSequence(std::vector<Node> &contour_nodes, std::vector<int> &edges, int  start_node_id, std::vector<std::vector<int>> &node_sequences)
		{
			for (int i = 0; i < edges.size(); i=i+2)
			{
				if (edges[i] == start_node_id&&!contour_nodes[edges[i + 1]].has_connected)
					node_sequences.push_back(GetNodeSequence(contour_nodes, edges, start_node_id, edges[i + 1]));
				if (edges[i + 1] == start_node_id&&!contour_nodes[edges[i]].has_connected)
					node_sequences.push_back(GetNodeSequence(contour_nodes, edges, start_node_id, edges[i]));
			}
			for (int i = 0; i < node_sequences.size(); i++)
			{
				for (int j = 0; j < node_sequences[i].size()-1; j++)
					contour_nodes[node_sequences[i][j]].has_connected = true;
			}
		}

		//Get all degree=2 sequences
		//For each sequence, generate one Fermat spiral
		//contour_nodes: all nodes
		//edges: indicate the edges among these nodes
		//node_sequences: output of all sequences
		//start_node_id: start searching node id
		static void GetAllSequences(std::vector<Node> &contour_nodes, std::vector<int> &edges, std::vector<std::vector<int>> &node_sequences, int  start_node_id=0)
		{
			int nb_0 = node_sequences.size();
			GetNodeSequence(contour_nodes, edges, start_node_id, node_sequences);
			int nb_1 = node_sequences.size();
			for (int i = nb_0; i < nb_1; i++)
			{
				if (contour_nodes[node_sequences[i][node_sequences[i].size() - 1]].degree>2)
					GetAllSequences(contour_nodes, edges, node_sequences, node_sequences[i][node_sequences[i].size() - 1]);
			}
		}

		static int Degree3Number(int node_nb, std::vector<int> tree)
		{
			std::vector<int> degree(node_nb,0);
			for (int i = 0; i < tree.size(); i = i + 2)
			{
				degree[tree[i]]++;
				degree[tree[i+1]]++;
			}
			int nb = 0;
			for (int i = 0; i < degree.size(); i++)
				if (degree[i] >= 3)
					nb++;
			return nb;
		}


		static std::vector<std::vector<int>> ConnectedComponents(const int &node_nb, const std::vector<int> &tree)
		{
			std::vector<std::vector<int>> components;
			std::vector<int> index(node_nb, -1);
			int nb = 0;
			for (int i = 0; i < tree.size(); i = i + 2)
			{
				if (index[tree[i]] == -1 && index[tree[i+1]] == -1)
				{
					index[tree[i]] = nb;
					index[tree[i+1]] = nb;
					nb++;
				}
				if (index[tree[i]] == -1 && index[tree[i + 1]] != -1)
				{
					index[tree[i]] = index[tree[i + 1]];
				}
				if (index[tree[i]] != -1 && index[tree[i + 1]] == -1)
				{
					index[tree[i + 1]] = index[tree[i]];
				}
			}

			for (int i = 0; i < nb; i++)
			{
				components.emplace_back(std::vector<int>());
				for (int j = 0; j < index.size(); j++)
					if (index[j] == i)
						components.back().emplace_back(j);
			}

			for (int j = 0; j < index.size(); j++) if (index[j] == -1)components.emplace_back(std::vector<int>(1,j));

			return components;
		}
	};
}
#endif
