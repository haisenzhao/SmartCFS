#include "stdafx.h"
#include "cgalpackage.h"


#include "NewtonApple_hull3D.h"
#include "Wm5MeshCurvature.h"
#include "Wm5IntrRay3Sphere3.h"
#include "Wm5BSplineCurveFit.h"

#include "Wm5ContBox2.h"
#include "Wm5Box2.h"
#include "Wm5Matrix2.h"

#include "clipper.hpp"
#include "kdtree.h"

#include <CGAL/subdivision_method_3.h>

#include<CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include<CGAL/Polygon_with_holes_2.h>
#include<CGAL/create_offset_polygons_from_polygon_with_holes_2.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/bounding_box.h>
#include <CGAL/barycenter.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/squared_distance_2.h>
#include <CGAL/intersections.h>
#include<CGAL/create_straight_skeleton_from_polygon_with_holes_2.h>
#include <CGAL/Segment_Delaunay_graph_2.h>
#include <CGAL/Segment_Delaunay_graph_filtered_traits_2.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Triangulation_2.h>

#include <CGAL/Polyhedron_incremental_builder_3.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>
#include <CGAL/AABB_halfedge_graph_segment_primitive.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/Polygon_mesh_slicer.h>
#include <CGAL/mesh_segmentation.h>
#include <CGAL/property_map.h>

#include <cstdlib>
#include <iterator>
#include <CGAL/Random.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Surface_mesh_shortest_path.h>
#include <CGAL/algorithm.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/AABB_triangle_primitive.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/Polygon_set_2.h>

typedef CGAL::Simple_cartesian<double> KC;
typedef CGAL::Segment_Delaunay_graph_filtered_traits_without_intersections_2<KC> Gt;
typedef CGAL::Segment_Delaunay_graph_2<Gt>  SDG2;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef K::FT FT;
typedef K::Point_2                    Point_2;
typedef K::Line_2                     Line_2;
typedef K::Segment_2                  Segment_2;
typedef K::Ray_2                      Ray_2;
typedef CGAL::Polygon_2<K>            Polygon_2;
typedef K::Point_3                    Point_3;
typedef K::Line_3                     Line_3;
typedef K::Segment_3                  Segment_3;
typedef K::Ray_3                      Ray_3;
typedef K::Plane_3					  Plane_3;
typedef K::Vector_3                   Vector_3;
typedef K::Sphere_3                   Sphere_3;
typedef CGAL::Polygon_set_2<K, std::vector<K::Point_2>> Polygon_set_2;

bool DebugInformation()
{
	return true;
}

struct FaceInfo2
{
	FaceInfo2(){}
	int nesting_level;
	bool in_domain(){
		return nesting_level % 2 == 1;
	}
};

Vector3d PointVector3d(Point_3 p)
{
	return Vector3d(p[0], p[1], p[2]);
}
Point_3 VectorPoint3d(Vector3d p)
{
	return Point_3(p[0], p[1], p[2]);
}

Vector2d PointVector2d(Point_2 p)
{
	return Vector2d(p[0], p[1]);
}
Point_2 VectorPoint2d(Vector2d p)
{
	return Point_2(p[0], p[1]);
}


// ----------------------- A CGAL::Vertex with decoration ------------------
template < class Gt, class Vb = CGAL::Triangulation_vertex_base_2<Gt> >
class Vertex : public  Vb {
	typedef Vb superclass;
public:
	typedef typename Vb::Vertex_handle      Vertex_handle;
	typedef typename Vb::Point              Point;

	template < typename TDS2 >
	struct Rebind_TDS {
		typedef typename Vb::template Rebind_TDS<TDS2>::Other Vb2;
		typedef Vertex<Gt, Vb2> Other;
	};

public:
	Vertex() : superclass() {}
	Vertex(const Point & p) : superclass(p) {}
	int index;
};

typedef CGAL::Triangulation_2<K>         Triangulation;
typedef Vertex<K>                      Vb;
typedef CGAL::Triangulation_face_base_with_info_2<FaceInfo2, K>    Fbb;
typedef CGAL::Constrained_triangulation_face_base_2<K, Fbb>        Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb>               TDS;
typedef CGAL::Exact_predicates_tag                                Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, TDS, Itag>  CDT;


//typedef CDT::Point                                                Point;
//typedef CGAL::Polygon_2<K>                                        Polygon_2;

typedef CGAL::Polyhedron_3<K, CGAL::Polyhedron_items_with_id_3> Polyhedron_3;
typedef Polyhedron_3::Facet_iterator Poly_facet_iterator;
typedef Polyhedron_3::Point_3 Poly_point_3;
typedef Polyhedron_3::HalfedgeDS Poly_halfedgeds_3;
typedef Polyhedron_3::Halfedge_handle Halfedge_handle;
typedef Polyhedron_3::Vertex_handle Vertex_handle;
typedef Polyhedron_3::Halfedge_around_vertex_const_circulator  HV_circulator;

typedef CGAL::Surface_mesh_shortest_path_traits<K, Polyhedron_3> Traits;
typedef CGAL::Surface_mesh_shortest_path<Traits> Surface_mesh_shortest_path;
typedef boost::graph_traits<Polyhedron_3> Graph_traits;
typedef Graph_traits::vertex_iterator vertex_iterator;
typedef Graph_traits::face_iterator face_iterator;

typedef CGAL::Surface_mesh<K::Point_3> Mesh;
typedef std::vector<K::Point_3> Polyline_type;
typedef std::list< Polyline_type > Polylines;
typedef CGAL::AABB_halfedge_graph_segment_primitive<Mesh> HGSP;
typedef CGAL::AABB_traits<K, HGSP>    AABB_traits;
typedef CGAL::AABB_tree<AABB_traits>  AABB_tree;


typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron_3> Primitive;
typedef CGAL::AABB_traits<K, Primitive> Traits_poly;
typedef CGAL::AABB_tree<Traits_poly> Tree;
typedef Tree::Point_and_primitive_id Point_and_primitive_id;
typedef boost::optional< Tree::Intersection_and_primitive_id<Ray_3>::Type> Ray_intersection;

typedef CGAL::AABB_face_graph_triangle_primitive<Mesh> Mesh_Primitive;
typedef CGAL::AABB_traits<K, Mesh_Primitive> Mesh_Traits;
typedef CGAL::AABB_tree<Mesh_Traits> Mesh_Tree;
typedef boost::optional<Mesh_Tree::Intersection_and_primitive_id<Ray_3>::Type> Mesh_Ray_intersection;


typedef KC::Triangle_3 Triangle_3;
typedef std::list<Triangle_3>::iterator Iterator_3;
typedef CGAL::AABB_triangle_primitive<KC, Iterator_3> Primitive_3;
typedef CGAL::AABB_traits<KC, Primitive_3> AABB_triangle_traits_3;
typedef CGAL::AABB_tree<AABB_triangle_traits_3> Tree_3;


// A modifier creating a triangle with the incremental builder.
template<class HDS>
class polyhedron_builder : public CGAL::Modifier_base<HDS> {
public:

	std::vector<double> coords;
	std::vector<int>    tris;

	polyhedron_builder(std::vector<double>& c, std::vector<int>& t){
		coords = c;
		tris = t;
	}

	void operator()(HDS& hds) {
		typedef typename HDS::Vertex   Vertex;
		typedef typename Vertex::Point Point;

		// create a cgal incremental builder
		CGAL::Polyhedron_incremental_builder_3<HDS> B(hds, true);
		B.begin_surface(coords.size() / 3, tris.size() / 3);
		//B.begin_surface(3, 1, 6);

		// add the polyhedron vertices
		for (int i = 0; i<(int)coords.size(); i += 3){
			B.add_vertex(Point(coords[i + 0], coords[i + 1], coords[i + 2]));
		}

		// add the polyhedron triangles
		for (int i = 0; i<(int)tris.size(); i += 3){
			B.begin_facet();
			B.add_vertex_to_facet(tris[i + 0]);
			B.add_vertex_to_facet(tris[i + 1]);
			B.add_vertex_to_facet(tris[i + 2]);
			B.end_facet();
		}
		//finish up the surface
		//B.end_surface();
	}
};


void
mark_domains(CDT& ct,
CDT::Face_handle start,
int index,
std::list<CDT::Edge>& border)
{
	if (start->info().nesting_level != -1){
		return;
	}
	std::list<CDT::Face_handle> queue;
	queue.push_back(start);
	while (!queue.empty()){
		CDT::Face_handle fh = queue.front();
		queue.pop_front();
		if (fh->info().nesting_level == -1){
			fh->info().nesting_level = index;
			for (int i = 0; i < 3; i++){
				CDT::Edge e(fh, i);
				CDT::Face_handle n = fh->neighbor(i);
				if (n->info().nesting_level == -1){
					if (ct.is_constrained(e)) border.push_back(e);
					else queue.push_back(n);
				}
			}
		}
	}
}
//explore set of facets connected with non constrained edges,
//and attribute to each such set a nesting level.
//We start from facets incident to the infinite vertex, with a nesting
//level of 0. Then we recursively consider the non-explored facets incident 
//to constrained edges bounding the former set and increase the nesting level by 1.
//Facets in the domain are those with an odd nesting level.
void
mark_domains(CDT& cdt)
{
	for (CDT::All_faces_iterator it = cdt.all_faces_begin(); it != cdt.all_faces_end(); ++it){
		it->info().nesting_level = -1;
	}
	std::list<CDT::Edge> border;
	mark_domains(cdt, cdt.infinite_face(), 0, border);
	while (!border.empty()){
		CDT::Edge e = border.front();
		border.pop_front();
		CDT::Face_handle n = e.first->neighbor(e.second);
		if (n->info().nesting_level == -1){
			mark_domains(cdt, n, e.first->info().nesting_level + 1, border);
		}
	}
}

Point_2 point_to_2d(const Point_3& p, Plane_3& pl)
{
	Vector_3 basis[2];
	//auto pop = pl.point();
	auto pop = pl.projection(Point_3(0.0, 0.0, 0.0));
	const Vector_3 vpop(pop.x(), pop.y(), pop.z());
	basis[0] = pl.base1() / CGAL::sqrt(pl.base1().squared_length());
	basis[1] = pl.base2() / CGAL::sqrt(pl.base2().squared_length());
	const Vector_3 ter(pop, p);
	auto x = ter * basis[0];
	auto y = ter * basis[1];
	return Point_2(x, y);
}

Point_3 point_to_3d(const Point_2& p, Plane_3& pl)
{
	Vector_3 basis[2];
	//auto pop = pl.point();
	auto pop = pl.projection(Point_3(0.0, 0.0, 0.0));
	const Vector_3 vpop(pop.x(), pop.y(), pop.z());
	basis[0] = pl.base1() / CGAL::sqrt(pl.base1().squared_length());
	basis[1] = pl.base2() / CGAL::sqrt(pl.base2().squared_length());
	Vector_3 nr(pl.a(), pl.b(), pl.c());
	const Point_3 vi = pop + (p.x() * basis[0] + p.y() * basis[1]);
	return vi;
}

void  Construct_Polyhedron(Polyhedron_3 &polyhedron, Vector3d1 &vecs, std::vector<int> &face_id_0, std::vector<int> &face_id_1, std::vector<int> &face_id_2)
{
	std::vector<double> coords;
	std::vector<int>    tris;
	for (int i = 0; i < vecs.size(); i++)
	{
		coords.push_back(vecs[i][0]);
		coords.push_back(vecs[i][1]);
		coords.push_back(vecs[i][2]);
	}
	for (int i = 0; i < face_id_0.size(); i++)
	{
		tris.push_back(face_id_0[i]);
		tris.push_back(face_id_1[i]);
		tris.push_back(face_id_2[i]);
	}

	polyhedron_builder<Poly_halfedgeds_3> builder(coords, tris);
	polyhedron.delegate(builder);

	CGAL::set_halfedgeds_items_id(polyhedron);
	std::size_t facet_id = 0;
	for (Polyhedron_3::Facet_iterator facet_it = polyhedron.facets_begin();
		facet_it != polyhedron.facets_end(); ++facet_it, ++facet_id) {
		facet_it->id() = facet_id;
	}

	std::vector<double>().swap(coords);
	std::vector<int>().swap(tris);
}

void  Construct_Polyhedron(Polyhedron_3 &polyhedron, std::string path)
{
	if (path.substr(path.size() - 3, path.size()) == "off")
	{
		std::ifstream input(path);
		input >> polyhedron;
		input.close();

		CGAL::set_halfedgeds_items_id(polyhedron);
		std::size_t facet_id = 0;
		for (Polyhedron_3::Facet_iterator facet_it = polyhedron.facets_begin();
			facet_it != polyhedron.facets_end(); ++facet_it, ++facet_id) {
			facet_it->id() = facet_id;
		}
	}
	if (path.substr(path.size() - 3, path.size()) == "obj")
	{
		Vector3d1 vecs;
		std::vector<int> face_id_0;
		std::vector<int> face_id_1;
		std::vector<int> face_id_2;
		CGAL_3D_Read_Triangle_Mesh(path, vecs, face_id_0, face_id_1, face_id_2);
		Construct_Polyhedron(polyhedron, vecs, face_id_0, face_id_1, face_id_2);
	}
}

void  Construct_Polyhedron(Polyhedron_3 &polyhedron, std::string path, 
	Vector3d1 &vecs, std::vector<int> &face_id_0, std::vector<int> &face_id_1, std::vector<int> &face_id_2)
{
	if (path.substr(path.size() - 3, path.size()) == "off")
	{
		std::ifstream input(path);
		input >> polyhedron;
		input.close();

		CGAL::set_halfedgeds_items_id(polyhedron);
		std::size_t facet_id = 0;
		for (Polyhedron_3::Facet_iterator facet_it = polyhedron.facets_begin();
			facet_it != polyhedron.facets_end(); ++facet_it, ++facet_id) {
			facet_it->id() = facet_id;
		}

		for (Polyhedron_3::Vertex_iterator iter = polyhedron.vertices_begin();iter != polyhedron.vertices_end(); iter++)
		{
			Poly_point_3 p = iter->point();
			vecs.push_back(Vector3d(p[0], p[1], p[2]));
		}

		for (Polyhedron_3::Face_iterator iter = polyhedron.facets_begin(); iter != polyhedron.facets_end(); iter++)
		{
			face_id_0.push_back(iter->halfedge()->next()->next()->vertex()->id());
			face_id_1.push_back(iter->halfedge()->vertex()->id());
			face_id_2.push_back(iter->halfedge()->next()->vertex()->id());
		}
	}
	if (path.substr(path.size() - 3, path.size()) == "obj")
	{
		CGAL_3D_Read_Triangle_Mesh(path, vecs, face_id_0, face_id_1, face_id_2);
		Construct_Polyhedron(polyhedron, vecs, face_id_0, face_id_1, face_id_2);
	}
}



/***************************************************************************************************/

//static const double DOUBLE_EPSILON = 1.0E-05;
//inline bool Math::IsAlmostZero(double value) {
//	return value < 10.0 * DOUBLE_EPSILON && value > -10.0 * DOUBLE_EPSILON;
//}

//all kinds of distance computing in 2D
/***************************************************************************************************/
double CGAL_2D_Distance_Point_Point(Vector2d p_0, Vector2d p_1)
{
	return sqrt((double)CGAL::squared_distance(Point_2(p_0[0], p_0[1]), Point_2(p_1[0], p_1[1])));
}

double CGAL_2D_Distance_Point_Segment(Vector2d v, Vector2d s_0, Vector2d s_1)
{
	return sqrt((double)CGAL::squared_distance(Point_2(v[0], v[1]), Segment_2(Point_2(s_0[0], s_0[1]), Point_2(s_1[0], s_1[1]))));
}

double CGAL_2D_Distance_Point_Line(double p_x, double p_y, double l_s_x, double l_s_y, double l_e_x, double l_e_y)
{
	return sqrt((double)CGAL::squared_distance(Point_2(p_x, p_y), Line_2(Point_2(l_s_x, l_s_y), Point_2(l_e_x, l_e_y))));
}
double CGAL_2D_Distance_Segment_Segment(double s_0_s_x, double s_0_s_y, double s_0_e_x, double s_0_e_y, double s_1_s_x, double s_1_s_y, double s_1_e_x, double s_1_e_y)
{
	return sqrt((double)CGAL::squared_distance(Segment_2(Point_2(s_0_s_x, s_0_s_y), Point_2(s_0_e_x, s_0_e_y)), Segment_2(Point_2(s_1_s_x, s_1_s_y), Point_2(s_1_e_x, s_1_e_y))));
}

void CGAL_2D_Projection_Point_Segment(double p_x, double p_y, double s_s_x, double s_s_y, double s_e_x, double s_e_y, double &o_x, double &o_y)
{
	Line_2 l(Point_2(s_s_x, s_s_y), Point_2(s_e_x, s_e_y));
	Point_2 m_p = l.projection(Point_2(p_x, p_y));

	double d_m_s = CGAL_2D_Distance_Point_Point(Vector2d(m_p[0], m_p[1]), Vector2d(s_s_x, s_s_y));
	double d_m_e = CGAL_2D_Distance_Point_Point(Vector2d(m_p[0], m_p[1]), Vector2d(s_e_x, s_e_y));
	double d_s_e = CGAL_2D_Distance_Point_Point(Vector2d(s_s_x, s_s_y), Vector2d(s_e_x, s_e_y));

	if (d_m_s >= d_s_e)
	{
		o_x = s_e_x;
		o_y = s_e_y;
		return;
	}

	if (d_m_e >= d_s_e)
	{
		o_x = s_s_x;
		o_y = s_s_y;
		return;
	}

	o_x = m_p[0];
	o_y = m_p[1];
}
/***************************************************************************************************/

//many kinds of distance computing in 3D
/***************************************************************************************************/
double CGAL_3D_Distance_Point_Point(double p_0_x, double p_0_y, double p_0_z, double p_1_x, double p_1_y, double p_1_z)
{
	return sqrt((double)CGAL::squared_distance(Point_3(p_0_x, p_0_y, p_0_z), Point_3(p_1_x, p_1_y, p_1_z)));
}

double CGAL_3D_Distance_Point_Point(Vector3d v0, Vector3d v1)
{
	return sqrt((double)CGAL::squared_distance(Point_3(v0[0], v0[1], v0[2]), Point_3(v1[0],v1[1],v1[2])));
}


//double CGAL_3D_Distance_Point_Segment(Vector3d p, Vector3d s_s, Vector3d s_e)
//{
//	return sqrt((double)CGAL::squared_distance(VectorPoint3d(p), Segment_3(VectorPoint3d(s_s), VectorPoint3d(s_e))));
//}

 double CGAL_3D_Distance_Point_Segment(Vector3d p, Vector3d s_s, Vector3d s_e)
{
	return sqrt((double)CGAL::squared_distance(VectorPoint3d(p), Segment_3(VectorPoint3d(s_s), VectorPoint3d(s_e))));
}

double CGAL_3D_Distance_Point_Line(Vector3d p, Vector3d l_s, Vector3d l_e)
{
	return sqrt((double)CGAL::squared_distance(VectorPoint3d(p), Line_3(VectorPoint3d(l_s), VectorPoint3d(l_e))));
}

Vector3d CGAL_3D_Projection_Point_Segment(Vector3d p, Vector3d s_s, Vector3d s_e)
{
	Line_3 l(VectorPoint3d(s_s), VectorPoint3d(s_e));
	Point_3 m_p = l.projection(VectorPoint3d(p));
	double d_m_s = CGAL_3D_Distance_Point_Point(m_p[0], m_p[1], m_p[2], s_s[0], s_s[1], s_s[2]);
	double d_m_e = CGAL_3D_Distance_Point_Point(m_p[0], m_p[1], m_p[2], s_e[0], s_e[1], s_e[2]);
	double d_s_e = CGAL_3D_Distance_Point_Point(s_s[0], s_s[1], s_s[2], s_e[0], s_e[1], s_e[2]);

	if (d_m_s >= d_s_e)
		return s_e;
	if (d_m_e >= d_s_e)
		return s_s;
	return PointVector3d(m_p);
}

Vector3d CGAL_3D_Projection_Point_Line(Vector3d p, Vector3d l_s, Vector3d l_e)
{
	Line_3 l(VectorPoint3d(l_s), VectorPoint3d(l_e));
	Point_3 o_p = l.projection(VectorPoint3d(p));
	return PointVector3d(o_p);
}

double CGAL_3D_Distance_Segment_Segment(Vector3d s_0_s, Vector3d s_0_e, Vector3d s_1_s, Vector3d s_1_e)
{
	return sqrt((double)CGAL::squared_distance(Segment_3(VectorPoint3d(s_0_s), VectorPoint3d(s_0_e)), Segment_3(VectorPoint3d(s_1_s), VectorPoint3d(s_1_e))));
}

double CGAL_3D_Distance_Point_Triangle(Vector3d p, Vector3d t_0, Vector3d t_1, Vector3d t_2)
{
	KC::Point_3 a(t_0[0], t_0[1], t_0[2]);
	KC::Point_3 b(t_1[0], t_1[1], t_1[2]);
	KC::Point_3 c(t_2[0], t_2[1], t_2[2]);

	std::list<Triangle_3> triangles;
	triangles.push_back(Triangle_3(a, b, c));

	Tree_3 tree(triangles.begin(), triangles.end());
	KC::Point_3 point_query(p[0], p[1], p[2]);

	return sqrt(tree.squared_distance(point_query));
}

double CGAL_3D_Distance_Point_Triangles(Vector3d p, Vector3d1& vecs, std::vector<int>& face_id_0, std::vector<int>& face_id_1, std::vector<int>& face_id_2)
{
	//build polyhedron
	Polyhedron_3 polyhedron;
	Construct_Polyhedron(polyhedron, vecs, face_id_0, face_id_1, face_id_2);

	//build tree
	Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);
	tree.accelerate_distance_queries();
	
	return sqrt(tree.squared_distance(Point_3(p[0], p[1], p[2])));
}

Vector3d CGAL_3D_Nearest_Point_Triangles(Vector3d p, Vector3d1& vecs, std::vector<int>& face_id_0, std::vector<int>& face_id_1, std::vector<int>& face_id_2)
{
	//build polyhedron
	Polyhedron_3 polyhedron;
	Construct_Polyhedron(polyhedron, vecs, face_id_0, face_id_1, face_id_2);

	//build tree
	Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);
	tree.accelerate_distance_queries();

	Point_3 c_p = tree.closest_point(Point_3(p[0], p[1], p[2]));
	return Vector3d(c_p[0],c_p[1],c_p[2]);
}

void CGAL_3D_Distance_Point_Mesh(std::string path, std::vector<double>& xs, std::vector<double>& ys, std::vector<double>& zs, std::vector<double>& ds)
{
	std::cout << "CGAL_3D_Distance_Point_Mesh" << std::endl;

	Polyhedron_3 polyhedron;
	Construct_Polyhedron(polyhedron, path);

	Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);
	tree.accelerate_distance_queries();


	for (int i = 0; i < xs.size(); i++)
	{
		if (i % (xs.size() / 10) == 0)
		{
			std::cout << (double)i / (double)xs.size() << std::endl;
		}
		ds.push_back(sqrt(tree.squared_distance(Point_3(xs[i], ys[i], zs[i]))));
	}
}

void CGAL_3D_Distance_Point_Mesh(std::string path, Vector3d1 &query_points, std::vector<double>& distances)
{
	std::cout << "CGAL_3D_Distance_Point_Mesh" << std::endl;

	Polyhedron_3 polyhedron;
	Construct_Polyhedron(polyhedron, path);

	Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);
	tree.accelerate_distance_queries();


	for (int i = 0; i < query_points.size(); i++)
	{
		if (i % (query_points.size() / 10) == 0)
		{
			std::cout << (double)i / (double)query_points.size() << std::endl;
		}
		distances.push_back(sqrt(tree.squared_distance(VectorPoint3d(query_points[i]))));
	}
}


void CGAL_3D_Neareast_Point_Mesh(std::string path, Vector3d1 &ves, Vector3d1 &ners)
{
	std::cout << "CGAL_3D_Distance_Point_Mesh" << std::endl;

	Polyhedron_3 polyhedron;
	Construct_Polyhedron(polyhedron, path);
	Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);
	tree.accelerate_distance_queries();

	for (int i = 0; i < ves.size(); i++)
	{
		if (i % (ves.size() / 10) == 0)
		{
			std::cout << (double)i / (double)ves.size() << std::endl;
		}
		Point_3 p = tree.closest_point(Point_3(ves[i][0], ves[i][1], ves[i][2]));

		ners.push_back(Vector3d(p[0], p[1], p[2]));
	}
}

double CGAL_3D_Distance_Point_Plane(Vector3d v, Vector3d plane_p, Vector3d plane_n)
{
	return sqrt((double)CGAL::squared_distance(VectorPoint3d(v), Plane_3(VectorPoint3d(plane_p), Vector_3(plane_n[0], plane_n[1], plane_n[2]))));
}


void CGAL_3D_Plane_3D_to_2D_Points(Vector3d &plane_p, Vector3d &plane_n,
	Vector3d1 &points_3d,
	Vector2d1 &points_2d) {
	Plane_3 plane(VectorPoint3d(plane_p), Vector_3(plane_n[0], plane_n[1], plane_n[2]));
	for (int i = 0; i < points_3d.size(); i++)
		points_2d.push_back(PointVector2d(point_to_2d(VectorPoint3d(points_3d[i]), plane)));
}

/***************************************************************************************************/

//search the neareast triangels within distance of d 
void CGAL_3D_Mesh_Near_Triangles(Vector3d1 &vecs, std::vector<int> &face_id_0, std::vector<int> &face_id_1, std::vector<int> &face_id_2,
	Vector3d1 &points, double d, std::vector<std::vector<int>> &triangles)
{
	//build polyhedron

	Polyhedron_3 polyhedron;
	Construct_Polyhedron(polyhedron, vecs, face_id_0, face_id_1, face_id_2);

	//build tree
	Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);
	tree.accelerate_distance_queries();

	for (int i = 0; i < points.size(); i++)
	{
		std::vector<int> triangle;

		Poly_point_3 query(points[i][0], points[i][1], points[i][2]);
		Point_and_primitive_id pp = tree.closest_point_and_primitive(query);

		std::priority_queue<Polyhedron_3::Facet_handle> facets;
		std::vector<int> save_index;
		facets.push(pp.second);
		save_index.push_back(pp.second->id());

		while (facets.size()!=0)
		{
			Polyhedron_3::Facet_handle fh = facets.top();
			triangle.push_back(fh->id());
			facets.pop();

			std::vector<Polyhedron_3::Facet_handle> neighbors;

			if (!fh->halfedge()->is_border_edge())
				neighbors.push_back(fh->halfedge()->opposite()->face());
			if (!fh->halfedge()->next()->is_border_edge())
				neighbors.push_back(fh->halfedge()->next()->opposite()->face());
			if (!fh->halfedge()->next()->next()->is_border_edge())
				neighbors.push_back(fh->halfedge()->next()->next()->opposite()->face());

			for (int j = 0; j < neighbors.size(); j++)
			{
				Polyhedron_3::Facet_handle n_fh = neighbors[j];
				std::vector<Poly_point_3> n_fh_vecs;
				n_fh_vecs.push_back(n_fh->halfedge()->vertex()->point());
				n_fh_vecs.push_back(n_fh->halfedge()->next()->vertex()->point());
				n_fh_vecs.push_back(n_fh->halfedge()->next()->next()->vertex()->point());

				bool add_bool = false;
				for (int k = 0; k < 3; k++)
				{
					double distance = sqrt((double)CGAL::squared_distance(pp.first, n_fh_vecs[k]));
					if (distance < d){
						add_bool = true;
						break;
					}
				}
				if (add_bool&&!(std::find(save_index.begin(), save_index.end(), n_fh->id()) != save_index.end()))
				{
					facets.push(n_fh);
					save_index.push_back(n_fh->id());
				}
			}
		}

		triangles.push_back(triangle);
	}
}

void CGAL_3D_Points_inside_Triangles(Vector3d1 &vecs, std::vector<int> &face_id_0, std::vector<int> &face_id_1, std::vector<int> &face_id_2,
	Vector3d1 &points, std::vector<bool> &insides)
{
	//build polyhedron
	Polyhedron_3 polyhedron;
	Construct_Polyhedron(polyhedron, vecs, face_id_0, face_id_1, face_id_2);

	CGAL::Side_of_triangle_mesh<Polyhedron_3, K> inside(polyhedron);

	for (int i = 0; i < points.size(); i++)
	{
		CGAL::Bounded_side res = inside(Point_3(points[i][0], points[i][1], points[i][2]));
		if (res == CGAL::ON_BOUNDED_SIDE)
			insides.push_back(true);
		else
			insides.push_back(false);
	}
}



void CGAL_3D_Points_inside_Triangles(std::string path, Vector3d1 &points, std::vector<bool> &insides)
{
	Polyhedron_3 polyhedron;
	Construct_Polyhedron(polyhedron, path);

	CGAL::Side_of_triangle_mesh<Polyhedron_3, K> inside(polyhedron);

	for (int i = 0; i < points.size(); i++)
	{
		CGAL::Bounded_side res = inside(Point_3(points[i][0], points[i][1], points[i][2]));
		if (res == CGAL::ON_BOUNDED_SIDE)
			insides.push_back(true);
		else
			insides.push_back(false);
	}
}

//Vector_3 RelatedFaceNormal(Polyhedron_3 &polyhedron, Tree &tree, double source_x, double source_y, double source_z)
//{
//	Poly_point_3 query(source_x, source_y, source_z);
//	Point_and_primitive_id pp = tree.closest_point_and_primitive(query);
//
//	Poly_point_3 p0 = pp.second->halfedge()->next()->next()->vertex()->point();
//	Poly_point_3 p1 = pp.second->halfedge()->vertex()->point();
//	Poly_point_3 p2 = pp.second->halfedge()->next()->vertex()->point();
//
//	Vector_3   n = CGAL::cross_product(p2 - p1, p0 - p1);
//
//	return n / std::sqrt(n*n);
//}

//many kinds of polygon computing in 2D
/***************************************************************************************************/

bool CGAL_2D_Detect_Polygon_Inside(Vector2d1 outside_py, Vector2d p)
{
	Vector2d outside_max, outside_min;
	CGAL_2d_Polygon_Boundingbox(outside_py, outside_min, outside_max);

	if (outside_max[0] < p[0])return false;
	if (outside_max[1] < p[1])return false;
	if (p[0] < outside_min[0])return false;
	if (p[1] < outside_min[1])return false;

	Polygon_2 poly;
	for (int i = 0; i < outside_py.size(); i++)
		poly.push_back(Point_2(outside_py[i][0], outside_py[i][1]));

	if (poly.bounded_side(Point_2(p[0], p[1])) == CGAL::ON_UNBOUNDED_SIDE)
	{
		double dis = CGAL_2D_Distance_Point_Polygon(p, outside_py);
		if (dis>0.1)
			return false;
	}

	return true;
}

bool CGAL_2D_Detect_Polygon_Inside(Vector2d1 outside_py, Vector2d1 inside_py)
{
	Vector2d outside_max, outside_min, inside_max, inside_min;
	CGAL_2d_Polygon_Boundingbox(outside_py, outside_min, outside_max);
	CGAL_2d_Polygon_Boundingbox(inside_py, inside_min, inside_max);

	if (outside_max[0] < inside_min[0])return false;
	if (outside_max[1] < inside_min[1])return false;
	if (inside_max[0] < outside_min[0])return false;
	if (inside_max[1] < outside_min[1])return false;

	Polygon_2 poly;
	for (int i = 0; i < outside_py.size(); i++)
		poly.push_back(Point_2(outside_py[i][0], outside_py[i][1]));

	for (auto p : inside_py)
	{
		if (poly.bounded_side(Point_2(p[0], p[1]))==CGAL::ON_UNBOUNDED_SIDE)
		{
			double dis = CGAL_2D_Distance_Point_Polygon(p, outside_py);
			if (dis>0.1)	
				return false;
		}
	}
	return true;
}

bool CGAL_2D_Detect_Polygon_Inside(Vector2d2 outside_pys, Vector2d pppp)
{
	if (!CGAL_2D_Detect_Polygon_Inside(outside_pys[0], pppp)) return false;

	for (int i = 1; i < outside_pys.size(); i++)
	{
		Polygon_2 poly;
		for (auto p : outside_pys[i]) poly.push_back(Point_2(p[0], p[1]));

		if (poly.bounded_side(Point_2(pppp[0], pppp[1])) == CGAL::ON_BOUNDED_SIDE)
		{
			double dis = CGAL_2D_Distance_Point_Polygon(pppp, outside_pys[i]);
			if (dis>0.1)
				return false;
		}
	}

	return true;
}

bool CGAL_2D_Detect_Polygon_Inside(Vector2d2 outside_pys, Vector2d1 inside_py)
{
	if (!CGAL_2D_Detect_Polygon_Inside(outside_pys[0],inside_py)) return false;

	for (int i = 1; i < outside_pys.size(); i++)
	{
		Polygon_2 poly;
		for (auto p : outside_pys[i]) poly.push_back(Point_2(p[0], p[1]));

		for (auto p : inside_py)
		{
			if (poly.bounded_side(Point_2(p[0], p[1])) == CGAL::ON_BOUNDED_SIDE)
			{
				double dis = CGAL_2D_Distance_Point_Polygon(p, outside_pys[i]);
				if (dis>0.1)
					return false;
			}
		}
	}

	return true;
}

bool CGAL_2D_Detect_Polygon_Inside(Vector2d2 outside_pys, Vector2d2 inside_pys)
{
	for (auto inside_py : inside_pys)
	{
		if (!CGAL_2D_Detect_Polygon_Inside(outside_pys, inside_py))
		{
			return false;
		}
	}
	return true;
}

//Check the location relationship between a point "p" and a 2d polygon "py"
//p: a 2d point
//py: a 2d polygon
//return true: p is in py
//return false: p is outside of py
bool CGAL_2D_Location_Point_Polygon(Vector2d p, Vector2d1 py)
{
	Polygon_2 poly;
	for (int i = 0; i < py.size(); i++)
		poly.push_back(Point_2(py[i][0], py[i][1]));
	return poly.bounded_side(Point_2(p[0], p[1])) == CGAL::ON_BOUNDED_SIDE;
}

bool CGAL_2D_Polygon_Is_Clockwise_Oriented(Vector2d1 &ps)
{
	Polygon_2 poly;
	for (int i = 0; i < ps.size(); i++)
		poly.push_back(Point_2(ps[i][0], ps[i][1]));

	return poly.is_clockwise_oriented();
}

double CGAL_2D_Polygon_Area(Vector2d1 py)
{
	Polygon_2 poly;
	for (int i = 0; i < py.size(); i++)
		poly.push_back(Point_2(py[i][0], py[i][1]));
	return abs(poly.area());
}

void CGAL_2D_Polygon_Dart_Sampling(Vector2d1 &py, double d, Vector2d1 &sampling_points)
{
	Polygon_2 poly;
	for (int i = 0; i < py.size(); i++)
		poly.push_back(Point_2(py[i][0], py[i][1]));

	double xmin = poly.bbox().xmin();
	double ymin = poly.bbox().ymin();
	double xmax = poly.bbox().xmax();
	double ymax = poly.bbox().ymax();

	int run = 0;
	Vector2d1 insert_points;

	while (run<10000)
	{
		run++;
		double x = rand() / double(RAND_MAX);
		double y = rand() / double(RAND_MAX);
		x = (xmax - xmin)*x + xmin;
		y = (ymax - ymin)*y + ymin;
		
		if (poly.bounded_side(Point_2(x, y)) == CGAL::ON_BOUNDED_SIDE)
		{
			double distance = CGAL_IA_MAX_DOUBLE;
			for (int i = 0; i < insert_points.size(); i++)
				distance = std::min(distance, CGAL_2D_Distance_Point_Point(insert_points[i], Vector2d(x, y)));

			//for (int i = 0; i < py.size(); i++)
			//distance = std::min(distance, CGAL_2D_Distance_Point_Point(py[i], Vector2d(x, y)));

			if (distance > d)
			{
				insert_points.push_back(Vector2d(x, y));
				run = 0;
			}
		}
	}


	for (int i = 0; i < insert_points.size(); i++)
	{
		double distance = CGAL_IA_MAX_DOUBLE;
		for (int j = 0; j < py.size(); j++)
			distance = std::min(distance, CGAL_2D_Distance_Point_Point(py[j], insert_points[i]));
		if (distance > d / 2.0)
			sampling_points.push_back(insert_points[i]);
	}

}

//void CGAL_2D_Polygon_Inside_Point(std::vector<double> &py_x, std::vector<double> &py_y, double &inside_x, double &inside_y)
//{

//}

double CGAL_2D_Distance_Point_Polygon(Vector2d p, Vector2d1 py) {
	double distance = 1000000000000.0;
	for (int i = 0; i < py.size(); i++)
		distance = std::min(distance, CGAL_2D_Distance_Point_Segment(p, py[i], py[(i + 1) % py.size()]));
	return distance;
}

double CGAL_2D_Distance_Point_Polygon(Vector2d p, Vector2d2 pys)
{
	double distance = 1000000000000.0;
	for (int i = 0; i < pys.size(); i++)
		distance = std::min(distance, CGAL_2D_Distance_Point_Polygon(p,pys[i]));
	return distance;
}

bool CGAL_2D_Polygon_Inside_Point(Vector2d2 polys, Vector2d &inner_vec)
{
	auto CheckValid = [](std::vector<Polygon_2> &cgal_polys, double p_x,double p_y)
	{
		if (cgal_polys.front().bounded_side(Point_2(p_x, p_y)) == CGAL::ON_BOUNDED_SIDE)
		{
			for (int i = 1; i < cgal_polys.size(); i++)
			{
				if (cgal_polys[i].bounded_side(Point_2(p_x, p_y)) == CGAL::ON_BOUNDED_SIDE)
				{
					return false;
				}
			}
			return true;
		}
		return false;
	};

	std::vector<Polygon_2> cgal_polys;
	for (int i = 0; i < polys.size(); i++)
	{
		Polygon_2 poly;
		for (int j = 0; j < polys[i].size(); j++)
			poly.push_back(VectorPoint2d(polys[i][j]));
		if (poly.is_clockwise_oriented()) poly.reverse_orientation();
		cgal_polys.emplace_back(poly);
	}

	auto outside_poly = cgal_polys.front();

	double xmin = outside_poly.bbox().xmin();
	double ymin = outside_poly.bbox().ymin();
	double xmax = outside_poly.bbox().xmax();
	double ymax = outside_poly.bbox().ymax();
	int inter = 0;
	int success_iter = 0;
	double dis_success;
	while (true)
	{
		double p_x = rand() / double(RAND_MAX)*(xmax - xmin) + xmin;
		double p_y = rand() / double(RAND_MAX)*(ymax - ymin) + ymin;

		if (CheckValid(cgal_polys,p_x,p_y))
		{
			if (success_iter == 0)
			{
				dis_success = CGAL_2D_Distance_Point_Polygon(Vector2d(p_x, p_y), polys);
				inner_vec[0] = p_x;
				inner_vec[1] = p_y;
			}
			else
			{
				//double distance = CGAL_2D_Distance_Point_Polygon(Vector2d(p_x, p_y), polys);
				double distance = p_y;
				if (distance > dis_success)
				{
					dis_success = distance;
					inner_vec[0] = p_x;
					inner_vec[1] = p_y;
				}
			}
			success_iter++;
		}

		if (inter > 10000 || success_iter>50 )
		{
			break;
		}
		inter++;
	}

	return success_iter!=0;
}

void CGAL_2d_Polygon_Boundingbox(Vector2d1 &ps, Vector2d &min_corner, Vector2d &max_corner)
{
	min_corner = ps[0];
	max_corner = ps[0];
	for (int i = 0; i < ps.size(); i++)
	{
		min_corner[0] = std::min(min_corner[0], ps[i][0]);
		min_corner[1] = std::min(min_corner[1], ps[i][1]);
		max_corner[0] = std::max(max_corner[0], ps[i][0]);
		max_corner[1] = std::max(max_corner[1], ps[i][1]);
	}
}

Vector2d CGAL_2D_Polygon_Inside_Point(Vector2d1 py)
{
	double inside_x, inside_y;
	Polygon_2 poly;
	for (int i = 0; i < py.size(); i++)
		poly.push_back(VectorPoint2d(py[i]));

	if (poly.is_clockwise_oriented()) poly.reverse_orientation();

	double max_dis = -10000.0;
	for (int i = 0; i < py.size(); i++){
		for (int j = i + 2; j < py.size(); j++){
			double p_x = (py[i][0] + py[j][0]) / 2.0;
			double p_y = (py[i][0] + py[j][0]) / 2.0;
			if (poly.bounded_side(Point_2(p_x, p_y)) == CGAL::ON_BOUNDED_SIDE)
			{
				double dis = CGAL_2D_Distance_Point_Point(py[i], py[j]);
				if (dis > max_dis)
				{
					max_dis = dis;
					inside_x = p_x;
					inside_y = p_y;
				}
			}
		}
	}
	return Vector2d(inside_x, inside_y);
}

static void RemoveClosePoints(std::vector<double> &xs, std::vector<double> &ys)
{
	std::vector<int> remove_int;
	if (xs.size() > 2)
	{
		for (int i = 0; i < xs.size() - 1; i++)
		{
			double d = CGAL_2D_Distance_Point_Point(Vector2d(xs[i], ys[i]), Vector2d(xs[(i + 1) % xs.size()], ys[(i + 1) % xs.size()]));
			if (d < 0.00001) remove_int.push_back((i + 1) % xs.size());
		}

		for (int i = remove_int.size() - 1; i >= 0; i--)
		{
			xs.erase(xs.begin() + remove_int[i]);
			ys.erase(ys.begin() + remove_int[i]);
		}
	}
}

void CGAL_2D_Polygon_Simple_0(Vector2d1 points_2d)
{
	for (int i = 0; i < points_2d.size(); i++)
	{
		//i//i+1
		auto s0 = points_2d[i];
		auto e0 = points_2d[(i + 1) % points_2d.size()];

		for (int j = i + 2; j < points_2d.size() + i - 1; j++)
		{
			auto s1 = points_2d[(j) % points_2d.size()];
			auto e1 = points_2d[(j + 1) % points_2d.size()];
			Vector2d inter;
			if (CGAL_2D_Intersection_Segment_Segment(s0, e0, s1, e1, inter))
			{
				std::cerr << i << " " << (j) % points_2d.size() << std::endl;
			}
		}
	}
}

bool CGAL_2D_Polygon_Simple(Vector2d2 poly)
{
	for (auto & py : poly)
		if (!CGAL_2D_Polygon_Simple(py))
			return false;
	return true;
}

bool CGAL_2D_Polygon_Simple(Vector2d1 py)
{
	Polygon_2 poly;
	for (int i = 0; i < py.size(); i++)
		poly.push_back(VectorPoint2d(py[i]));

	return poly.is_simple();
}

bool CGAL_2D_Polygon_Simple_Inter(const Vector2d1 &poly)
{
	for (int i = 0; i < poly.size(); i++)
	{
		auto s_0 = poly[i];
		auto e_0 = poly[(i + 1) % poly.size()];

		for (int j = 0; j < poly.size(); j++)
		{
			auto s_1 = poly[j];
			auto e_1 = poly[(j + 1) % poly.size()];

			if (i != j&&i != (j + 1) % poly.size() && (i + 1) % poly.size() != j && (i + 1) % poly.size() != (j + 1) % poly.size())
			{
				Vector2d inter;
				bool b = CGAL_2D_Intersection_Segment_Segment(s_0, e_0, s_1, e_1, inter);

				if (b)
				{
					return false;
				}
			}
		}
	}
	return true;
}



double CGAL_2D_Two_Polygons_Intersection(Vector2d1 poly_0, Vector2d1 poly_1)
{
	double scale = 1000000.0;

	ClipperLib::Paths subj(1);
	for (int i = 0; i < poly_0.size(); i++)
		subj[0] << ClipperLib::IntPoint(poly_0[i][0] * scale, poly_0[i][1] * scale);

	ClipperLib::Paths cliper(1);
	for (int i = 0; i < poly_1.size(); i++)
		cliper[0] << ClipperLib::IntPoint(poly_1[i][0] * scale, poly_1[i][1] * scale);

	ClipperLib::Paths solution;
	ClipperLib::Clipper c;
	c.AddPaths(subj, ClipperLib::ptSubject, true);
	c.AddPaths(cliper, ClipperLib::ptClip, true);
	c.Execute(ClipperLib::ctIntersection, solution, ClipperLib::pftNonZero, ClipperLib::pftNonZero);

	double area = 0.0;

	for (int i = 0; i < solution.size(); i++)
	{
		Polygon_2 poly_2;
		for (int j = 0; j < solution[i].size(); j++)
		{
			poly_2.push_back(Point_2(((double)solution[i][j].X) / scale, ((double)solution[i][j].Y) / scale));
		}
		area += poly_2.area();
	}

	return area;
}

double CGAL_2D_Two_Polygons_Intersection(Vector2d1 poly_0, Vector2d1 poly_1, Vector2d2 &inter_polygons)
{
	double scale = 1000000.0;

	ClipperLib::Paths subj(1);
	for (int i = 0; i < poly_0.size(); i++)
		subj[0] << ClipperLib::IntPoint(poly_0[i][0] * scale, poly_0[i][1] * scale);

	ClipperLib::Paths cliper(1);
	for (int i = 0; i < poly_1.size(); i++)
		cliper[0] << ClipperLib::IntPoint(poly_1[i][0] * scale, poly_1[i][1] * scale);

	ClipperLib::Paths solution;
	ClipperLib::Clipper c;
	c.AddPaths(subj, ClipperLib::ptSubject, true);
	c.AddPaths(cliper, ClipperLib::ptClip, true);
	c.Execute(ClipperLib::ctIntersection, solution, ClipperLib::pftNonZero, ClipperLib::pftNonZero);

	double area = 0.0;

	for (int i = 0; i < solution.size(); i++)
	{
		Polygon_2 poly_2;

		std::vector<double> xs;
		std::vector<double> ys;

		Vector2d1 polygon;
		
		for (int j = 0; j < solution[i].size(); j++)
		{
			poly_2.push_back(Point_2(((double)solution[i][j].X) / scale, ((double)solution[i][j].Y) / scale));
			polygon.push_back(Vector2d(((double)solution[i][j].X) / scale, ((double)solution[i][j].Y) / scale));
		}

		if (poly_2.area() > 0.0)
		{
			inter_polygons.push_back(polygon);
		}

		area += poly_2.area();
	}

	return area;
}

void insert_polygon(CDT& cdt, const Polygon_2& polygon, std::vector<int> &indexInt){
	if (polygon.is_empty()) return;
	int index = 0;

	CDT::Vertex_handle v_prev = cdt.insert(*CGAL::cpp11::prev(polygon.vertices_end()));

	for (Polygon_2::Vertex_iterator vit = polygon.vertices_begin();
		vit != polygon.vertices_end(); ++vit)
	{
		CDT::Vertex_handle vh = cdt.insert(*vit);
		vh->index = indexInt[index];
		index++;
		cdt.insert_constraint(vh, v_prev);
		v_prev = vh;
	}
}

std::vector<std::vector<int>> CGAL_2D_Polygon_Triangulation(const Vector2d1 &poly)
{
	Vector2d2 polys(1,poly);
	return CGAL_2D_Polygon_Triangulation(polys);
}

//bool CGAL_2D_Intersection_Segment_Segment(Vector2d s_0_s, Vector2d s_0_e, Vector2d s_1_s, Vector2d s_1_e, Vector2d &inter)

std::vector<std::vector<int>> CGAL_2D_Polygon_Triangulation(const Vector2d2 &polys)
{
	int nb = 0;
	CDT cdt;
	for (int i = 0; i < polys.size(); i++)
	{
		Polygon_2 polygon;
		std::vector<int> indexInt;
		for (int j = 0; j < polys[i].size(); j++)
		{
			polygon.push_back(Point_2(polys[i][j][0], polys[i][j][1]));
			indexInt.emplace_back(j+nb);
		}
		nb += polys[i].size();
		insert_polygon(cdt, polygon, indexInt);
	}

	//Mark facets that are inside the domain bounded by the polygon
	mark_domains(cdt);

	std::vector<std::vector<int>> faces;
	for (CDT::Finite_faces_iterator fit = cdt.finite_faces_begin();
		fit != cdt.finite_faces_end(); ++fit)
		if (fit->info().in_domain())
			faces.emplace_back(std::vector<int>{fit->vertex(2)->index, fit->vertex(1)->index, fit->vertex(0)->index});

	return faces;
}

void CGAL_2D_Polygon_Triangulation(Vector2d1 &p, std::string output_file, std::string path)
{
	//construct two non-intersecting nested polygons
	Vector3d1 vecs;
	Polygon_2 polygon;
	std::vector<int> indexInt;
	for (int i = 0; i < p.size(); i++)
	{
		polygon.push_back(Point_2(p[i][0], p[i][1]));
		indexInt.push_back(i);
		vecs.push_back(Vector3d(p[i][0], 0.0, p[i][1]));
	}

	//Insert the polygons into a constrained triangulation
	CDT cdt;
	insert_polygon(cdt, polygon, indexInt);
	//Mark facets that are inside the domain bounded by the polygon
	mark_domains(cdt);

	int count = 0;
	std::vector<int> face_id_0;
	std::vector<int> face_id_1;
	std::vector<int> face_id_2;
	for (CDT::Finite_faces_iterator fit = cdt.finite_faces_begin();
		fit != cdt.finite_faces_end(); ++fit)
	{
		if (fit->info().in_domain())
		{
			++count;
			face_id_0.push_back(fit->vertex(2)->index);
			face_id_1.push_back(fit->vertex(1)->index);
			face_id_2.push_back(fit->vertex(0)->index);
		}
	}
	//Output file as Obj or Off format
	if (output_file == "OFF" || output_file == "Off" || output_file == "off")
		CGAL_Output_Off(path,vecs,face_id_0, face_id_1, face_id_2);
	if (output_file == "OBJ" || output_file == "Obj" || output_file == "obj")
		CGAL_Output_Obj(path, vecs, face_id_0, face_id_1, face_id_2);
}


void CGAL_2D_Polygon_Triangulation(Vector2d1 &p, Vector2d1 &inside_p,
	std::vector<double> heights, std::string output_file, std::string path)
{
	//construct two non-intersecting nested polygons
	Vector3d1 vecs;
	Polygon_2 polygon;
	std::vector<int> indexInt;
	for (int i = 0; i < p.size(); i++)
	{
		polygon.push_back(Point_2(p[i][0], p[i][1]));
		indexInt.push_back(i);
		vecs.push_back(Vector3d(p[i][0], heights[i], p[i][1]));
	}

	//Insert the polygons into a constrained triangulation
	CDT cdt;
	insert_polygon(cdt, polygon, indexInt);

	for (int i = 0; i < inside_p.size(); i++)
	{
		vecs.push_back(Vector3d(inside_p[i][0], heights[i+p.size()], inside_p[i][1]));
		CDT::Vertex_handle vh = cdt.insert(Point_2(inside_p[i][0], inside_p[i][1]));
		vh->index = i + p.size();
	}

	//Mark facets that are inside the domain bounded by the polygon
	mark_domains(cdt);

	int count = 0;
	std::vector<int> face_id_0;
	std::vector<int> face_id_1;
	std::vector<int> face_id_2;
	for (CDT::Finite_faces_iterator fit = cdt.finite_faces_begin();
		fit != cdt.finite_faces_end(); ++fit)
	{
		if (fit->info().in_domain())
		{
			++count;
			face_id_0.push_back(fit->vertex(2)->index);
			face_id_1.push_back(fit->vertex(1)->index);
			face_id_2.push_back(fit->vertex(0)->index);
		}
	}

	//Output file as Obj or Off format
	if (output_file == "OFF" || output_file == "Off" || output_file == "off")
		CGAL_Output_Off(path, vecs, face_id_0, face_id_1, face_id_2);
	if (output_file == "OBJ" || output_file == "Obj" || output_file == "obj")
		CGAL_Output_Obj(path, vecs, face_id_0, face_id_1, face_id_2);
}



void CGAL_2D_Polygon_One_Offsets(const std::vector<Vector2d>& xys, double d, std::vector<std::vector<Vector2d>> &offsets_xys)
{
	CGAL_2D_Polygon_One_Offsets(std::vector<std::vector<Vector2d>>(1,xys), d, offsets_xys);
}

void CGAL_2D_Polygon_One_Offsets(const std::vector<std::vector<Vector2d>> &xys, double d, std::vector<std::vector<Vector2d>> &offsets_xys)
{
	std::vector<std::vector<double>> xs, ys;
	for (int i = 0; i < xys.size(); i++)
	{
		xs.emplace_back(std::vector<double>());
		ys.emplace_back(std::vector<double>());
		for (int j = 0; j < xys[i].size(); j++)
		{
			xs.back().emplace_back(xys[i][j][0]);
			ys.back().emplace_back(xys[i][j][1]);
		}
	}

	std::vector<std::vector<double>> offsets_xs, offsets_ys;
	CGAL_2D_Polygon_One_Offsets(xs, ys, d, offsets_xs, offsets_ys);

	for (int i = 0; i < offsets_xs.size(); i++)
	{
		offsets_xys.emplace_back(std::vector<Vector2d>());
		for (int j = 0; j < offsets_xs[i].size(); j++)
			offsets_xys.back().emplace_back(offsets_xs[i][j], offsets_ys[i][j]);
	}
}
void CGAL_2D_Polygon_One_Offsets(std::vector<std::vector<double>> xs, std::vector<std::vector<double>> ys, double d,
	std::vector<std::vector<double>> &offsets_xs, std::vector<std::vector<double>> &offsets_ys)
{
	CGAL_2D_Polygon_Offsets(xs, ys, d, offsets_xs, offsets_ys);
}

void CGAL_2D_Polygon_Offsets(std::vector<std::vector<double>> xs, std::vector<std::vector<double>> ys, double d,
	std::vector<std::vector<double>> &offsets_xs, std::vector<std::vector<double>> &offsets_ys)
{
	if (!(xs.size()>0 && xs[0].size()>0)) return;

	double scale = 1000000.0;

	ClipperLib::ClipperOffset co;
	co.ArcTolerance = co.ArcTolerance*scale / 1000.0;

	ClipperLib::Path subj;
	ClipperLib::Paths solution;

	//build the most outside path
	for (int i = 0; i < xs[0].size(); i++)
		subj << ClipperLib::IntPoint(xs[0][i] * scale, ys[0][i] * scale);
	co.AddPath(subj, ClipperLib::jtRound, ClipperLib::etClosedPolygon);

	//build the following pathes
	for (int i = 1; i < xs.size(); i++)
	{
		subj.clear();
		for (int j = 0; j < xs[i].size(); j++)
			subj << ClipperLib::IntPoint(xs[i][j] * scale, ys[i][j] * scale);
		co.AddPath(subj, ClipperLib::jtRound, ClipperLib::etClosedPolygon);
	}

	// execute
	co.Execute(solution, -d*scale);

	//output
	for (int i = 0; i < solution.size(); i++)
	{
		std::vector<double> one_offset_xs;
		std::vector<double> one_offset_ys;

		Polygon_2 poly_2;
		for (int j = 0; j < solution[i].size(); j++)
		{
			one_offset_xs.push_back(((double)solution[i][j].X) / scale);
			one_offset_ys.push_back(((double)solution[i][j].Y) / scale);

			poly_2.push_back(Point_2(((double)solution[i][j].X) / scale, ((double)solution[i][j].Y) / scale));
		}

		if (poly_2.is_clockwise_oriented())
		{
			std::reverse(one_offset_xs.begin(), one_offset_xs.end());
			std::reverse(one_offset_ys.begin(), one_offset_ys.end());
		}

		//remove closed points
		RemoveClosePoints(one_offset_xs, one_offset_ys);

		offsets_xs.push_back(one_offset_xs);
		offsets_ys.push_back(one_offset_ys);

		std::vector<double>().swap(one_offset_xs);
		std::vector<double>().swap(one_offset_ys);
	}
}

double CGAL_2D_Distance_Polygon_Polygon(Vector2d1 poly_0, Vector2d1 poly_1)
{
	double min_d = 1000000000.0;
	for (int i = 0; i < poly_0.size(); i++)
	{
		Vector2d v = CGAL_2D_Nearest_Point_Polygon(poly_0[i], poly_1);
		double l = CGAL_2D_Distance_Point_Point(poly_0[i],v);
		if (l < min_d)
		{
			min_d = l;
		}
	}
	return min_d;
}

double CGAL_2D_Distance_Polygon_Polygon(Vector2d2 poly_0, Vector2d2 poly_1)
{
	double min_d = 1000000000.0;
	for (auto poly_ : poly_0)
	{
		for (auto poly__ : poly_1)
		{
			double l = CGAL_2D_Distance_Polygon_Polygon(poly_, poly__);
			if (l < min_d)
			{
				min_d = l;
			}
		}
	}
	return min_d;
}

Vector2d CGAL_2D_Nearest_Point_Polygon(Vector2d v, Vector2d1 poly)
{
	double min_d = 1000000000.0;
	int min_i = -1;

	for (int i = 0; i < poly.size(); i++)
	{
		double l = CGAL_2D_Distance_Point_Segment(v, poly[i], poly[(i+1)%poly.size()]);

		if (l < min_d)
		{
			min_d = l;
			min_i = i;
		}
	}

	double o_x, o_y;
	CGAL_2D_Projection_Point_Segment(v[0], v[1], poly[min_i][0], poly[min_i][1], poly[(min_i + 1) % poly.size()][0], poly[(min_i + 1) % poly.size()][1], o_x, o_y);
	return Vector2d(o_x, o_y);
}

void CGAL_2D_Nearest_Point_Polygon(Vector2d v, Vector2d1 poly, Vector2d &p, double &min_d)
{
	min_d = 1000000000.0;
	int min_i = -1;

	for (int i = 0; i < poly.size(); i++)
	{
		double l = CGAL_2D_Distance_Point_Segment(v, poly[i], poly[(i + 1) % poly.size()]);

		if (l < min_d)
		{
			min_d = l;
			min_i = i;
		}
	}

	double o_x, o_y;
	CGAL_2D_Projection_Point_Segment(v[0], v[1], poly[min_i][0], poly[min_i][1], poly[(min_i + 1) % poly.size()][0], poly[(min_i + 1) % poly.size()][1], o_x, o_y);
	p[0] = o_x;
	p[1] = o_y;
}

Vector2d CGAL_2D_Nearest_Point_Polygon(Vector2d v, Vector2d2 polys)
{
	Vector2d result;
	double min_d = 1000000000.0;
	for (int i = 0; i < polys.size(); i++)
	{
		Vector2d p;
		double p_d;
		CGAL_2D_Nearest_Point_Polygon(v, polys[i], p, p_d);

		if (p_d < min_d)
		{
			min_d = p_d;
			result = p;
		}
	}
	return result;
}

/***************************************************************************************************/


//many kinds of intersection in 2D
/***************************************************************************************************/
bool CGAL_2D_Intersection_Segment_Line(Vector2d s_s, Vector2d s_e, Vector2d l_s, Vector2d l_e, Vector2d &inter)
{
	CGAL::Object result = CGAL::intersection(Segment_2(Point_2(s_s[0], s_s[1]), Point_2(s_e[0], s_e[1])),
		Line_2(Point_2(l_s[0], l_s[1]), Point_2(l_e[0], l_e[1])));

	if (const Point_2 *ipoint = CGAL::object_cast<Point_2>(&result))
	{
		inter[0] = ipoint->x();
		inter[1] = ipoint->y();
		return true;
	}
	else
	{
		return false;
	}
}

bool CGAL_2D_Intersection_Segment_Segment(Vector2d s_0_s, Vector2d s_0_e, Vector2d s_1_s, Vector2d s_1_e, Vector2d &inter)
{
	CGAL::Object result = CGAL::intersection(Segment_2(Point_2(s_0_s[0], s_0_s[1]), Point_2(s_0_e[0], s_0_e[1])),
		Segment_2(Point_2(s_1_s[0], s_1_s[1]), Point_2(s_1_e[0], s_1_e[1])));

	if (const Point_2 *ipoint = CGAL::object_cast<Point_2>(&result))
	{
		inter[0] = ipoint->x();
		inter[1] = ipoint->y();
		return true;
	}
	else
	{
		return false;
	}
}

bool CGAL_2D_Intersection_Line_Line(double s_s_x, double s_s_y, double s_e_x, double s_e_y
	, double l_s_x, double l_s_y, double l_e_x, double l_e_y, double& i_x, double& i_y)
{
	CGAL::Object result = CGAL::intersection(Line_2(Point_2(s_s_x, s_s_y), Point_2(s_e_x, s_e_y)),
		Line_2(Point_2(l_s_x, l_s_y), Point_2(l_e_x, l_e_y)));

	if (const Point_2 *ipoint = CGAL::object_cast<Point_2>(&result))
	{
		i_x = ipoint->x();
		i_y = ipoint->y();
		return true;
	}
	else
	{
		return false;
	}
}

bool CGAL_2D_Intersection_Segment_Polygon(Vector2d s_s, Vector2d s_e, Vector2d1 &p)
{
	for (int i = 0; i < p.size(); i++)
	{
		Vector2d inter;
		if (CGAL_2D_Intersection_Segment_Segment(s_s, s_e, p[i], p[(i + 1) % p.size()], inter))
		{
			return true;
		}
	}

	return false;

	/*
	Polygon_2 poly_2;
	for (int i = 0; i < p.size(); i++)
		poly_2.push_back(Point_2(p[i][0], p[i][1]));

	Polygon_set_2 ps;
	Polygon_2     line; // line is a polygon defined by 2 points
	line.push_back(Point_2(s_s[0], s_s[1]));
	line.push_back(Point_2(s_e[0], s_e[1]));

	ps.insert(poly_2);
	return ps.do_intersect(line);
	*/

}

/***************************************************************************************************/


//many kinds of intersection in 3D
/***************************************************************************************************/

bool CGAL_3D_Intersection_Segment_Line(Vector3d s_s, Vector3d s_e, Vector3d l_s, Vector3d l_e, Vector3d& inter)
{
	double d = CGAL_3D_Distance_Point_Point(s_s[0], s_s[1], s_s[2], s_e[0], s_e[1], s_e[2]);

	if (!Math::IsAlmostZero(d))
	{
		CGAL::Object result = CGAL::intersection(Segment_3(VectorPoint3d(s_s), VectorPoint3d(s_e)),
			Line_3(VectorPoint3d(l_s), VectorPoint3d(l_e)));

		if (const Point_3 *ipoint = CGAL::object_cast<Point_3>(&result))
		{
			inter[0] = ipoint->x();
			inter[1] = ipoint->y();
			inter[2] = ipoint->z();
			return true;
		}
		else
		{
			return false;
		}
	}
	else
	{
		return false;
	}
}



bool CGAL_3D_Intersection_Segment_Segment(Vector3d s_0_s, Vector3d s_0_e, Vector3d s_1_s, Vector3d s_1_e,Vector3d &iter)
{
	double d0 = CGAL_3D_Distance_Point_Point(s_0_s[0], s_0_s[1], s_0_s[2], s_0_e[0], s_0_e[1], s_0_e[2]);
	double d1 = CGAL_3D_Distance_Point_Point(s_1_s[0], s_1_s[1], s_1_s[2], s_1_e[0], s_1_e[1], s_1_e[2]);

	if (!Math::IsAlmostZero(d0) && !Math::IsAlmostZero(d1))
	{
		CGAL::Object result = CGAL::intersection(Segment_3(VectorPoint3d(s_0_s), VectorPoint3d(s_0_e)),
			Segment_3(VectorPoint3d(s_1_s), VectorPoint3d(s_1_e)));

		if (const Point_3 *ipoint = CGAL::object_cast<Point_3>(&result))
		{
			iter = PointVector3d(*ipoint);
			return true;
		}
		else
		{
			return false;
		}
	}
	else
	{
		return false;
	}
}

bool CGAL_3D_Intersection_Segment_Plane(Vector3d s_s, Vector3d s_e, Vector3d plane_p, Vector3d plane_n, Vector3d& inter)
{
	Segment_3 s(VectorPoint3d(s_s), VectorPoint3d(s_e));
	Plane_3 p(VectorPoint3d(plane_p), Vector_3(plane_n[0], plane_n[1], plane_n[2]));
	CGAL::Object result = CGAL::intersection(s, p);
	if (const Point_3 *ipoint = CGAL::object_cast<Point_3>(&result))
	{
		inter[0] = ipoint->x();
		inter[1] = ipoint->y();
		inter[2] = ipoint->z();
		return true;
	}
	else
	{
		return false;
	}
}


bool CGAL_3D_Intersection_Line_Plane(Vector3d l_s, Vector3d l_e, Vector3d plane_p, Vector3d plane_n, Vector3d& inter)
{
	Line_3 s(VectorPoint3d(l_s), VectorPoint3d(l_e));
	Plane_3 p(VectorPoint3d(plane_p), Vector_3(plane_n[0], plane_n[1], plane_n[2]));

	CGAL::Object result = CGAL::intersection(s, p);

	if (const Point_3 *ipoint = CGAL::object_cast<Point_3>(&result))
	{
		inter[0] = ipoint->x();
		inter[1] = ipoint->y();
		inter[2] = ipoint->z();
		return true;
	}
	else
	{
		return false;
	}
}


bool CGAL_3D_Intersection_Sphere_Ray(double center_x, double center_y, double center_z, double radius,
	double ray_origin_x, double ray_origin_y, double ray_origin_z, double ray_direction_x, double ray_direction_y, double ray_direction_z,
	std::vector<double>& i_x, std::vector<double>& i_y, std::vector<double>& i_z)
{
	Wm5::Sphere3<double> sphere(Wm5::Vector3d(center_x, center_y, center_z), radius);
	Wm5::Ray3<double> ray(Wm5::Vector3d(ray_origin_x, ray_origin_y, ray_origin_z), Wm5::Vector3d(ray_direction_x, ray_direction_y, ray_direction_z));

	Wm5::IntrRay3Sphere3d intr(ray,sphere);

	intr.Test();
	intr.Find();

	int nb = intr.GetQuantity();

	for (int i = 0; i < nb; i++)
	{
		Wm5::Vector3d p = intr.GetPoint(i);
		i_x.push_back(p[0]);
		i_y.push_back(p[1]);
		i_z.push_back(p[2]);
	}

	return nb>0;
}

bool CGAL_3D_Intersection_Ray_Triangle(Vector3d p, Vector3d n, Vector3d p0, Vector3d p1, Vector3d p2)
{
	Ray_3 ray(VectorPoint3d(p), Vector_3(n[0], n[1], n[2]));

	K::Triangle_3 tri(VectorPoint3d(p0), VectorPoint3d(p1), VectorPoint3d(p2));
	CGAL::Object result = CGAL::intersection(ray, tri);
	if (const Point_2 *ipoint = CGAL::object_cast<Point_2>(&result))
	{
		return true;
	}
	else
	{
		return false;
	}
}


bool CGAL_3D_Intersection_Ray_Mesh(Vector3d p, Vector3d n, std::string path)
{
	std::cout << "CGAL_3D_Intersection_Ray_Mesh..." << std::endl;

	Polyhedron_3 polyhedron;
	Construct_Polyhedron(polyhedron, path);

	Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);
	tree.accelerate_distance_queries();

	Ray_3 ray(Point_3(p[0], p[1], p[2]), Vector_3(n[0], n[1], n[2]));
	
	if (tree.do_intersect(ray))
		return true;
	else
		return false;
}

void CGAL_3D_Intersection_Ray_Mesh(Vector3d1 ps, Vector3d1 ns, std::string path, Vector3d1 &inters)
{

	std::ifstream input(path);
	Mesh mesh;
	input >> mesh;
	Mesh_Tree tree(faces(mesh).first, faces(mesh).second, mesh);

	for (int i = 0; i < ps.size(); i++)
	{
		Ray_3 ray(Point_3(ps[i][0], ps[i][1], ps[i][2]), Vector_3(ns[i][0], ns[i][1], ns[i][2]));

		std::list<Mesh_Ray_intersection> intersections;
		tree.all_intersections(ray, std::back_inserter(intersections));
		
		bool b = false;

		double min_d = 100000000000.0;

		Vector3d near_p;
		for (auto iter = intersections.begin(); iter != intersections.end(); iter++)
		{
			const Point_3* p = boost::get<Point_3>(&(iter->value().first));

			if (p)
			{
				Vector3d v(p->x(), p->y(), p->z());
				double d = CGAL_3D_Distance_Point_Point(v, ps[i]);

				if (min_d > d)
				{
					min_d = d;
					near_p = v;
				}

				b = true;
			}

		}

	

		if (!b)
		{
			inters.push_back(ps[i]);
		}
		else
		{
			inters.push_back(near_p);
		}

	}
}

/***************************************************************************************************/


//subdivision the input mesh
/***************************************************************************************************/
void CGAL_Mesh_Subdivision(std::string in_path, std::string sub, int step, std::string out_path)
{
	//load the input obj
	Polyhedron_3 polyhedron;
	Construct_Polyhedron(polyhedron, in_path);
	//subdivision

	if (sub == "Loop" || sub == "loop" || sub == "l" || sub == "L")
		CGAL::Subdivision_method_3::Loop_subdivision(polyhedron, step);

	if (sub == "Sqrt" || sub == "sqrt" || sub == "s" || sub == "S")
		CGAL::Subdivision_method_3::Sqrt3_subdivision(polyhedron, step);

	Vector3d1 vecs;
	std::vector<int> face_id_0;
	std::vector<int> face_id_1;
	std::vector<int> face_id_2;

	for (Polyhedron_3::Vertex_iterator iter = polyhedron.vertices_begin();
		iter != polyhedron.vertices_end(); iter++)
	{
		Poly_point_3 p = iter->point();
		vecs.push_back(Vector3d(p[0], p[1], p[2]));
	}

	for (Polyhedron_3::Face_iterator iter = polyhedron.facets_begin(); iter != polyhedron.facets_end(); iter++)
	{
		face_id_0.push_back(iter->halfedge()->next()->next()->vertex()->id());
		face_id_1.push_back(iter->halfedge()->vertex()->id());
		face_id_2.push_back(iter->halfedge()->next()->vertex()->id());
	}

	CGAL_Output_Obj(out_path, vecs, face_id_0, face_id_1, face_id_2);
}

void CGAL_Mesh_Loop_Subdivision_Own_Version(std::string in_path, int step, std::string out_path, int laplace_nb)
{
	Vector3d1 vecs;
	std::vector<int> face_id_0;
	std::vector<int> face_id_1;
	std::vector<int> face_id_2;

	CGAL_3D_Read_Triangle_Mesh(in_path, vecs, face_id_0, face_id_1, face_id_2);

	for (int i = 0; i < step; i++)
	{
		CGAL_Mesh_Loop_Subdivision_One_Step(vecs, face_id_0, face_id_1, face_id_2);
		CGAL_Mesh_Laplace_Smooth(vecs, face_id_0, face_id_1, face_id_2, laplace_nb);
	}

	CGAL_Output_Obj(out_path, vecs, face_id_0, face_id_1, face_id_2);

	Vector3d1().swap(vecs);
	std::vector<int>().swap(face_id_0);
	std::vector<int>().swap(face_id_1);
	std::vector<int>().swap(face_id_2);
}

void CGAL_Mesh_Laplace_Smooth(std::string in_path, std::string out_path, int laplace_nb)
{
	Vector3d1 vecs;
	std::vector<int> face_id_0;
	std::vector<int> face_id_1;
	std::vector<int> face_id_2;
	CGAL_3D_Read_Triangle_Mesh(in_path, vecs, face_id_0, face_id_1, face_id_2);
	CGAL_Mesh_Laplace_Smooth(vecs, face_id_0, face_id_1, face_id_2, laplace_nb);
	CGAL_Output_Obj(out_path, vecs, face_id_0, face_id_1, face_id_2);
}

Vector3d LaplaceMeshSmoothForOnePoint(int id, Vector3d1 &vecs, Vector3d1 &min_curvature_direction, Vector3d1 &vecs_normals, std::vector<std::vector<std::vector<int>>> &surface_vectices_to_neighbor_edges)
{
	Vector3d normal = vecs_normals[id];
	Vector3d position = vecs[id];
	Vector3d curvature_direction = min_curvature_direction[id];
	std::vector<std::vector<int>> &edges = surface_vectices_to_neighbor_edges[id];

	Vector3d product_vector = Math::GetCrossproduct(normal, curvature_direction);

	Vector3d1 inters;

	for (int i = 0; i < edges.size(); i++)
	{
		Vector3d inter;

		if (CGAL_3D_Intersection_Segment_Plane(vecs[edges[i][0]], vecs[edges[i][1]], position, product_vector, inter))
		{
			inters.push_back(inter);
		}
	}

	if (inters.size() == 2)
	{
		Vector3d surface_normal = Math::GetCrossproduct(Math::GetCrossproduct(normal, inters[1] - inters[0]), inters[1] - inters[0]);
		Vector3d inter;
		if (CGAL_3D_Intersection_Line_Plane(position, position + normal, inters[0], surface_normal, inter))
			return inter;
		else
			return position;
	}
	return 	position;
}

void CGAL_Mesh_Laplace_Smooth_by_Curvature(Vector3d1 &vecs, std::vector<int> &face_id_0, std::vector<int> &face_id_1, std::vector<int> &face_id_2, double low_curvature)
{
	std::vector<double> max_curvature;
	std::vector<double> min_curvature;
	Vector3d1 max_curvature_direction;
	Vector3d1 min_curvature_direction;

	std::vector<bool> boundary;
	CGAL_3D_Triangle_Mesh_Boundary(vecs, face_id_0, face_id_1, face_id_2, boundary);

	Vector3d1 vecs_normals;

	int last_nb=0;

	std::vector<std::vector<int>> vecs_neighbors;
	CGAL_3D_Triangle_Mesh_Vecs_Neighbors(vecs, face_id_0, face_id_1, face_id_2, vecs_neighbors);

	std::vector<std::vector<std::vector<int>>> surface_vectices_to_neighbor_edges;
	CGAL_3D_Triangle_Mesh_Vecs_Neighbor_Edges(vecs, face_id_0, face_id_1, face_id_2, surface_vectices_to_neighbor_edges);

	int stop = 0;

	double par_0 = 0.1;
	double par_1 = 0.6;
	double par_2 = 0.3;

	for (int iter = 0; iter < 500;iter++)
	//while (true)
	{
		//compute vertices curvature
		std::vector<double>().swap(max_curvature);
		std::vector<double>().swap(min_curvature);
		Vector3d1().swap(max_curvature_direction);
		Vector3d1().swap(min_curvature_direction);
		Vector3d1().swap(vecs_normals);

		CGAL_3D_Mesh_Curvature(vecs, face_id_0, face_id_1, face_id_2, max_curvature, min_curvature, max_curvature_direction, min_curvature_direction, vecs_normals);

		int nb = 0;
		double minimal_cur = 100000.0;
		for (int i = 0; i < vecs.size(); i++)
		{
			if (min_curvature[i] < low_curvature&&!boundary[i])
			{
				nb++;
				minimal_cur = std::min(minimal_cur, min_curvature[i]);
			}
		}
	

		//terminal condition

		if (nb < 5)break;

		if (abs(last_nb-nb)<2)
			stop++;
		else
			stop = 0;

		if (stop == 60)
		{
			//break;
			par_0 = 0.10;
			par_1 = 0.65;
			par_2 = 0.25;
		}

		if (stop == 100)
		{
			par_0 = 0.10;
			par_1 = 0.70;
			par_2 = 0.20;
		}

		if (stop == 300)
		{
			break;
		}

		std::cout << "Current low curvature points number: " << stop <<" "<< nb << " " << minimal_cur << std::endl;

		last_nb = nb;

		//one iteration
		Vector3d1 iteration_vecs;

		for (int i = 0; i < vecs.size(); i++)
		{
			if (boundary[i])
			{
				iteration_vecs.push_back(vecs[i]);
			}
			else
			{
				bool run = min_curvature[i] < low_curvature + 0.1;

				for (int j = 0; j < vecs_neighbors[i].size() && !run; j++)
				{
					if (min_curvature[vecs_neighbors[i][j]] < low_curvature)
					{
						run = true;
					}
					
					if (true)
					{
						int vertice_id = vecs_neighbors[i][j];

						for (int k = 0; k < vecs_neighbors[vertice_id].size() && !run; k++)
						{
							if (min_curvature[vecs_neighbors[vertice_id][k]] < low_curvature)
							{
								run = true;
							}
						}
					}
				}

				if (run)
				{
					Vector3d cur_v = vecs[i] + Math::SetVectorLength(vecs_normals[i], 0.001*min_curvature[i] / low_curvature);

					Vector3d smooth_v(0.0,0.0,0.0);

					double weight = 0.0;
					for (int j = 0; j < vecs_neighbors[i].size(); j++)
					{
						double l = CGAL_3D_Distance_Point_Point(vecs[vecs_neighbors[i][j]],vecs[i]);
						smooth_v += vecs[vecs_neighbors[i][j]]*(float)l;
						weight += l;
					}
					smooth_v = smooth_v / (float)weight;

					if (min_curvature[i]<0.0&&max_curvature[i]>0.0)
					{
						//smooth_v = vecs[i];
						iteration_vecs.push_back((float)par_0*vecs[i] + (float)par_1*cur_v + (float)par_2*smooth_v);
					}
					else
						iteration_vecs.push_back((float)par_0*vecs[i] + (float)par_1*cur_v + (float)par_2*smooth_v);

					//iteration_vecs.push_back(LaplaceMeshSmoothForOnePoint(i, vecs, min_curvature_direction, vecs_normals, surface_vectices_to_neighbor_edges));
				}
				else
				{
					iteration_vecs.push_back(vecs[i]);
				}
			}
		}

		Vector3d1().swap(vecs);
		vecs = iteration_vecs;
		Vector3d1().swap(iteration_vecs);
	}

}

void CGAL_Mesh_Laplace_Smooth(Vector3d1 &vecs, std::vector<int> &face_id_0, std::vector<int> &face_id_1, std::vector<int> &face_id_2, int laplace_nb)
{
	std::vector<bool> vertices_boundary;
	CGAL_3D_Triangle_Mesh_Boundary(vecs, face_id_0, face_id_1, face_id_2, vertices_boundary);
	std::vector<std::vector<int>> neighs;
	CGAL_3D_Triangle_Mesh_Vecs_Neighbors(vecs, face_id_0, face_id_1, face_id_2, neighs);
	
	for (int iter = 0; iter < laplace_nb; iter++)
	{
		Vector3d1 new_vecs;
		for (int i = 0; i < vecs.size(); i++)
		{
			if (!vertices_boundary[i])
			{
				Vector3d v(0.0,0.0,0.0);
				double w = 0.0;
				for (int j = 0; j < neighs[i].size(); j++)
				{
					double d = CGAL_3D_Distance_Point_Point(vecs[neighs[i][j]], vecs[i]);
					v += vecs[neighs[i][j]]*(float)(1.0/d);
					w += (1.0 / d);
				}
				v = vecs[i]*(float)0.5 + (float)0.5*v / (float)w;

				new_vecs.push_back(v);
			}
			else
			{
				new_vecs.push_back(vecs[i]);
			}
		
		}
		vecs = new_vecs;
	}

}

void CGAL_Mesh_Loop_Subdivision_One_Step(Vector3d1 &vecs, std::vector<int> &face_id_0, std::vector<int> &face_id_1, std::vector<int> &face_id_2)
{
	Vector3d1 loop_vecs = vecs;
	std::vector<int> loop_face_id_0;
	std::vector<int> loop_face_id_1;
	std::vector<int> loop_face_id_2;
	
	//edges
	std::vector<Edge> edges;

	std::vector<std::vector<int>> vecs_neighbors(vecs.size(), std::vector<int>());
	std::vector<std::vector<int>> vecs_neighbors_labels(vecs.size(), std::vector<int>());

	for (int i = 0; i < face_id_0.size(); i++)
	{
		int index_0 = face_id_0[i];
		int index_1 = face_id_1[i];
		int index_2 = face_id_2[i];
		edges.push_back(Edge(index_0, index_1));
		edges.push_back(Edge(index_1, index_0));
		edges.push_back(Edge(index_0, index_2));
		edges.push_back(Edge(index_2, index_0));
		edges.push_back(Edge(index_1, index_2));
		edges.push_back(Edge(index_2, index_1));
	}

	for (int i = 0; i < edges.size(); i++)
	{
		int source = edges[i].source;
		int end = edges[i].end;

		bool b = false;
		for (int j = 0; j < vecs_neighbors[source].size(); j++)
		{
			if (vecs_neighbors[source][j] == end)
			{
				b = true;
				break;
			}
		}
		if (!b)
		{
			vecs_neighbors[source].push_back(end);
			vecs_neighbors_labels[source].push_back(-1);
		}
	}
	std::vector<Edge>().swap(edges);

	for (int i = 0; i < vecs.size(); i++)
	{
		int  source = i;
		for (int j = 0; j < vecs_neighbors[i].size(); j++)
		{
			int end = vecs_neighbors[i][j];
			if (vecs_neighbors_labels[i][j] < 0)
			{
				edges.push_back(Edge(source,end));
				vecs_neighbors_labels[i][j] = edges.size() - 1;

				for (int k = 0; k < vecs_neighbors[end].size(); k++)
				{
					if (vecs_neighbors[end][k] == source)
					{
						vecs_neighbors_labels[end][k] = edges.size() - 1;
					}
				}
			}
		}
	}

	//loop_vecs
	Vector3d1 edge_middle_points;
	for (int i = 0; i < edges.size(); i++)
		edge_middle_points.push_back((vecs[edges[i].source] + vecs[edges[i].end])/(float)2.0);
	for (int i = 0; i< edge_middle_points.size(); i++) loop_vecs.push_back(edge_middle_points[i]);

	//loop faces
	std::vector<Edge> face_edges;
	std::vector<int> face_edges_id;

	for (int i = 0; i < face_id_0.size(); i++){
		int index_0 = face_id_0[i];
		int index_1 = face_id_1[i];
		int index_2 = face_id_2[i];
		face_edges.push_back(Edge(index_0, index_1));
		face_edges.push_back(Edge(index_1, index_2));
		face_edges.push_back(Edge(index_2, index_0));
	}
	for (int i = 0; i < face_edges.size(); i++)
	{
		int source = face_edges[i].source;
		int end = face_edges[i].end;
	
		for (int j = 0; j < vecs_neighbors[source].size(); j++)
		{
			if (vecs_neighbors[source][j] == end)
			{
				face_edges_id.push_back(vecs_neighbors_labels[source][j]);
				break;
			}
		}
	}

	//     0
	//    2 0
	//  2  1  1
	for (int i = 0; i < face_edges.size(); i=i+3)
	{
		int index_0 = face_edges[i].source;
		int index_1 = face_edges[i+1].source;
		int index_2 = face_edges[i+2].source;

		int edge_id_0 = face_edges_id[i]+vecs.size();
		int edge_id_1 = face_edges_id[i + 1] + vecs.size();
		int edge_id_2 = face_edges_id[i + 2] + vecs.size();

		loop_face_id_0.push_back(index_0);
		loop_face_id_1.push_back(edge_id_0);
		loop_face_id_2.push_back(edge_id_2);

		loop_face_id_0.push_back(edge_id_0);
		loop_face_id_1.push_back(index_1);
		loop_face_id_2.push_back(edge_id_1);

		loop_face_id_0.push_back(edge_id_2);
		loop_face_id_1.push_back(edge_id_1);
		loop_face_id_2.push_back(index_2);

		loop_face_id_0.push_back(edge_id_2);
		loop_face_id_1.push_back(edge_id_0);
		loop_face_id_2.push_back(edge_id_1);
	}

	//release
	Vector3d1().swap(vecs);
	std::vector<int>().swap(face_id_0);
	std::vector<int>().swap(face_id_1);
	std::vector<int>().swap(face_id_2);

	vecs = loop_vecs;
	face_id_0 = loop_face_id_0;
	face_id_1 = loop_face_id_1;
	face_id_2 = loop_face_id_2;

	Vector3d1().swap(loop_vecs);
	std::vector<int>().swap(loop_face_id_0);
	std::vector<int>().swap(loop_face_id_1);
	std::vector<int>().swap(loop_face_id_2);

	Vector3d1().swap(edge_middle_points);
	std::vector<Edge>().swap(face_edges);
	std::vector<int>().swap(face_edges_id);
	std::vector<Edge>().swap(edges);
}

/***************************************************************************************************/


//slice the input mesh
/***************************************************************************************************/

void CGAL_Slicer_Mesh(std::string path, double normal_x, double normal_y, double normal_z, std::vector<double>  plane_d,
	std::vector<std::vector<std::vector<double>>>& slice_x, std::vector<std::vector<std::vector<double>>>& slice_y, std::vector<std::vector<std::vector<double>>>& slice_z)
{

	std::ifstream input(path.c_str());
	Mesh mesh;
	if (!input || !(input >> mesh) || mesh.is_empty()) {
		std::cerr << "Not a valid off file." << std::endl;
		return;
	}
	// Slicer constructor from the mesh
	CGAL::Polygon_mesh_slicer<Mesh, K> slicer(mesh);
	Polylines polylines;

	// Use the Slicer constructor from a pre-built AABB_tree
	AABB_tree tree(edges(mesh).first, edges(mesh).second, mesh);


	CGAL::Polygon_mesh_slicer<Mesh, K> slicer_aabb(mesh, tree);

	for (int i = 0; i < plane_d.size(); i++)
	{
		slicer_aabb(K::Plane_3(normal_x, normal_y, normal_z, -plane_d[i]), std::back_inserter(polylines));

		std::vector<std::vector<double>> xs;
		std::vector<std::vector<double>> ys;
		std::vector<std::vector<double>> zs;

		Polylines::iterator iter;
		for (iter = polylines.begin(); iter != polylines.end(); iter++)
		{
			std::vector<double> x;
			std::vector<double> y;
			std::vector<double> z;
			Polyline_type::iterator p_iter;
			for (p_iter = iter->begin(); p_iter != iter->end(); p_iter++)
			{
				x.push_back(p_iter->x());
				y.push_back(p_iter->y());
				z.push_back(p_iter->z());
			}
			x.pop_back();
			y.pop_back();
			z.pop_back();

			xs.push_back(x);
			ys.push_back(y);
			zs.push_back(z);

			std::vector<double>().swap(x);
			std::vector<double>().swap(y);
			std::vector<double>().swap(z);
		}

		polylines.clear();

		slice_x.push_back(xs);
		slice_y.push_back(ys);
		slice_z.push_back(zs);

		std::vector<std::vector<double>>().swap(xs);
		std::vector<std::vector<double>>().swap(ys);
		std::vector<std::vector<double>>().swap(zs);
	}
}

void CGAL_Slicer_Mesh(std::string path, Vector3d plane_normal, std::vector<double> plane_d,
	Vector3d3& offsetses, Vector3d2& offsets)
{
	std::ifstream input(path.c_str());
	Mesh mesh;
	if (!input || !(input >> mesh) || mesh.is_empty()) {
		std::cerr << "Not a valid off file." << std::endl;
		return;
	}
	// Slicer constructor from the mesh
	CGAL::Polygon_mesh_slicer<Mesh, K> slicer(mesh);
	Polylines polylines;

	// Use the Slicer constructor from a pre-built AABB_treen
	AABB_tree tree(edges(mesh).first, edges(mesh).second, mesh);
	
	CGAL::Polygon_mesh_slicer<Mesh, K> slicer_aabb(mesh, tree);

	for (int i = 0; i < plane_d.size(); i++)
	{
		std::cout << i << "/" << plane_d.size() << std::endl;

		slicer_aabb(K::Plane_3(plane_normal[0], plane_normal[1], plane_normal[2], -plane_d[i]), std::back_inserter(polylines));

		std::vector<std::vector<double>> xs;
		std::vector<std::vector<double>> ys;
		std::vector<std::vector<double>> zs;

		Vector3d2 circles;

		Polylines::iterator iter;
		for (iter = polylines.begin(); iter != polylines.end(); iter++)
		{
			Vector3d1 one_circle;

			Polyline_type::iterator p_iter;
			for (p_iter = iter->begin(); p_iter != iter->end(); p_iter++)
			{
				one_circle.push_back(Vector3d(p_iter->x(), p_iter->y(), p_iter->z()));
			}
			circles.push_back(one_circle);
			offsets.push_back(one_circle);
		}
		polylines.clear();
		offsetses.push_back(circles);
	}
}

/***************************************************************************************************/

//shortest geodesic path
/***************************************************************************************************/
void CGAL_Shortest_Geodesic_Path(std::string path, std::vector<double> &x, std::vector<double> &y, std::vector<double> &z)
{
	// read input polyhedron
	Polyhedron_3 polyhedron;
	Construct_Polyhedron(polyhedron,path);
	// pick up a random face
	const size_t randSeed = 7915421;
	CGAL::Random rand(randSeed);
	const int target_face_index = rand.get_int(0, num_faces(polyhedron));
	face_iterator face_it = faces(polyhedron).first;
	std::advance(face_it, target_face_index);
	// ... and define a barycentric coordinate inside the face
	Traits::Barycentric_coordinate face_location = { { 0.25, 0.5, 0.25 } };
	// construct a shortest path query object and add a source point
	Surface_mesh_shortest_path shortest_paths(polyhedron);
	shortest_paths.add_source_point(*face_it, face_location);

	vertex_iterator vit = polyhedron.vertices_begin();
	std::vector<Traits::Point_3> points;
	shortest_paths.shortest_path_points_to_source_points(*vit, std::back_inserter(points));

	for (int i = 0; i < points.size(); i++)
	{
		x.push_back(points[i][0]);
		y.push_back(points[i][1]);
		z.push_back(points[i][2]);
	}

}

/***************************************************************************************************/

Poly_point_3 Minus(Poly_point_3 a, Poly_point_3 b)
{
	return Poly_point_3(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}
float Dot(Poly_point_3 a, Poly_point_3 b)
{
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

// Compute barycentric coordinates (u, v, w) for
// point p with respect to triangle (a, b, c)
void Barycentric(Poly_point_3 p, Poly_point_3 a, Poly_point_3 b, Poly_point_3 c, double &u, double &v, double &w)
{
	Poly_point_3 v0 = Minus(b, a), v1 = Minus(c, a), v2 = Minus(p,a);

	double d00 = Dot(v0, v0);
	double d01 = Dot(v0, v1);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
	double d11 = Dot(v1, v1);
	double d20 = Dot(v2, v0);
	double d21 = Dot(v2, v1);
	double denom = d00 * d11 - d01 * d01;
	v = (d11 * d20 - d01 * d21) / denom;
	w = (d00 * d21 - d01 * d20) / denom;
	u = 1.0f - v - w;
}

void CGAL_Barycentric(Vector3d p, Vector3d a, Vector3d b, Vector3d c, double &u, double &v, double &w)
{
	Barycentric(Poly_point_3(p[0], p[1], p[2]), Poly_point_3(a[0], a[1], a[2]), Poly_point_3(b[0], b[1], b[2]), Poly_point_3(c[0], c[1], c[2]), u, v, w);
}


void RelatedFaceAndBarycentric(Polyhedron_3 &polyhedron, Tree &tree, 
	double source_x, double source_y, double source_z, double &u, double &v, double &w, Poly_point_3 &nearest_point, face_iterator &face_it)
{
	Poly_point_3 query(source_x, source_y, source_z);
	Point_and_primitive_id pp = tree.closest_point_and_primitive(query);
	nearest_point = pp.first;
	face_it = pp.second;
	
	Poly_point_3 p0 = pp.second->halfedge()->next()->next()->vertex()->point();
	Poly_point_3 p1 = pp.second->halfedge()->vertex()->point();
	Poly_point_3 p2 = pp.second->halfedge()->next()->vertex()->point();
	
	Barycentric(pp.first, p0, p1, p2, u, v, w);
}


void CGAL_Shortest_Geodesic_Path(Polyhedron_3 &polyhedron, Tree &tree,
	double source_x, double source_y, double source_z, double target_x, double target_y, double target_z, 
	std::vector<double> &x, std::vector<double> &y, std::vector<double> &z)
{
	////////////////////////////////////////
	face_iterator source_face, target_face;
	double source_x_w, source_y_w, source_z_w;
	double target_x_w, target_y_w, target_z_w;
	Poly_point_3 source_nearest_point, target_nearest_point;

	RelatedFaceAndBarycentric(polyhedron, tree, source_x, source_y, source_z, source_x_w, source_y_w, source_z_w, source_nearest_point, source_face);
	RelatedFaceAndBarycentric(polyhedron, tree, target_x, target_y, target_z, target_x_w, target_y_w, target_z_w, target_nearest_point, target_face);

	Traits::Barycentric_coordinate source_face_location = { { source_x_w, source_y_w, source_z_w } };
	Traits::Barycentric_coordinate target_face_location = { { target_x_w, target_y_w, target_z_w } };
	//////////////////////////////////////////////////////////////

	Surface_mesh_shortest_path shortest_paths(polyhedron);
	shortest_paths.add_source_point(*source_face, source_face_location);

	std::vector<Traits::Point_3> points;
	shortest_paths.shortest_path_points_to_source_points(*target_face, target_face_location, std::back_inserter(points));

	for (int i = points.size() - 1; i >= 0; i--)
	{
		x.push_back(points[i][0]);
		y.push_back(points[i][1]);
		z.push_back(points[i][2]);
	}
}

void CGAL_Shortest_Geodesic_Path(std::string path, Vector3d source, Vector3d target, Vector3d1 &output)
{
	Polyhedron_3 polyhedron;
	Construct_Polyhedron(polyhedron, path);

	Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);
	tree.accelerate_distance_queries();

	//////////////////////////////////////////////////////////////
	face_iterator source_face, target_face;
	double source_x_w, source_y_w, source_z_w;
	double target_x_w, target_y_w, target_z_w;
	Poly_point_3 source_nearest_point, target_nearest_point;

	RelatedFaceAndBarycentric(polyhedron, tree, source[0], source[1], source[2], source_x_w, source_y_w, source_z_w, source_nearest_point, source_face);
	RelatedFaceAndBarycentric(polyhedron, tree, target[0], target[1], target[2], target_x_w, target_y_w, target_z_w, target_nearest_point, target_face);

	Traits::Barycentric_coordinate source_face_location = { { source_x_w, source_y_w, source_z_w } };
	Traits::Barycentric_coordinate target_face_location = { { target_x_w, target_y_w, target_z_w } };
	//////////////////////////////////////////////////////////////

	Surface_mesh_shortest_path shortest_paths(polyhedron);
	shortest_paths.add_source_point(*source_face, source_face_location);

	std::vector<Traits::Point_3> points;
	shortest_paths.shortest_path_points_to_source_points(*target_face, target_face_location, std::back_inserter(points));

	for (int i = points.size() - 1; i >= 0; i--)
	{
		output.push_back(Vector3d(points[i][0], points[i][1], points[i][2]));
	}
}

void CGAL_Shortest_Geodesic_Path(std::string path, std::vector<double> sources_x, std::vector<double> sources_y, std::vector<double> sources_z, 
	std::vector<double> targets_x, std::vector<double> targets_y, std::vector<double> targets_z,
	std::vector<std::vector<double>> &xs, std::vector<std::vector<double>> &ys, std::vector<std::vector<double>> &zs)
{
	Polyhedron_3 polyhedron;

	Construct_Polyhedron(polyhedron, path);

	Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);
	tree.accelerate_distance_queries();

	std::cout << "Start to compute the geodesic path..." << std::endl;

	for (int i = 0; i < sources_x.size(); i++)
	{
		std::cout << "Path: " << i << std::endl;

		double source_x = sources_x[i];
		double source_y = sources_y[i];
		double source_z = sources_z[i];

		double target_x = targets_x[i];
		double target_y = targets_y[i];
		double target_z = targets_z[i];

		std::vector<double> x, y, z;
		CGAL_Shortest_Geodesic_Path(polyhedron, tree, source_x, source_y, source_z, target_x, target_y, target_z, x, y, z);

		xs.push_back(x);
		ys.push_back(y);
		zs.push_back(z);

		std::vector<double>().swap(x);
		std::vector<double>().swap(y);
		std::vector<double>().swap(z);
	}
}



/***************************************************************************************************/


//shortest geodesic distance
/***************************************************************************************************/
double CGAL_Geodesic_Distance(std::string path, Vector3d source, Vector3d target)
{
	std::cout << "one time geodesic computing.." << std::endl;

	Polyhedron_3 polyhedron;
	Construct_Polyhedron(polyhedron, path);

	Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);
	tree.accelerate_distance_queries();

	//////////////////////////////////////////////////////////////
	face_iterator source_face, target_face;
	double source_x_w, source_y_w, source_z_w;
	double target_x_w, target_y_w, target_z_w;
	Poly_point_3 source_nearest_point, target_nearest_point;

	RelatedFaceAndBarycentric(polyhedron, tree, source[0], source[1], source[2], source_x_w, source_y_w, source_z_w, source_nearest_point, source_face);
	RelatedFaceAndBarycentric(polyhedron, tree, target[0], target[1], target[2], target_x_w, target_y_w, target_z_w, target_nearest_point, target_face);

	Traits::Barycentric_coordinate source_face_location = { { source_x_w, source_y_w, source_z_w } };
	Traits::Barycentric_coordinate target_face_location = { { target_x_w, target_y_w, target_z_w } };
	//////////////////////////////////////////////////////////////

	Surface_mesh_shortest_path shortest_paths(polyhedron);
	shortest_paths.add_source_point(*source_face, source_face_location);

	return shortest_paths.shortest_distance_to_source_points(*target_face, target_face_location).first;
}


/***************************************************************************************************/

//normal computation
/***************************************************************************************************/


//Vector_3 RelatedFaceNormal(Polyhedron_3 &polyhedron, Tree &tree, double source_x, double source_y, double source_z)
//{
//	Poly_point_3 query(source_x, source_y, source_z);
//	Point_and_primitive_id pp = tree.closest_point_and_primitive(query);
//
//	Poly_point_3 p0 = pp.second->halfedge()->next()->next()->vertex()->point();
//	Poly_point_3 p1 = pp.second->halfedge()->vertex()->point();
//	Poly_point_3 p2 = pp.second->halfedge()->next()->vertex()->point();
//
//	Vector_3   n = CGAL::cross_product(p2 - p1, p0 - p1);
//
//	return n / std::sqrt(n*n);
//}

Vector3d RelatedFaceNormal(Polyhedron_3 &polyhedron, Tree &tree,
	Vector3d1 &normals,Vector3d source)
{
	Poly_point_3 query(source[0], source[1], source[2]);
	Point_and_primitive_id pp = tree.closest_point_and_primitive(query);

	Poly_point_3 p0 = pp.second->halfedge()->next()->next()->vertex()->point();
	Poly_point_3 p1 = pp.second->halfedge()->vertex()->point();
	Poly_point_3 p2 = pp.second->halfedge()->next()->vertex()->point();

	int point_id_0 = pp.second->halfedge()->next()->next()->vertex()->id();
	int point_id_1 = pp.second->halfedge()->vertex()->id();
	int point_id_2 = pp.second->halfedge()->next()->vertex()->id();

	double u, v, w;
	Barycentric(query, p0, p1, p2, u, v, w);

	return (float)u*normals[point_id_0] + (float)v*normals[point_id_1] + (float)w*normals[point_id_2];
}

Vector3d CGAL_Face_Normal(Vector3d source, Vector3d tri_0, Vector3d tri_1, Vector3d tri_2, 
	                                        Vector3d normal_0, Vector3d normal_1, Vector3d normal_2)
{
	double u, v, w;
	CGAL_Barycentric(source, tri_0, tri_1, tri_2, u, v, w);
	return (float)u*normal_0 + (float)v*normal_1 + (float)w*normal_2;
}

void CGAL_Normal_Mesh(std::string path, Vector3d2 &mesh_pointses, Vector3d2 &mesh_normalses)
{
	std::cout << "CGAL_Normal_Mesh.." << std::endl;
	if (path.substr(path.size() - 3, path.size()) == "off")
	{
		Polyhedron_3 polyhedron;
		Construct_Polyhedron(polyhedron, path);

		Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);
		tree.accelerate_distance_queries();

		//get mesh vertices and surface
		std::vector<double> mesh_xs, mesh_ys, mesh_zs;
		std::vector<int> mesh_face_id_0, mesh_face_id_1, mesh_face_id_2;
		for (Polyhedron_3::Vertex_iterator iter = polyhedron.vertices_begin();
			iter != polyhedron.vertices_end(); iter++)
		{
			Poly_point_3 p = iter->point();
			mesh_xs.push_back(p[0]);
			mesh_ys.push_back(p[1]);
			mesh_zs.push_back(p[2]);
		}
		for (Polyhedron_3::Face_iterator iter = polyhedron.facets_begin(); iter != polyhedron.facets_end(); iter++)
		{
			mesh_face_id_0.push_back(iter->halfedge()->next()->next()->vertex()->id());
			mesh_face_id_1.push_back(iter->halfedge()->vertex()->id());
			mesh_face_id_2.push_back(iter->halfedge()->next()->vertex()->id());
		}

		//compute surface normals
		int verticeSize = mesh_xs.size();
		int faceindiceSize = mesh_face_id_0.size() * 3;

		Wm5::Vector3<double> *points = new Wm5::Vector3<double>[verticeSize];
		int *indices = new int[faceindiceSize];

		for (int i = 0; i<verticeSize; i++)
		{
			points[i].X() = mesh_xs[i];
			points[i].Y() = mesh_ys[i];
			points[i].Z() = mesh_zs[i];
		}

		for (int i = 0; i<mesh_face_id_0.size(); i++)
		{
			indices[3 * i] = mesh_face_id_0[i];
			indices[3 * i + 1] = mesh_face_id_1[i];
			indices[3 * i + 2] = mesh_face_id_2[i];
		}
		Wm5::MeshCurvature<double> meshCurvature(verticeSize, points, faceindiceSize / 3, indices);

		Vector3d1 normals;
		for (int i = 0; i < verticeSize; i++)
		{
			Wm5::Vector3<double> normal = meshCurvature.GetNormals()[i];
			normals.push_back(Vector3d(normal[0], normal[1], normal[2]));
		}

		for (int i = 0; i < mesh_pointses.size(); i++)
		{
			Vector3d1 mesh_normals;
			for (int j = 0; j < mesh_pointses[i].size(); j++)
			{
				Vector3d n = RelatedFaceNormal(polyhedron, tree, normals, mesh_pointses[i][j]);
				mesh_normals.push_back(n);
			}
			mesh_normalses.push_back(mesh_normals);
			Vector3d1().swap(mesh_normals);
		}
	}

	if (path.substr(path.size() - 3, path.size()) == "obj")
	{
		Vector3d1 vecs;
		std::vector<int> face_id_0;
		std::vector<int> face_id_1;
		std::vector<int> face_id_2;
		CGAL_3D_Read_Triangle_Mesh(path, vecs, face_id_0, face_id_1, face_id_2);

		Polyhedron_3 polyhedron;

		Construct_Polyhedron(polyhedron, vecs, face_id_0, face_id_1, face_id_2);

		Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);
		tree.accelerate_distance_queries();

		//compute normals
		/**********************************************************/
		Vector3d1 normals(vecs.size(), Vector3d(0.0, 0.0, 0.0));
		std::vector<double> areas(vecs.size(), 0.0);
		for (int i = 0; i < face_id_0.size(); i++)
		{
			double area = CGAL_3D_One_Triangle_Area(vecs[face_id_0[i]], vecs[face_id_1[i]], vecs[face_id_2[i]]);
			Poly_point_3 p0(vecs[face_id_0[i]][0], vecs[face_id_0[i]][1], vecs[face_id_0[i]][2]);
			Poly_point_3 p1(vecs[face_id_1[i]][0], vecs[face_id_1[i]][1], vecs[face_id_1[i]][2]);
			Poly_point_3 p2(vecs[face_id_2[i]][0], vecs[face_id_2[i]][1], vecs[face_id_2[i]][2]);
			Vector_3   n = CGAL::cross_product(p2 - p1, p0 - p1);

			Vector3d n0(n[0], n[1], n[2]);

			normals[face_id_0[i]] += (float)area*n0;
			normals[face_id_1[i]] += (float)area*n0;
			normals[face_id_2[i]] += (float)area*n0;
			areas[face_id_0[i]] += area;
			areas[face_id_1[i]] += area;
			areas[face_id_2[i]] += area;
		}

		for (int i = 0; i < vecs.size(); i++)
		{
			normals[i] = normals[i] / (float)areas[i];
		}
		/**********************************************************/

		for (int i = 0; i < mesh_pointses.size(); i++)
		{
			Vector3d1 mesh_normals;
			for (int j = 0; j < mesh_pointses[i].size(); j++)
			{
				Vector3d n = RelatedFaceNormal(polyhedron, tree, normals, mesh_pointses[i][j]);
				mesh_normals.push_back(n);
			}
			mesh_normalses.push_back(mesh_normals);
			Vector3d1().swap(mesh_normals);
		}

	}
}

void CGAL_Normal_Mesh(std::string path, Vector3d1 &mesh_points, Vector3d1 &mesh_normals)
{
	std::cout << "CGAL_Normal_Mesh.." << std::endl;
	if (path.substr(path.size() - 3, path.size()) == "off")
	{
		Polyhedron_3 polyhedron;
		Construct_Polyhedron(polyhedron, path);

		Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);
		tree.accelerate_distance_queries();

		//get mesh vertices and surface
		std::vector<double> mesh_xs, mesh_ys, mesh_zs;
		std::vector<int> mesh_face_id_0, mesh_face_id_1, mesh_face_id_2;
		for (Polyhedron_3::Vertex_iterator iter = polyhedron.vertices_begin();
			iter != polyhedron.vertices_end(); iter++)
		{
			Poly_point_3 p = iter->point();
			mesh_xs.push_back(p[0]);
			mesh_ys.push_back(p[1]);
			mesh_zs.push_back(p[2]);
		}
		for (Polyhedron_3::Face_iterator iter = polyhedron.facets_begin(); iter != polyhedron.facets_end(); iter++)
		{
			mesh_face_id_0.push_back(iter->halfedge()->next()->next()->vertex()->id());
			mesh_face_id_1.push_back(iter->halfedge()->vertex()->id());
			mesh_face_id_2.push_back(iter->halfedge()->next()->vertex()->id());
		}

		//compute surface normals
		int verticeSize = mesh_xs.size();
		int faceindiceSize = mesh_face_id_0.size() * 3;

		Wm5::Vector3<double> *points = new Wm5::Vector3<double>[verticeSize];
		int *indices = new int[faceindiceSize];

		for (int i = 0; i<verticeSize; i++)
		{
			points[i].X() = mesh_xs[i];
			points[i].Y() = mesh_ys[i];
			points[i].Z() = mesh_zs[i];
		}

		for (int i = 0; i<mesh_face_id_0.size(); i++)
		{
			indices[3 * i] = mesh_face_id_0[i];
			indices[3 * i + 1] = mesh_face_id_1[i];
			indices[3 * i + 2] = mesh_face_id_2[i];
		}
		Wm5::MeshCurvature<double> meshCurvature(verticeSize, points, faceindiceSize / 3, indices);

		Vector3d1 normals;
		for (int i = 0; i < verticeSize; i++)
		{
			Wm5::Vector3<double> normal = meshCurvature.GetNormals()[i];
			normals.push_back(Vector3d(normal[0], normal[1], normal[2]));
		}

		for (int i = 0; i < mesh_points.size(); i++)
		{
			Vector3d n = RelatedFaceNormal(polyhedron, tree, normals, mesh_points[i]);
			mesh_normals.push_back(n);
		}
	}

	if (path.substr(path.size() - 3, path.size()) == "obj")
	{
		Vector3d1 vecs;
		std::vector<int> face_id_0;
		std::vector<int> face_id_1;
		std::vector<int> face_id_2;
		CGAL_3D_Read_Triangle_Mesh(path, vecs, face_id_0, face_id_1, face_id_2);

		Polyhedron_3 polyhedron;

		Construct_Polyhedron(polyhedron, vecs, face_id_0, face_id_1, face_id_2);

		Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);
		tree.accelerate_distance_queries();

		//compute normals
		/**********************************************************/
		Vector3d1 normals(vecs.size(),Vector3d(0.0,0.0,0.0));
		std::vector<double> areas(vecs.size(),0.0);
		for (int i = 0; i < face_id_0.size(); i++)
		{
			double area = CGAL_3D_One_Triangle_Area(vecs[face_id_0[i]], vecs[face_id_1[i]], vecs[face_id_2[i]]);
			Poly_point_3 p0(vecs[face_id_0[i]][0], vecs[face_id_0[i]][1], vecs[face_id_0[i]][2]);
			Poly_point_3 p1(vecs[face_id_1[i]][0], vecs[face_id_1[i]][1], vecs[face_id_1[i]][2]);
			Poly_point_3 p2(vecs[face_id_2[i]][0], vecs[face_id_2[i]][1], vecs[face_id_2[i]][2]);
			Vector_3   n = CGAL::cross_product(p2 - p1, p0 - p1);

			Vector3d n0(n[0],n[1],n[2]);

			normals[face_id_0[i]] += (float)area*n0;
			normals[face_id_1[i]] += (float)area*n0;
			normals[face_id_2[i]] += (float)area*n0;
			areas[face_id_0[i]] += area;
			areas[face_id_1[i]] += area;
			areas[face_id_2[i]] += area;
		}
		
		for (int i = 0; i < vecs.size(); i++)
		{
			normals[i] = normals[i] / (float)areas[i];
		}
		/**********************************************************/

		for (int i = 0; i < mesh_points.size(); i++)
		{
			Vector3d n = RelatedFaceNormal(polyhedron, tree, normals, mesh_points[i]);
			mesh_normals.push_back(n);
		}
	}
}

void CGAL_3D_Read_Triangle_Mesh(std::string path, Vector3d1& vecs,
	std::vector<int>& face_id_0, std::vector<int>& face_id_1, std::vector<int>& face_id_2)
{
	if (path.substr(path.size() - 3, path.size()) == "off")
	{
		Polyhedron_3 polyhedron;
		Construct_Polyhedron(polyhedron, path);

		for (Polyhedron_3::Vertex_iterator iter = polyhedron.vertices_begin();
			iter != polyhedron.vertices_end(); iter++)
		{
			Poly_point_3 p = iter->point();
			vecs.push_back(Vector3d(p[0], p[1], p[2]));
		}

		for (Polyhedron_3::Face_iterator iter = polyhedron.facets_begin(); iter != polyhedron.facets_end(); iter++)
		{
			//Poly_point_3 p0 = iter->halfedge()->next()->next()->vertex()->point();
			//Poly_point_3 p1 = iter->halfedge()->vertex()->point();
			//Poly_point_3 p2 = iter->halfedge()->next()->vertex()->point();
			face_id_0.push_back(iter->halfedge()->next()->next()->vertex()->id());
			face_id_1.push_back(iter->halfedge()->vertex()->id());
			face_id_2.push_back(iter->halfedge()->next()->vertex()->id());
		}
	}

	if (path.substr(path.size() - 3, path.size()) == "obj")
	{
		std::vector<double> coords;
		std::vector<int>    tris;
		CGAL_Load_Obj(path.c_str(), coords, tris);
		if (coords.size() == 0)
			return;

		std::cout << "Size: " << coords.size() / 3 << " " << tris.size() / 3 << std::endl;

		for (int i = 0; i<(int)coords.size(); i += 3){
			vecs.push_back(Vector3d(coords[i + 0], coords[i + 1], coords[i + 2]));
		}

		for (int i = 0; i<(int)tris.size(); i += 3){
			face_id_0.push_back(tris[i + 0]);
			face_id_1.push_back(tris[i + 1]);
			face_id_2.push_back(tris[i + 2]);
		}
		/*********************************************************************************/
	}
}

void CGAL_3D_Read_Triangle_Mesh(std::string path, Vector3d1& vecs, std::vector<std::vector<int>>& face_ids)
{
	std::vector<int> face_id_0;
	std::vector<int> face_id_1;
	std::vector<int> face_id_2;
	CGAL_3D_Read_Triangle_Mesh(path, vecs, face_id_0, face_id_1, face_id_2);
	for (int i = 0; i < face_id_0.size(); i++)
	{
		std::vector<int> face;
		face.push_back(face_id_0[i]);
		face.push_back(face_id_1[i]);
		face.push_back(face_id_2[i]);
		face_ids.push_back(face);
	}
}

double CGAL_3D_One_Triangle_Area(Vector3d v0, Vector3d v1, Vector3d v2)
{
	double a = CGAL_3D_Distance_Point_Point(v0[0], v0[1], v0[2], v1[0], v1[1], v1[2]);
	double b = CGAL_3D_Distance_Point_Point(v2[0], v2[1], v2[2], v1[0], v1[1], v1[2]);
	double c = CGAL_3D_Distance_Point_Point(v0[0], v0[1], v0[2], v2[0], v2[1], v2[2]);
	double p = (a + b + c) / 2.0;
	return sqrt(p*(p - a)*(p - b)*(p - c));
}

//compute the triangle mesh area
double CGAL_3D_Triangle_Mesh_Area(Vector3d1 &vecs, std::vector<int> &face_id_0, std::vector<int> &face_id_1, std::vector<int> &face_id_2)
{
	double area = 0.0;

	for (int i = 0; i < face_id_0.size(); i++){
		int index_0 = face_id_0[i];
		int index_1 = face_id_1[i];
		int index_2 = face_id_2[i];

		Vector3d v0 = vecs[index_0];
		Vector3d v1 = vecs[index_1];
		Vector3d v2 = vecs[index_2];
		area += CGAL_3D_One_Triangle_Area(v0,v1,v2);
	}

	return area;
}

void CGAL_3D_Triangle_Mesh_Vecs_Neighbors(Vector3d1 &vecs, std::vector<int> &face_id_0, std::vector<int> &face_id_1, std::vector<int> &face_id_2,std::vector<std::vector<int>> &neighs)
{
	for (int i = 0; i < vecs.size(); i++)
		neighs.push_back(std::vector<int>());

	for (int i = 0; i < face_id_0.size(); i++)
	{
		int id_0 = face_id_0[i];
		int id_1 = face_id_1[i];
		int id_2 = face_id_2[i];

		neighs[id_0].push_back(id_1);
		neighs[id_0].push_back(id_2);

		neighs[id_1].push_back(id_0);
		neighs[id_1].push_back(id_2);

		neighs[id_2].push_back(id_0);
		neighs[id_2].push_back(id_1);
	}

	for (int i = 0; i < neighs.size(); i++)
	{
		sort(neighs[i].begin(), neighs[i].end());
		neighs[i].erase(unique(neighs[i].begin(), neighs[i].end()), neighs[i].end());
	}
}

void CGAL_3D_Triangle_Mesh_Vecs_Faces(Vector3d1 &vecs, std::vector<int> &face_id_0, std::vector<int> &face_id_1, std::vector<int> &face_id_2, 
	std::vector<std::vector<int>> &surface_vectices_to_face)
{
	surface_vectices_to_face = std::vector<std::vector<int>>(vecs.size(), std::vector<int>());

	std::vector<std::vector<int>> sets(vecs.size(), std::vector<int>());
	for (int i = 0; i < face_id_0.size(); i++)
	{
		//surface_vectices_to_face
		sets[face_id_0[i]].emplace_back(i);
		sets[face_id_1[i]].emplace_back(i);
		sets[face_id_2[i]].emplace_back(i);
	}

	for (int i = 0; i < vecs.size(); i++)
	{
		set<int>s(sets[i].begin(), sets[i].end());
		vector<int> vec;
		vec.assign(s.begin(), s.end());
		surface_vectices_to_face[i] = vec;
	}

}

void CGAL_3D_Triangle_Mesh_Vecs_Neighbor_Edges(Vector3d1 &vecs, std::vector<int> &face_id_0, std::vector<int> &face_id_1, std::vector<int> &face_id_2, 
	std::vector<std::vector<std::vector<int>>> &surface_vectices_to_neighbor_edges)
{
	std::vector<std::vector<int>> surface_vectices_to_face;
	CGAL_3D_Triangle_Mesh_Vecs_Faces(vecs, face_id_0, face_id_1, face_id_2, surface_vectices_to_face);

	for (int i = 0; i < vecs.size(); i++)
	{
		int vertice_id = i;

		std::vector<std::vector<int>> edges;

		for (int j = 0; j < surface_vectices_to_face[i].size(); j++)
		{
			int surface_id = surface_vectices_to_face[i][j];

			std::vector<int> face;
			face.push_back(face_id_0[surface_id]);
			face.push_back(face_id_1[surface_id]);
			face.push_back(face_id_2[surface_id]);

			for (int k = 0; k < face.size(); k++)
			{
				if (face[k] == vertice_id)
				{
					int vertice_id_0 = face[(k + 1) % 3];
					int vertice_id_1 = face[(k + 2) % 3];

					std::vector<int> edge;
					edge.push_back(vertice_id_0);
					edge.push_back(vertice_id_1);
					edges.push_back(edge);
					break;
				}
			}
		}
		surface_vectices_to_neighbor_edges.push_back(edges);
	}
}


//search for the mesh boundary
void CGAL_3D_Triangle_Mesh_Boundary(Vector3d1 &vecs,
	std::vector<int> &face_id_0, std::vector<int> &face_id_1, std::vector<int> &face_id_2, 
	std::vector<bool> &bools)
{
	std::vector<bool>().swap(bools);

	std::vector<std::vector<int>> vecs_neigbor(vecs.size(),std::vector<int>());
	std::vector<std::vector<int>> vecs_neigbor_lable(vecs.size(), std::vector<int>());
	std::vector<int> edges;
	for (int i = 0; i < face_id_0.size(); i++){
		int index_0 = face_id_0[i];
		int index_1 = face_id_1[i];
		int index_2 = face_id_2[i];
		edges.push_back(index_0);
		edges.push_back(index_1);
		edges.push_back(index_1);
		edges.push_back(index_2);
		edges.push_back(index_2);
		edges.push_back(index_0);

		edges.push_back(index_1);
		edges.push_back(index_0);
		edges.push_back(index_2);
		edges.push_back(index_1);
		edges.push_back(index_0);
		edges.push_back(index_2);
	}

	for (int i = 0; i < edges.size(); i = i + 2){
		int index_0 = edges[i];
		int index_1 = edges[i+1];

		int search_0 = -1;
		for (int j = 0; j < vecs_neigbor[index_0].size() && search_0<0; j++){
			if (vecs_neigbor[index_0][j] == index_1){
				search_0 = j;
				vecs_neigbor_lable[index_0][j]++;
			}
		}

		if (search_0 < 0){
			vecs_neigbor[index_0].push_back(index_1);
			vecs_neigbor_lable[index_0].push_back(1);
		}
	}

	for (int i = 0; i < vecs.size(); i++){
		bool b = false;
		for (int j = 0; j < vecs_neigbor_lable[i].size() &!b; j++){
			if (vecs_neigbor_lable[i][j] == 1){
				b = true;
			}
		}
		bools.push_back(b);
	}

	std::vector<std::vector<int>>().swap(vecs_neigbor);
	std::vector<std::vector<int>>().swap(vecs_neigbor_lable);
	std::vector<int>().swap(edges);
}

void CGAL_3D_Triangle_Mesh_Boundary(std::string path, std::vector<bool> &bools)
{
	Vector3d1 vecs;
	std::vector<int> face_id_0;
	std::vector<int> face_id_1;
	std::vector<int> face_id_2;
	CGAL_3D_Read_Triangle_Mesh(path, vecs, face_id_0, face_id_1, face_id_2);
	CGAL_3D_Triangle_Mesh_Boundary(vecs, face_id_0, face_id_1, face_id_2, bools);
}

void CGAL_3D_Triangle_Mesh_Boundary(Vector3d1& vecs, std::vector<std::vector<int>>& face_ids, Vector3d2 &boundaries)
{
	std::vector<int> face_id_0;
	std::vector<int> face_id_1;
	std::vector<int> face_id_2;

	for (auto &face : face_ids)
	{
		face_id_0.emplace_back(face[0]);
		face_id_1.emplace_back(face[1]);
		face_id_2.emplace_back(face[2]);
	}

	CGAL_3D_Triangle_Mesh_Boundary(vecs, face_id_0, face_id_1, face_id_2, boundaries);
}

//Get the mesh boundaries
//vecs, face_id_0/1/2: points of mesh
//boundaries: output boundaries
//inside: return a inside vector among the input surface
void CGAL_3D_Triangle_Mesh_Boundary(Vector3d1 &vecs, std::vector<int> &face_id_0, std::vector<int> &face_id_1, std::vector<int> &face_id_2, Vector3d2 &boundaries, Vector3d &inside)
{
	Vector3d2 segments;

	std::vector<std::vector<int>> vecs_neigbor(vecs.size(), std::vector<int>());
	std::vector<std::vector<int>> vecs_neigbor_lable(vecs.size(), std::vector<int>());
	std::vector<int> edges;
	for (int i = 0; i < face_id_0.size(); i++){
		int index_0 = face_id_0[i];
		int index_1 = face_id_1[i];
		int index_2 = face_id_2[i];
		edges.push_back(index_0);
		edges.push_back(index_1);
		edges.push_back(index_1);
		edges.push_back(index_2);
		edges.push_back(index_2);
		edges.push_back(index_0);

		edges.push_back(index_1);
		edges.push_back(index_0);
		edges.push_back(index_2);
		edges.push_back(index_1);
		edges.push_back(index_0);
		edges.push_back(index_2);
	}

	for (int i = 0; i < edges.size(); i = i + 2){
		int index_0 = edges[i];
		int index_1 = edges[i + 1];

		int search_0 = -1;
		for (int j = 0; j < vecs_neigbor[index_0].size() && search_0<0; j++){
			if (vecs_neigbor[index_0][j] == index_1){
				search_0 = j;
				vecs_neigbor_lable[index_0][j]++;
			}
		}

		if (search_0 < 0){
			vecs_neigbor[index_0].push_back(index_1);
			vecs_neigbor_lable[index_0].push_back(1);
		}
	}

	for (int i = 0; i < vecs.size(); i++)
	{
		int index_0 = i;

		for (int j = 0; j < vecs_neigbor_lable[i].size(); j++)
		{
			if (vecs_neigbor_lable[i][j] == 1)
			{
				int index_1 = vecs_neigbor[i][j];

				Vector3d1 segment;
				segment.push_back(vecs[index_0]);
				segment.push_back(vecs[index_1]);

				segments.push_back(segment);

				//delete
				for (int k = 0; k < vecs_neigbor[index_1].size(); k++)
				{
					if (vecs_neigbor[index_1][k] == index_0)
					{
						vecs_neigbor_lable[index_1][k] = 0;
					}
				}

			}
		}
	}

	CGAL_3D_Connecting_Segments(segments, boundaries);

	//CGAL_3D_Triangel_Mesh_Most_Inside_Point(vecs,face_id_0,face_id_1,face_id_2,inside);

	for (int i = 0; i < segments.size(); i++)
		Vector3d1().swap(segments[i]);
	Vector3d2().swap(segments);
	std::vector<std::vector<int>>().swap(vecs_neigbor);
	std::vector<std::vector<int>>().swap(vecs_neigbor_lable);
	std::vector<int>().swap(edges);
}

void CGAL_3D_Triangel_Mesh_Most_Inside_Point(Vector3d1 &vecs, std::vector<int> &face_id_0, std::vector<int> &face_id_1, std::vector<int> &face_id_2, Vector3d &inside)
{
	//inside
	std::vector<bool> vec_boundary;
	CGAL_3D_Triangle_Mesh_Boundary(vecs, face_id_0, face_id_1, face_id_2, vec_boundary);

	kdtree *tree = kd_create(3);

	for (int i = 0; i <vecs.size(); i++)
	{
		if (vec_boundary[i])
		{
			void* val = &vecs[i];
			kd_insert3(tree, vecs[i][0], vecs[i][1], vecs[i][2], val);
		}
	}

	double max_d = 0.0;
	for (int i = 0; i <vecs.size(); i++)
	{
		if (!vec_boundary[i])
		{
			double *pos = new double[3];
			pos[0] = vecs[i][0];
			pos[1] = vecs[i][1];
			pos[2] = vecs[i][2];
			struct kdres *r = kd_nearest(tree, pos);
			double position[3];
			*(int*)kd_res_item(r, position);
			double d = CGAL_3D_Distance_Point_Point(pos[0], pos[1], pos[2], position[0], position[1], position[2]);
			if (d > max_d)
			{
				max_d = d;
				inside = vecs[i];
			}
		}
	}

	std::vector<bool>().swap(vec_boundary);
}

//Get the mesh boundaries
//path: input obj/off path (Note that input model should have open boundaries.)
//boundaries: output boundaries
//inside: return a inside vector among the input surface
void CGAL_3D_Triangle_Mesh_Boundary(std::string path, Vector3d2 &boundaries, Vector3d &inside)
{
	Vector3d1 vecs;
	std::vector<int> face_id_0;
	std::vector<int> face_id_1;
	std::vector<int> face_id_2;
	CGAL_3D_Read_Triangle_Mesh(path, vecs, face_id_0, face_id_1, face_id_2);
	CGAL_3D_Triangle_Mesh_Boundary(vecs, face_id_0, face_id_1, face_id_2, boundaries,inside);
}


//compute 3d Convex Hulls
/***************************************************************************************************/

void CGAL_2D_Convex_Hulls(Vector2d1& vec, Vector2d1& hull_points)
{
	std::vector<Point_2> points;
	std::vector<Point_2> results;

	for (int i = 0; i < vec.size(); i++)
		points.push_back(VectorPoint2d(vec[i]));

	CGAL::convex_hull_2(points.begin(), points.end(), std::back_inserter(results));
	std::cout << results.size() << " points on the convex hull" << std::endl;

	for (int i = 0; i < results.size(); i++)
		hull_points.push_back(PointVector2d(results[i]));
}

void CGAL_2D_OBB_Box(Vector2d1& vec,
	Vector2d &center, Vector2d &axis_0, Vector2d &axis_1, double &entent_0, double &entent_1)
{
	Wm5::Vector2<float> *points = new Wm5::Vector2<float>[vec.size()];

	for (int i = 0; i<vec.size(); i++)
	{
		points[i].X() = vec[i][0];
		points[i].Y() = vec[i][1];;
	}
	Wm5::Box2<float> b = Wm5::ContOrientedBox(vec.size(), points);
	
	center[0] = b.Center.X();
	center[1] = b.Center.Y();

	axis_0[0] = b.Axis[0].X();
	axis_0[1] = b.Axis[0].Y();

	axis_1[0] = b.Axis[1].X();
	axis_1[1] = b.Axis[1].Y();

	entent_0 = b.Extent[0];
	entent_1 = b.Extent[1];
}

/***************************************************************************************************/


//compute 3d Convex Hulls
/***************************************************************************************************/

//a functor computing the plane containing a triangular facet
struct Plane_from_facet {
	Polyhedron_3::Plane_3 operator()(Polyhedron_3::Facet& f) {
		Halfedge_handle h = f.halfedge();
		return Polyhedron_3::Plane_3(h->vertex()->point(), h->next()->vertex()->point(),h->opposite()->vertex()->point());
	}
};


void CGAL_3D_Convex_Hulls(Vector3d1 &vec, Vector3d1 &hull_points)
{
	std::vector<Point_3> points;
	for (int i = 0; i < vec.size(); i++)
		points.push_back(VectorPoint3d(vec[i]));

	if (points.size() <= 3)
	{
		hull_points = vec;
		return;
	}
	// define polyhedron to hold convex hull
	Polyhedron_3 poly;
	// compute convex hull of non-collinear points
	CGAL::convex_hull_3(points.begin(), points.end(), poly);
	std::cout << "The convex hull contains " << poly.size_of_vertices() << " vertices" << std::endl;
	for (Polyhedron_3::Vertex_iterator iter = poly.vertices_begin(); iter != poly.vertices_end(); iter++)
	{
		Point_3 p = iter->point();
		hull_points.push_back(PointVector3d(p));
	}
}

void CGAL_3D_Convex_Hulls(Vector3d1 &vec, Vector3d1 &hull_points,
	std::vector<int>& hulls_surface_0, std::vector<int>& hulls_surface_1, std::vector<int>& hulls_surface_2)
{
	std::vector<R3> pts, pts2;
	R3 pt;
	pts.clear();

	for (int i = 0; i < vec.size(); i++)
	{
		//points.push_back(Point_3(xs[i], ys[i], zs[i]));
		pt.id = i;
		pt.r = vec[i][0];
		pt.c = vec[i][1];
		pt.z = vec[i][2];
		pts.push_back(pt);
	}

	std::vector<Tri> tris;

	std::vector<int> outx;
	int nx = de_duplicateR3(pts, outx, pts2);
	pts.clear();
	int ts = NewtonApple_hull_3D(pts2, tris);

	for (int i = 0; i < (int)tris.size(); i++)
	{
		pts.push_back(pts2[tris[i].a]);
		pts.push_back(pts2[tris[i].b]);
		pts.push_back(pts2[tris[i].c]);
	}
	pts2.clear();
	outx.clear();
	tris.clear();
	nx = de_duplicateR3(pts, outx, pts2);
	ts = NewtonApple_hull_3D(pts2, tris);


	for (int i = 0; i < pts2.size(); i++)
	{
		hull_points.push_back(Vector3d(pts2[i].r, pts2[i].c, pts2[i].z));
	}

	for (int i = 0; i < (int)tris.size(); i++)
	{
		hulls_surface_0.push_back(tris[i].a);
		hulls_surface_1.push_back(tris[i].b);
		hulls_surface_2.push_back(tris[i].c);
	}
	pts2.clear();

}


void CGAL_3D_Convex_Hulls(Vector3d1 &vec, Vector3d1 &hull_points,
	Vector3d1 &plane_p, Vector3d1 &plane_n)
{
	std::vector<R3> pts, pts2;
	R3 pt;
	pts.clear();

	for (int i = 0; i < vec.size(); i++)
	{
		//points.push_back(Point_3(xs[i], ys[i], zs[i]));
		pt.id = i;
		pt.r = vec[i][0];
		pt.c = vec[i][1];
		pt.z = vec[i][2];
		pts.push_back(pt);
	}

	std::vector<Tri> tris;

	std::vector<int> outx;
	int nx = de_duplicateR3(pts, outx, pts2);
	pts.clear();

	//int ts = NewtonApple_Delaunay( pts2, tris);
	int ts = NewtonApple_hull_3D(pts2, tris);

	int nr = (int)tris.size();

	for (int i = 0; i < nr; i++)
	{
		pts.push_back(pts2[tris[i].a]);
		pts.push_back(pts2[tris[i].b]);
		pts.push_back(pts2[tris[i].c]);

		Point_3 p_0 = Point_3(pts2[tris[i].a].r, pts2[tris[i].a].c, pts2[tris[i].a].z);
		Point_3 p_1 = Point_3(pts2[tris[i].b].r, pts2[tris[i].b].c, pts2[tris[i].b].z);
		Point_3 p_2 = Point_3(pts2[tris[i].c].r, pts2[tris[i].c].c, pts2[tris[i].c].z);

		Plane_3 plane(p_0, p_1, p_2);

		Point_3 p((p_0[0] + p_1[0] + p_2[0]) / 3.0, (p_0[1] + p_1[1] + p_2[1]) / 3.0, (p_0[2] + p_1[2] + p_2[2]) / 3.0);
		plane_p.push_back(PointVector3d(p));

		Vector_3 v = plane.orthogonal_vector();
		plane_n.push_back(Vector3d(v[0], v[1], v[2]));
	}
	pts2.clear();
	nx = de_duplicateR3(pts, outx, pts2);
}


void CGAL_3D_Convex_Hulls(Vector3d1 &vec, Vector3d1 &hull_points,
	std::vector<int>& hulls_surface_0, std::vector<int>& hulls_surface_1, std::vector<int>& hulls_surface_2,
	Vector3d1 &plane_p, Vector3d1 &plane_n)
{
	std::vector<R3> pts, pts2;
	R3 pt;
	pts.clear();

	for (int i = 0; i < vec.size(); i++)
	{
		//points.push_back(Point_3(xs[i], ys[i], zs[i]));
		pt.id = i;
		pt.r = vec[i][0];
		pt.c = vec[i][1];
		pt.z = vec[i][2];
		pts.push_back(pt);
	}

	std::vector<Tri> tris;

	std::vector<int> outx;
	int nx = de_duplicateR3(pts, outx, pts2);
	pts.clear();
	int ts = NewtonApple_hull_3D(pts2, tris);

	for (int i = 0; i < (int)tris.size(); i++)
	{
		pts.push_back(pts2[tris[i].a]);
		pts.push_back(pts2[tris[i].b]);
		pts.push_back(pts2[tris[i].c]);
	}
	pts2.clear();
	outx.clear();
	tris.clear();
	nx = de_duplicateR3(pts, outx, pts2);
	ts = NewtonApple_hull_3D(pts2, tris);


	for (int i = 0; i < pts2.size(); i++)
	{
		hull_points.push_back(Vector3d(pts2[i].r, pts2[i].c, pts2[i].z));
	}

	for (int i = 0; i < (int)tris.size(); i++)
	{
		hulls_surface_0.push_back(tris[i].a);
		hulls_surface_1.push_back(tris[i].b);
		hulls_surface_2.push_back(tris[i].c);

		Point_3 p_0 = Point_3(pts2[tris[i].a].r, pts2[tris[i].a].c, pts2[tris[i].a].z);
		Point_3 p_1 = Point_3(pts2[tris[i].b].r, pts2[tris[i].b].c, pts2[tris[i].b].z);
		Point_3 p_2 = Point_3(pts2[tris[i].c].r, pts2[tris[i].c].c, pts2[tris[i].c].z);

		Plane_3 plane(p_0, p_1, p_2);

		Point_3 p((p_0[0] + p_1[0] + p_2[0]) / 3.0, (p_0[1] + p_1[1] + p_2[1]) / 3.0, (p_0[2] + p_1[2] + p_2[2]) / 3.0);
		//Point_3 p = plane.point();

		plane_p.push_back(PointVector3d(p));

		Vector_3 v = plane.orthogonal_vector();
		plane_n.push_back(Vector3d(v[0],v[1],v[2]));
	}
	pts2.clear();
}


Vector3d CGAL_3D_Projection_Point_Plane(Vector3d p, Vector3d plane_p, Vector3d plane_n)
{
	Plane_3 plane(Point_3(plane_p[0], plane_p[1], plane_p[2]), Vector_3(plane_n[0], plane_n[1], plane_n[2]));
	Point_3 project = plane.projection(Point_3(p[0], p[1], p[1]));
	Vector3d result(project[0], project[1], project[2]);
	return result;
}
Vector3d CGAL_3D_Projection_Point_Plane(Vector3d p, Vector3d plane_p_0, Vector3d plane_p_1, Vector3d plane_p_2)
{
	Point_3 p0 = VectorPoint3d(plane_p_0);
	Point_3 p1 = VectorPoint3d(plane_p_1);
	Point_3 p2 = VectorPoint3d(plane_p_2);

	Plane_3 plane(p1, CGAL::cross_product(p2 - p1, p0 - p1));
	Point_3 project = plane.projection(VectorPoint3d(p));
	Vector3d result(project[0], project[1], project[2]);
	return result;

}

Vector2d CGAL_3D_Projection_3D_Point_Plane_2D(Vector3d p, Vector3d plane_p, Vector3d plane_n)
{
	Plane_3 plane(VectorPoint3d(plane_p), Vector_3(plane_n[0], plane_n[1], plane_n[2]));
	Point_2 r = point_to_2d(VectorPoint3d(p), plane);
	return Vector2d(r[0],r[1]);
}


void CGAL_3D_Plane_ABCD(Vector3d plane_p, Vector3d plane_n, double &a, double &b, double &c, double &d)
{
	Plane_3 plane(VectorPoint3d(plane_p), Vector_3(plane_n[0], plane_n[1], plane_n[2]));
	a = plane.a();
	b = plane.b();
	c = plane.c();
	d = plane.d();
}

Vector2d CGAL_3D_Projection_3D_Point_Plane_2D(Vector3d p, Vector3d plane_p_0, Vector3d plane_p_1, Vector3d plane_p_2)
{
	Point_3 p0 = VectorPoint3d(plane_p_0);
	Point_3 p1 = VectorPoint3d(plane_p_1);
	Point_3 p2 = VectorPoint3d(plane_p_2);

	Plane_3 plane(p1, CGAL::cross_product(p2 - p1, p0 - p1));
	Point_2 r = point_to_2d(VectorPoint3d(p), plane);
	return Vector2d(r[0], r[1]);
}

/***************************************************************************************************/

//3d plane relaed
/***************************************************************************************************/
Vector3d CGAL_3D_Plane_Base_1(Vector3d plane_p, Vector3d plane_n)
{
	Plane_3 plane(VectorPoint3d(plane_p), Vector_3(plane_n[0], plane_n[1], plane_n[2]));
	Vector_3 v = plane.base1();
	return Vector3d(v[0],v[1],v[2]);
}

/***************************************************************************************************/

//3D mesh curvature
/***************************************************************************************************/

void CGAL_3D_Mesh_Curvature(Vector3d1& vecs, std::vector<int>& face_id_0, std::vector<int>& face_id_1,
	std::vector<int>& face_id_2, std::vector<double> &max_curs, std::vector<double> &min_curs)
{
	int verticeSize = vecs.size();
	int faceindiceSize = face_id_0.size() * 3;

	Wm5::Vector3<double> *points = new Wm5::Vector3<double>[verticeSize];
	int *indices = new int[faceindiceSize];

	for (int i = 0; i<verticeSize; i++)
	{
		points[i].X() = vecs[i][0];
		points[i].Y() = vecs[i][1];
		points[i].Z() = vecs[i][2];
	}

	for (int i = 0; i<face_id_0.size(); i++)
	{
		indices[3 * i] = face_id_0[i];
		indices[3 * i + 1] = face_id_1[i];
		indices[3 * i + 2] = face_id_2[i];
	}
	Wm5::MeshCurvature<double> meshCurvature(verticeSize, points, faceindiceSize / 3, indices);

	for (int i = 0; i < verticeSize; i++)
	{
		max_curs.push_back(meshCurvature.GetMaxCurvatures()[i]);
		min_curs.push_back(meshCurvature.GetMinCurvatures()[i]);
	}
}

void CGAL_3D_Mesh_Curvature(Vector3d1& vecs, std::vector<std::vector<int>>& face_ids, std::vector<double> &max_curs, std::vector<double> &min_curs)
{
	std::vector<int> face_id_0, face_id_1, face_id_2;

	for (int i = 0; i < face_ids.size(); i++)
	{
		face_id_0.push_back(face_ids[i][0]);
		face_id_1.push_back(face_ids[i][1]);
		face_id_2.push_back(face_ids[i][2]);
	}

	CGAL_3D_Mesh_Curvature(vecs, face_id_0, face_id_1, face_id_2, max_curs, min_curs);

	std::vector<int>().swap(face_id_0);
	std::vector<int>().swap(face_id_1);
	std::vector<int>().swap(face_id_2);

}

void CGAL_3D_Mesh_Curvature(Vector3d1& vecs, std::vector<int>& face_id_0, std::vector<int>& face_id_1,
	std::vector<int>& face_id_2, std::vector<double> &max_curs, std::vector<double> &min_curs,
	Vector3d1 &max_curs_directions, Vector3d1 &min_curs_directions)
{
	int verticeSize = vecs.size();
	int faceindiceSize = face_id_0.size() * 3;

	Wm5::Vector3<double> *points = new Wm5::Vector3<double>[verticeSize];
	int *indices = new int[faceindiceSize];

	for (int i = 0; i<verticeSize; i++)
	{
		points[i].X() = vecs[i][0];
		points[i].Y() = vecs[i][1];
		points[i].Z() = vecs[i][2];
	}

	for (int i = 0; i<face_id_0.size(); i++)
	{
		indices[3 * i] = face_id_0[i];
		indices[3 * i + 1] = face_id_1[i];
		indices[3 * i + 2] = face_id_2[i];
	}

	Wm5::MeshCurvature<double> meshCurvature(verticeSize, points, faceindiceSize / 3, indices);

	for (int i = 0; i < verticeSize; i++)
	{
		max_curs.push_back(meshCurvature.GetMaxCurvatures()[i]);
		min_curs.push_back(meshCurvature.GetMinCurvatures()[i]);

		Wm5::Vector3<double> max_curs_direction = meshCurvature.GetMaxDirections()[i];
		Wm5::Vector3<double> min_curs_direction = meshCurvature.GetMinDirections()[i];

		max_curs_directions.push_back(Vector3d(max_curs_direction[0], max_curs_direction[1], max_curs_direction[2]));
		min_curs_directions.push_back(Vector3d(min_curs_direction[0], min_curs_direction[1], min_curs_direction[2]));
	}
}
void CGAL_3D_Mesh_Curvature(Vector3d1& vecs, std::vector<std::vector<int>>& face_ids, std::vector<double> &max_curs, std::vector<double> &min_curs,
	Vector3d1 &max_curs_directions, Vector3d1 &min_curs_directions)
{
	std::vector<int> face_id_0, face_id_1, face_id_2;

	for (int i = 0; i < face_ids.size(); i++)
	{
		face_id_0.push_back(face_ids[i][0]);
		face_id_1.push_back(face_ids[i][1]);
		face_id_2.push_back(face_ids[i][2]);
	}

	CGAL_3D_Mesh_Curvature(vecs, face_id_0, face_id_1, face_id_2, max_curs, min_curs, max_curs_directions, min_curs_directions);

	std::vector<int>().swap(face_id_0);
	std::vector<int>().swap(face_id_1);
	std::vector<int>().swap(face_id_2);
}


void CGAL_3D_Mesh_Curvature(Vector3d1& vecs, std::vector<int>& face_id_0, std::vector<int>& face_id_1,
	std::vector<int>& face_id_2, std::vector<double> &max_curs, std::vector<double> &min_curs,
	Vector3d1 &max_curs_directions, Vector3d1 &min_curs_directions, Vector3d1 &normals)
{
	int verticeSize = vecs.size();
	int faceindiceSize = face_id_0.size() * 3;

	Wm5::Vector3<double> *points = new Wm5::Vector3<double>[verticeSize];
	int *indices = new int[faceindiceSize];

	for (int i = 0; i<verticeSize; i++)
	{
		points[i].X() = vecs[i][0];
		points[i].Y() = vecs[i][1];
		points[i].Z() = vecs[i][2];
	}

	for (int i = 0; i<face_id_0.size(); i++)
	{
		indices[3 * i] = face_id_0[i];
		indices[3 * i + 1] = face_id_1[i];
		indices[3 * i + 2] = face_id_2[i];
	}
	Wm5::MeshCurvature<double> meshCurvature(verticeSize, points, faceindiceSize / 3, indices);

	for (int i = 0; i < verticeSize; i++)
	{
		max_curs.push_back(meshCurvature.GetMaxCurvatures()[i]);
		min_curs.push_back(meshCurvature.GetMinCurvatures()[i]);

		Wm5::Vector3<double> max_curs_direction = meshCurvature.GetMaxDirections()[i];
		Wm5::Vector3<double> min_curs_direction = meshCurvature.GetMinDirections()[i];

		max_curs_directions.push_back(Vector3d(max_curs_direction[0], max_curs_direction[1], max_curs_direction[2]));
		min_curs_directions.push_back(Vector3d(min_curs_direction[0], min_curs_direction[1], min_curs_direction[2]));

		Wm5::Vector3<double> normal = meshCurvature.GetNormals()[i];
		normals.push_back(Vector3d(normal[0], normal[1], normal[2]));
	}
}

void CGAL_3D_Mesh_Curvature(Vector3d1& vecs, std::vector<std::vector<int>>& face_ids, std::vector<double> &max_curs, std::vector<double> &min_curs,
	Vector3d1 &max_curs_directions, Vector3d1 &min_curs_directions, Vector3d1 &normals)
{
	std::vector<int> face_id_0, face_id_1, face_id_2;

	for (int i = 0; i < face_ids.size(); i++)
	{
		face_id_0.push_back(face_ids[i][0]);
		face_id_1.push_back(face_ids[i][1]);
		face_id_2.push_back(face_ids[i][2]);
	}

	CGAL_3D_Mesh_Curvature(vecs, face_id_0, face_id_1, face_id_2, max_curs, min_curs, max_curs_directions, min_curs_directions, normals);

	std::vector<int>().swap(face_id_0);
	std::vector<int>().swap(face_id_1);
	std::vector<int>().swap(face_id_2);
}

void CGAL_3D_Mesh_Normal(Vector3d1 &ps, std::vector<std::vector<int>> &face_ids, Vector3d1 &normals)
{
	int verticeSize = ps.size();
	int faceindiceSize = face_ids.size() * 3;

	Wm5::Vector3<double> *points = new Wm5::Vector3<double>[verticeSize];
	int *indices = new int[faceindiceSize];

	for (int i = 0; i<verticeSize; i++)
	{
		points[i].X() = ps[i][0];
		points[i].Y() = ps[i][1];
		points[i].Z() = ps[i][2];
	}
	for (int i = 0; i<face_ids.size(); i++)
	{
		indices[3 * i] = face_ids[i][0];
		indices[3 * i + 1] = face_ids[i][1];
		indices[3 * i + 2] = face_ids[i][2];
	}
	Wm5::MeshCurvature<double> meshCurvature(verticeSize, points, faceindiceSize / 3, indices);

	for (int i = 0; i < verticeSize; i++)
	{
		Wm5::Vector3<double> normal = meshCurvature.GetNormals()[i];
		normals.push_back(Vector3d(normal[0], normal[1], normal[2]));
	}
}

void CGAL_3D_Mesh_Normal(Vector3d1 &ps, std::vector<int>& face_id_0, std::vector<int>& face_id_1, std::vector<int>& face_id_2, Vector3d1 &normals)
{
	int verticeSize = ps.size();
	int faceindiceSize = face_id_0.size() * 3;

	Wm5::Vector3<double> *points = new Wm5::Vector3<double>[verticeSize];
	int *indices = new int[faceindiceSize];

	for (int i = 0; i<verticeSize; i++)
	{
		points[i].X() = ps[i][0];
		points[i].Y() = ps[i][1];
		points[i].Z() = ps[i][2];
	}
	for (int i = 0; i<face_id_0.size(); i++)
	{
		indices[3 * i] = face_id_0[i];
		indices[3 * i + 1] = face_id_1[i];
		indices[3 * i + 2] = face_id_2[i];
	}
	Wm5::MeshCurvature<double> meshCurvature(verticeSize, points, faceindiceSize / 3, indices);

	for (int i = 0; i < verticeSize; i++)
	{
		Wm5::Vector3<double> normal = meshCurvature.GetNormals()[i];
		normals.push_back(Vector3d(normal[0], normal[1], normal[2]));
	}
}

Vector3d CGAL_3D_Mesh_Center(const Vector3d2 &ps)
{
	Vector3d1 ps1;
	for (auto ps_ : ps)for (auto p : ps_)ps1.emplace_back(p);
	return CGAL_3D_Mesh_Center(ps1);
}

Vector3d CGAL_3D_Mesh_Center(Vector3d1 &ps)
{
	Vector3d center(0.0,0.0,0.0);
	
	for (int i = 0; i < ps.size(); i++)
	{
		center += ps[i];
	}

	center = center / (float)ps.size();
	return center;
}
void CGAL_3D_Mesh_Boundingbox(const Vector3d2 &ps, Vector3d &min_corner, Vector3d &max_corner)
{
	Vector3d1 ps1;
	for (auto ps_ : ps)for (auto p : ps_)ps1.emplace_back(p);
	CGAL_3D_Mesh_Boundingbox(ps1, min_corner, max_corner);
}

void CGAL_3D_Mesh_Boundingbox(Vector3d1 &ps, Vector3d &min_corner, Vector3d &max_corner)
{
	min_corner = ps[0];
	max_corner = ps[0];
	for (int i = 0; i < ps.size(); i++)
	{
		min_corner[0] = std::min(min_corner[0], ps[i][0]);
		min_corner[1] = std::min(min_corner[1], ps[i][1]);
		min_corner[2] = std::min(min_corner[2], ps[i][2]);
		max_corner[0] = std::max(max_corner[0], ps[i][0]);
		max_corner[1] = std::max(max_corner[1], ps[i][1]);
		max_corner[2] = std::max(max_corner[2], ps[i][2]);
	}
}

//
//void CGAL_3D_Mesh_Curvature(std::vector<double>& xs, std::vector<double>& ys, std::vector<double>& zs,
//	std::vector<int>& face_id_0, std::vector<int>& face_id_1, std::vector<int>& face_id_2, std::vector<double> &max_curs, std::vector<double> &min_curs,
//	std::vector<double> &max_curs_directions_x, std::vector<double> &max_curs_directions_y, std::vector<double> &max_curs_directions_z,
//	std::vector<double> &min_curs_directions_x, std::vector<double> &min_curs_directions_y, std::vector<double> &min_curs_directions_z,
//	std::vector<double> &normal_x, std::vector<double> &normal_y, std::vector<double> &normal_z)
//{
//
//	int verticeSize = xs.size();
//	int faceindiceSize = face_id_0.size() * 3;
//
//	Wm5::Vector3<double> *points = new Wm5::Vector3<double>[verticeSize];
//	int *indices = new int[faceindiceSize];
//
//	for (int i = 0; i<verticeSize; i++)
//	{
//		points[i].X() = xs[i];
//		points[i].Y() = ys[i];
//		points[i].Z() = zs[i];
//	}
//
//	for (int i = 0; i<face_id_0.size(); i++)
//	{
//		indices[3 * i] = face_id_0[i];
//		indices[3 * i + 1] = face_id_1[i];
//		indices[3 * i + 2] = face_id_2[i];
//	}
//	Wm5::MeshCurvature<double> meshCurvature(verticeSize, points, faceindiceSize / 3, indices);
//
//	for (int i = 0; i < verticeSize; i++)
//	{
//		max_curs.push_back(meshCurvature.GetMaxCurvatures()[i]);
//		min_curs.push_back(meshCurvature.GetMinCurvatures()[i]);
//
//		Wm5::Vector3<double> max_curs_direction = meshCurvature.GetMaxDirections()[i];
//		Wm5::Vector3<double> min_curs_direction = meshCurvature.GetMinDirections()[i];
//
//		Wm5::Vector3<double> normal = meshCurvature.GetNormals()[i];
//
//		normal_x.push_back(normal[0]);
//		normal_y.push_back(normal[1]);
//		normal_z.push_back(normal[2]);
//
//		max_curs_directions_x.push_back(max_curs_direction[0]);
//		max_curs_directions_y.push_back(max_curs_direction[1]);
//		max_curs_directions_z.push_back(max_curs_direction[2]);
//
//		min_curs_directions_x.push_back(min_curs_direction[0]);
//		min_curs_directions_y.push_back(min_curs_direction[1]);
//		min_curs_directions_z.push_back(min_curs_direction[2]);
//	}
//}

//Vector3d max_0 = mesh_max_curs_directions[point_id_0];
//Vector3d max_1 = mesh_max_curs_directions[point_id_1];
//Vector3d max_2 = mesh_max_curs_directions[point_id_2];
//
//Vector3d min_0 = mesh_min_curs_directions[point_id_0];
//Vector3d min_1 = mesh_min_curs_directions[point_id_1];
//Vector3d min_2 = mesh_min_curs_directions[point_id_2];

void Curvature_Direction_Adjusting(Vector3d &cur_0, Vector3d &cur_1, Vector3d &cur_2)
{
	int f_i = 0;
	int f_j = 0;
	int f_k = 0;

	double angle = 1000000000.0;
	for (int i = -1; i <= 1; i = i + 2)
	{
		for (int j = -1; j <= 1; j = j + 2)
		{
			for (int k = -1; k <= 1; k = k + 2)
			{
				double angle_0 = Math::GetAngleBetween(cur_0*(float)i, cur_1*(float)j);
				double angle_1 = Math::GetAngleBetween(cur_0*(float)i, cur_2*(float)k);
				double angle_2 = Math::GetAngleBetween(cur_2*(float)k, cur_1*(float)j);

				if (angle > angle_0 + angle_1 + angle_2)
				{
					angle = angle_0 + angle_1 + angle_2;
					f_i = i;
					f_j = j;
					f_k = k;
				}
			}
		}
	}

	cur_0 = cur_0*(float)f_i;
	cur_1 = cur_1*(float)f_j;
	cur_2 = cur_2*(float)f_k;
}

void CGAL_Mesh_Field_Query(std::string path, Vector3d1 &gradients, Vector3d1& input_points, Vector3d1 &points_gradients)
{
	Polyhedron_3 polyhedron;
	Construct_Polyhedron(polyhedron, path);
	Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);
	tree.accelerate_distance_queries();

	for (int i = 0; i < input_points.size(); i++)
	{
		Poly_point_3 query(input_points[i][0], input_points[i][1], input_points[i][2]);
		Point_and_primitive_id pp = tree.closest_point_and_primitive(query);

		Poly_point_3 p0 = pp.second->halfedge()->next()->next()->vertex()->point();
		Poly_point_3 p1 = pp.second->halfedge()->vertex()->point();
		Poly_point_3 p2 = pp.second->halfedge()->next()->vertex()->point();

		int point_id_0 = pp.second->halfedge()->next()->next()->vertex()->id();
		int point_id_1 = pp.second->halfedge()->vertex()->id();
		int point_id_2 = pp.second->halfedge()->next()->vertex()->id();

		double u, v, w;
		Barycentric(query, p0, p1, p2, u, v, w);
		points_gradients.push_back((float)u*gradients[point_id_0] + (float)v*gradients[point_id_1] + (float)w*gradients[point_id_2]);
	}
}

void CGAL_Mesh_Field_Query(std::string path, std::vector<double> &gradient_values, Vector3d1& input_points, std::vector<double> &points_gradient_values)
{
	Polyhedron_3 polyhedron;
	Construct_Polyhedron(polyhedron, path);
	Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);
	tree.accelerate_distance_queries();

	for (int i = 0; i < input_points.size(); i++)
	{
		Poly_point_3 query(input_points[i][0], input_points[i][1], input_points[i][2]);
		Point_and_primitive_id pp = tree.closest_point_and_primitive(query);

		Poly_point_3 p0 = pp.second->halfedge()->next()->next()->vertex()->point();
		Poly_point_3 p1 = pp.second->halfedge()->vertex()->point();
		Poly_point_3 p2 = pp.second->halfedge()->next()->vertex()->point();

		int point_id_0 = pp.second->halfedge()->next()->next()->vertex()->id();
		int point_id_1 = pp.second->halfedge()->vertex()->id();
		int point_id_2 = pp.second->halfedge()->next()->vertex()->id();

		double u, v, w;
		Barycentric(query, p0, p1, p2, u, v, w);
		points_gradient_values.push_back((float)u*gradient_values[point_id_0] + (float)v*gradient_values[point_id_1] + (float)w*gradient_values[point_id_2]);
	}
}
void CGAL_Mesh_Field_Query(std::string path, std::vector<double> &gradient_values, Vector3d2& input_point_es, std::vector<std::vector<double>> &points_gradient_value_es)
{
	Polyhedron_3 polyhedron;
	Construct_Polyhedron(polyhedron, path);
	Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);
	tree.accelerate_distance_queries();

	for (int i = 0; i < input_point_es.size(); i++)
	{
		std::vector<double> values;
		for (int j = 0; j < input_point_es[i].size(); j++)
		{
			Poly_point_3 query(input_point_es[i][j][0], input_point_es[i][j][1], input_point_es[i][j][2]);
			Point_and_primitive_id pp = tree.closest_point_and_primitive(query);

			Poly_point_3 p0 = pp.second->halfedge()->next()->next()->vertex()->point();
			Poly_point_3 p1 = pp.second->halfedge()->vertex()->point();
			Poly_point_3 p2 = pp.second->halfedge()->next()->vertex()->point();

			int point_id_0 = pp.second->halfedge()->next()->next()->vertex()->id();
			int point_id_1 = pp.second->halfedge()->vertex()->id();
			int point_id_2 = pp.second->halfedge()->next()->vertex()->id();

			double u, v, w;
			Barycentric(query, p0, p1, p2, u, v, w);
			values.push_back((float)u*gradient_values[point_id_0] + (float)v*gradient_values[point_id_1] + (float)w*gradient_values[point_id_2]);
		}
		points_gradient_value_es.push_back(values);
	}
}

void CGAL_Curvature_Mesh(std::string path, Vector3d1& input_points, std::vector<double> &max_curs, std::vector<double> &min_curs,
	Vector3d1 &max_curs_directions, Vector3d1 &min_curs_directions)
{
	std::cout << "CGAL_Curvature_Mesh.." << std::endl;

	if (path.substr(path.size() - 3, path.size()) == "off")
	{
		Polyhedron_3 polyhedron;
		Construct_Polyhedron(polyhedron, path);
		Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);
		tree.accelerate_distance_queries();

		//get mesh vertices and surface
		//std::vector<double> mesh_xs, mesh_ys, mesh_zs;
		Vector3d1 vecs;
		std::vector<int> mesh_face_id_0, mesh_face_id_1, mesh_face_id_2;

		//CGAL_3D_Read_Triangle_Mesh(path, vecs, mesh_face_id_0, mesh_face_id_1, mesh_face_id_2);
		
		//Polyhedron_3 polyhedron;
		//Construct_Polyhedron(polyhedron, vecs, mesh_face_id_0, mesh_face_id_1, mesh_face_id_2);
		//Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);
		//tree.accelerate_distance_queries();

		for (Polyhedron_3::Vertex_iterator iter = polyhedron.vertices_begin();
			iter != polyhedron.vertices_end(); iter++)
		{
			Poly_point_3 p = iter->point();
			vecs.push_back(Vector3d(p[0], p[1], p[2]));
		}
		for (Polyhedron_3::Face_iterator iter = polyhedron.facets_begin(); iter != polyhedron.facets_end(); iter++)
		{
			mesh_face_id_0.push_back(iter->halfedge()->next()->next()->vertex()->id());
			mesh_face_id_1.push_back(iter->halfedge()->vertex()->id());
			mesh_face_id_2.push_back(iter->halfedge()->next()->vertex()->id());
		}

		//compute surface normals
		int verticeSize = vecs.size();
		int faceindiceSize = mesh_face_id_0.size() * 3;

		Wm5::Vector3<double> *points = new Wm5::Vector3<double>[verticeSize];
		int *indices = new int[faceindiceSize];

		for (int i = 0; i<verticeSize; i++)
		{
			points[i].X() = vecs[i][0];
			points[i].Y() = vecs[i][1];
			points[i].Z() = vecs[i][2];
		}

		for (int i = 0; i<mesh_face_id_0.size(); i++)
		{
			indices[3 * i] = mesh_face_id_0[i];
			indices[3 * i + 1] = mesh_face_id_1[i];
			indices[3 * i + 2] = mesh_face_id_2[i];
		}
		Wm5::MeshCurvature<double> meshCurvature(verticeSize, points, faceindiceSize / 3, indices);

		std::vector<double> mesh_max_curs;
		std::vector<double> mesh_min_curs;

		Vector3d1 mesh_max_curs_directions;
		Vector3d1 mesh_min_curs_directions;

		for (int i = 0; i < verticeSize; i++)
		{
			mesh_max_curs.push_back(meshCurvature.GetMaxCurvatures()[i]);
			mesh_min_curs.push_back(meshCurvature.GetMinCurvatures()[i]);

			Wm5::Vector3<double> max_curs_direction = meshCurvature.GetMaxDirections()[i];
			Wm5::Vector3<double> min_curs_direction = meshCurvature.GetMinDirections()[i];

			mesh_max_curs_directions.push_back(Vector3d(max_curs_direction[0], max_curs_direction[1], max_curs_direction[2]));
			mesh_min_curs_directions.push_back(Vector3d(min_curs_direction[0], min_curs_direction[1], min_curs_direction[2]));
		}

		for (int i = 0; i < input_points.size(); i++)
		{
			Poly_point_3 query(input_points[i][0], input_points[i][1], input_points[i][2]);
			Point_and_primitive_id pp = tree.closest_point_and_primitive(query);

			Poly_point_3 p0 = pp.second->halfedge()->next()->next()->vertex()->point();
			Poly_point_3 p1 = pp.second->halfedge()->vertex()->point();
			Poly_point_3 p2 = pp.second->halfedge()->next()->vertex()->point();

			int point_id_0 = pp.second->halfedge()->next()->next()->vertex()->id();
			int point_id_1 = pp.second->halfedge()->vertex()->id();
			int point_id_2 = pp.second->halfedge()->next()->vertex()->id();

			double u, v, w;
			Barycentric(query, p0, p1, p2, u, v, w);

			max_curs.push_back(u*mesh_max_curs[point_id_0] + v*mesh_max_curs[point_id_1] + w*mesh_max_curs[point_id_2]);
			min_curs.push_back(u*mesh_min_curs[point_id_0] + v*mesh_min_curs[point_id_1] + w*mesh_min_curs[point_id_2]);

			Curvature_Direction_Adjusting(mesh_max_curs_directions[point_id_0], mesh_max_curs_directions[point_id_1], mesh_max_curs_directions[point_id_2]);
			Curvature_Direction_Adjusting(mesh_min_curs_directions[point_id_0], mesh_min_curs_directions[point_id_1], mesh_min_curs_directions[point_id_2]);

			Vector3d max_curs_direction = (float)u*mesh_max_curs_directions[point_id_0] + (float)v*mesh_max_curs_directions[point_id_1] + (float)w*mesh_max_curs_directions[point_id_2];
			Vector3d min_curs_direction = (float)u*mesh_min_curs_directions[point_id_0] + (float)v*mesh_min_curs_directions[point_id_1] + (float)w*mesh_min_curs_directions[point_id_2];

			max_curs_directions.push_back(max_curs_direction);
			min_curs_directions.push_back(min_curs_direction);
		}
	}
}

//grid 1/0 
/***************************************************************************************************/
struct GridIndex{
	int i;
	int j;
	GridIndex(int i0, int j0){
		i = i0;
		j = j0;
	}
	GridIndex(){
	}
};

struct GridEdge{
	int s_i, s_j;
	int e_i, e_j;
	int type;
	GridEdge(int s_i_0, int s_j_0, int e_i_0, int e_j_0, int t){
		s_i = s_i_0;
		s_j = s_j_0;
		e_i = e_i_0;
		e_j = e_j_0;
		type = t;
	}
};

struct GridEdgeRelation{
	std::vector<int> ids;
};

int GetIndex(std::vector<GridEdge> &grid_edges, GridIndex i_0, GridIndex i_1)
{
	int s_i = i_0.i;
	int s_j = i_0.j;
	int e_i = i_1.i;
	int e_j = i_1.j;
	for (int i = 0; i < grid_edges.size(); i++){
		if (grid_edges[i].s_i == s_i&&grid_edges[i].s_j == s_j&&grid_edges[i].e_i == e_i&&grid_edges[i].e_j == e_j){
			return i;
		}
		if (grid_edges[i].s_i == e_i&&grid_edges[i].s_j == e_j&&grid_edges[i].e_i == s_i&&grid_edges[i].e_j == s_j){
			return i;
		}
	}
	return -1;
}

void Decomposition_Mapping(std::vector<GridIndex> &one_boundary,
	std::vector<std::vector<double>> &boundary_xs, std::vector<std::vector<double>> &boundary_ys, bool remove_b=false)
{
	//remove multi grid index
	if (remove_b){
		std::vector<GridIndex> new_one_boundary(1, one_boundary[0]);
		for (int i = 0; i < one_boundary.size(); i++){
			GridIndex &last = new_one_boundary[new_one_boundary.size() - 1];
			if (!(one_boundary[i].i == last.i&&one_boundary[i].j == last.j)){
				new_one_boundary.push_back(one_boundary[i]);
			}
		}
		std::vector<GridIndex>().swap(one_boundary);
		one_boundary = new_one_boundary;
		std::vector<GridIndex>().swap(new_one_boundary);
	}

	if (one_boundary.size() <= 2) return;

	//find repeating index
	bool run = false;
	int index_0 = -1;
	int index_1 = -1;
	for (int i = 0; i < one_boundary.size()-1&&!run; i++){
		GridIndex &current = one_boundary[i];
		for (int j = i + 1; j < one_boundary.size() && !run; j++){
			GridIndex &next = one_boundary[j];
			if (current.i == next.i&&current.j == next.j){
				index_0 = i;
				index_1 = j;
				run = true;
			}
		}
	}

	if (!run){
		std::vector<double> xs;
		std::vector<double> ys;
		for (int i = 0; i < one_boundary.size(); i++){
			xs.push_back(one_boundary[i].i - 1.0);
			ys.push_back(one_boundary[i].j - 1.0);
		}
		boundary_xs.push_back(xs);
		boundary_ys.push_back(ys);
		std::vector<double>().swap(xs);
		std::vector<double>().swap(ys);
	}
	else{
		std::vector<GridIndex> one_boundary_0;
		std::vector<GridIndex> one_boundary_1;
		for (int i = index_0; i < index_1; i++){
			one_boundary_0.push_back(one_boundary[i]);
		}
		for (int i = index_1; i < one_boundary.size(); i++){
			one_boundary_1.push_back(one_boundary[i]);
		}
		for (int i = 0; i < index_0; i++){
			one_boundary_1.push_back(one_boundary[i]);
		}
		Decomposition_Mapping(one_boundary_0, boundary_xs, boundary_ys);
		Decomposition_Mapping(one_boundary_1, boundary_xs, boundary_ys);
	}

}

void CGAL_Image_Grid_Decomposition(std::vector<std::vector<int>> &image, std::vector<std::vector<double>> &boundary_xs, std::vector<std::vector<double>> &boundary_ys)
{
	//refine  lables
	std::vector<std::vector<int>> grid;
	for (int i = 0; i < image.size() + 2; i++){
		std::vector<int> one_grid_raw;
		for (int j = 0; j < image[0].size() + 2; j++){
			if (i >= 1 && i < image.size() + 1 && j >= 1 && j < image.size() + 1){
				one_grid_raw.push_back(image[i - 1][j - 1]);
			}
			else{
				one_grid_raw.push_back(0);
			}
		}
		grid.push_back(one_grid_raw);
	}


	//std::ofstream lable_file("D:\\123.lable");
	//for (int j = 0; j < grid.size(); j++)
	//{
	//	for (int k = 0; k < grid[j].size(); k++)
	//	{
	//		lable_file << grid[j][k] << " ";
	//	}
	//	lable_file << "" << std::endl;
	//}
	//lable_file.clear();
	//lable_file.close();


	//get grid edges
	std::vector<GridEdge> grid_edges;
	for (int i = 0; i < grid.size(); i++){
		for (int j = 0; j < grid[i].size(); j++){
			if (i + 1 < grid.size()){
				if (grid[i][j] != grid[i + 1][j]){
					grid_edges.push_back(GridEdge(i, j, i + 1, j, 0));
				}
			}
			if (j + 1 < grid[i].size()){
				if (grid[i][j] != grid[i][j + 1]){
					grid_edges.push_back(GridEdge(i, j, i, j + 1, 1));
				}
			}
		}
	}

	//connect cut edges
	//GridIndex 1 2
	//GridIndex s e
	//GridIndex 3 4
	std::vector<GridEdgeRelation> grid_relations;
	std::vector<bool> grid_edges_used;
	for (int i = 0; i < grid_edges.size(); i++){
		grid_edges_used.push_back(false);
		GridIndex index_1;
		GridIndex index_2;
		GridIndex index_3;
		GridIndex index_4;
		GridIndex index_s(grid_edges[i].s_i, grid_edges[i].s_j);
		GridIndex index_e(grid_edges[i].e_i, grid_edges[i].e_j);

		if (grid_edges[i].type == 0){
			index_1.i = grid_edges[i].s_i;
			index_1.j = grid_edges[i].s_j - 1;
			index_2.i = grid_edges[i].e_i;
			index_2.j = grid_edges[i].e_j - 1;
			index_3.i = grid_edges[i].s_i;
			index_3.j = grid_edges[i].s_j + 1;
			index_4.i = grid_edges[i].e_i;
			index_4.j = grid_edges[i].e_j + 1;
		}

		if (grid_edges[i].type == 1){
			index_1.i = grid_edges[i].s_i - 1;
			index_1.j = grid_edges[i].s_j;
			index_2.i = grid_edges[i].e_i - 1;
			index_2.j = grid_edges[i].e_j;
			index_3.i = grid_edges[i].s_i + 1;
			index_3.j = grid_edges[i].s_j;
			index_4.i = grid_edges[i].e_i + 1;
			index_4.j = grid_edges[i].e_j;
		}

		int lable_0 = GetIndex(grid_edges, index_1, index_2);
		int lable_1 = GetIndex(grid_edges, index_1, index_s);
		int lable_2 = GetIndex(grid_edges, index_2, index_e);
		int lable_3 = GetIndex(grid_edges, index_3, index_4);
		int lable_4 = GetIndex(grid_edges, index_3, index_s);
		int lable_5 = GetIndex(grid_edges, index_4, index_e);

		GridEdgeRelation gr;
		if (lable_1 >= 0 && (grid[index_s.i][index_s.j] == 0 || (grid[index_s.i][index_s.j] == 1 && grid[index_2.i][index_2.j] == 0))){
			gr.ids.push_back(lable_1);
		}
		if (lable_2 >= 0 && (grid[index_e.i][index_e.j] == 0 || (grid[index_e.i][index_e.j] == 1 && grid[index_1.i][index_1.j] == 0))){
			gr.ids.push_back(lable_2);
		}
		if (lable_4 >= 0 && (grid[index_s.i][index_s.j] == 0 || (grid[index_s.i][index_s.j] == 1 && grid[index_4.i][index_4.j] == 0))){
			gr.ids.push_back(lable_4);
		}
		if (lable_5 >= 0 && (grid[index_e.i][index_e.j] == 0 || (grid[index_e.i][index_e.j] == 1 && grid[index_3.i][index_3.j] == 0))){
			gr.ids.push_back(lable_5);
		}
		if (lable_0 >= 0 && lable_1<0 && lable_2<0){
			gr.ids.push_back(lable_0);
		}
		if (lable_3 >= 0 && lable_4<0 && lable_5<0){
			gr.ids.push_back(lable_3);
		}
		grid_relations.push_back(gr);
		//if ((grid_edges[i].s_i == 11+1 && grid_edges[i].s_j == 15+1) ||
		//	(grid_edges[i].e_i == 11+1 && grid_edges[i].e_j == 15+1))
		//{
		//	std::cout << "Edge: (" << grid_edges[i].s_i - 1 << ", " << grid_edges[i].s_j - 1 << ")_("
		//		<< grid_edges[i].e_i - 1 << ", " << grid_edges[i].e_j - 1 << ")" << std::endl;
		//	for (int j = 0; j < gr.ids.size(); j++)
		//	{
		//		int index = gr.ids[j];
		//		std::cout << "(" << grid_edges[index].s_i - 1 << ", " << grid_edges[index].s_j - 1 << ")_("
		//			<< grid_edges[index].e_i - 1 << ", " << grid_edges[index].e_j - 1 << ")" << std::endl;
		//	}
		//}
	}

	//std::vector<GridEdgeRelation> grid_relations;
	//std::vector<bool> grid_edges_used;
	//for (int i = 0; i < grid_edges.size(); i++)
	std::vector<std::vector<int>> grid_boundaries;
	while (true){
		int start_edge_index = -1;

		for (int i = 0; i < grid_edges_used.size(); i++){
			if (!grid_edges_used[i]){
				start_edge_index = i;
				break;
			}
		}
		if (start_edge_index < 0) break;
		grid_edges_used[start_edge_index] = true;
		std::vector<int> one_boundary(1, start_edge_index);
		while (true){
			int next_edge_index = -1;
			for (int i = 0; i < grid_relations[start_edge_index].ids.size(); i++){
				if (!grid_edges_used[grid_relations[start_edge_index].ids[i]]){
					next_edge_index = grid_relations[start_edge_index].ids[i];
					break;
				}
			}
			if (next_edge_index < 0) break;
			grid_edges_used[next_edge_index] = true;
			one_boundary.push_back(next_edge_index);
			start_edge_index = next_edge_index;
		}
		grid_boundaries.push_back(one_boundary);
	}

	//transform back

	#if 0
	{
		for (int i = 0; i < grid_boundaries.size(); i++){
			std::vector<double> xs;
			std::vector<double> ys;
			for (int j = 0; j < grid_boundaries[i].size(); j++){
				GridIndex index_s(grid_edges[grid_boundaries[i][j]].s_i, grid_edges[grid_boundaries[i][j]].s_j);
				GridIndex index_e(grid_edges[grid_boundaries[i][j]].e_i, grid_edges[grid_boundaries[i][j]].e_j);
				xs.push_back((index_s.i + index_e.i) / 2.0 - 1.0);
				ys.push_back((index_s.j + index_e.j) / 2.0 - 1.0);
			}
			boundary_xs.push_back(xs);
			boundary_ys.push_back(ys);
		}
	}
	#else
	{
		for (int i = 0; i < grid_boundaries.size(); i++){
			std::vector<GridIndex> one_boundary;
			for (int j = 0; j < grid_boundaries[i].size(); j++){
				GridIndex index_s(grid_edges[grid_boundaries[i][j]].s_i, grid_edges[grid_boundaries[i][j]].s_j);
				GridIndex index_e(grid_edges[grid_boundaries[i][j]].e_i, grid_edges[grid_boundaries[i][j]].e_j);
				if (grid[index_s.i][index_s.j] == 1){
					one_boundary.push_back(index_s);
				}
				else{
					one_boundary.push_back(index_e);
				}
			}
			Decomposition_Mapping(one_boundary, boundary_xs, boundary_ys,true);
			std::vector<GridIndex>().swap(one_boundary);
		}
	}
	#endif

	std::vector<std::vector<int>>().swap(grid);
	std::vector<GridEdge>().swap(grid_edges);
	std::vector<GridEdgeRelation>().swap(grid_relations);
	std::vector<bool>().swap(grid_edges_used);
	std::vector<std::vector<int>>().swap(grid_boundaries);
}


//there are bugs in this function
void CGAL_Image_Grid_Decomposition_Conservative(std::vector<std::vector<int>> &image, std::vector<std::vector<double>> &boundary_xs, std::vector<std::vector<double>> &boundary_ys)
{
	std::vector<std::vector<Index>> boundaries;

	//#pragma omp parallel 
	{
		//2
		for (int i = 0; i < image.size(); i++)
		{
			for (int j = 0; j < image[i].size(); j++)
			{
				if (image[i][j] == 1)
				{
					bool run = true;
					for (int k = -1; k < 2 && run; k++)
					{
						for (int m = -1; m < 2 && run; m++)
						{
							int index_0 = i + k;
							int index_1 = j + m;
							if (index_0 < 0 || index_1 < 0 || index_0 >= image.size() || index_1 >= image[i].size() || image[index_0][index_1] == 0)
							{
								image[i][j] = 2;
								run = false;
							}
						}
					}
				}
			}
		}

		//3
		std::vector<Index> lables_3;
		std::vector<bool> lables_3_used;
		for (int i = 0; i < image.size(); i++)
		{
			for (int j = 0; j < image[i].size(); j++)
			{
				if (image[i][j] == 2)
				{
					bool run = true;
					for (int k = -1; k < 2 && run; k++)
					{
						for (int m = -1; m < 2 && run; m++)
						{
							int index_0 = i + k;
							int index_1 = j + m;
							if (!(index_0 < 0 || index_1 < 0 || index_0 >= image.size() || index_1 >= image[i].size()) && image[index_0][index_1] == 1)
							{
								image[i][j] = 3 + lables_3.size();
								lables_3.push_back(Index(i, j));
								lables_3_used.push_back(true);
								run = false;
							}
						}
					}

				}
			}
		}

		//order
		while (true)
		{
			//start_index
			Index start_index(-1, -1);
			for (int i = 0; i < lables_3_used.size(); i++)
			{
				if (lables_3_used[i])
				{
					start_index.x = lables_3[i].x;
					start_index.y = lables_3[i].y;
					break;
				}
			}

			if (start_index.x < 0 || start_index.y < 0) break;

			std::vector<Index> one_boundary;
			one_boundary.push_back(start_index);
			lables_3_used[image[start_index.x][start_index.y] - 3] = false;
			image[start_index.x][start_index.y] = 2;


			while (true)
			{
				Index next_index(-1, -1);

				//next_index
				bool run = true;
				for (int k = -1; k < 2 && run; k++)
				{
					for (int m = -1; m < 2 && run; m++)
					{
						if (abs(k) + abs(m) == 1)
						{
							int index_0 = start_index.x + k;
							int index_1 = start_index.y + m;
							if (!(index_0 < 0 || index_1 < 0 || index_0 >= image.size() || index_1 >= image[0].size()))
							{
								if (image[index_0][index_1] >= 3)
								{
									next_index.x = index_0;
									next_index.y = index_1;
									run = false;
								}
							}
						}

					}
				}

				if (next_index.x < 0 || next_index.y < 0) break;

				one_boundary.push_back(next_index);
				lables_3_used[image[next_index.x][next_index.y] - 3] = false;
				image[next_index.x][next_index.y] = 2;

				start_index.x = next_index.x;
				start_index.y = next_index.y;
			}

			if (abs(one_boundary[0].x - one_boundary[one_boundary.size() - 1].x) <= 1 && abs(one_boundary[0].y - one_boundary[one_boundary.size() - 1].y) <= 1 && one_boundary.size() >= 3)
				boundaries.push_back(one_boundary);
		}
	}

	for (int i = 0; i < boundaries.size(); i++)
	{
		std::vector<double> xs;
		std::vector<double> ys;
		for (int j = 0; j < boundaries[i].size(); j++)
		{
			xs.push_back(boundaries[i][j].x);
			ys.push_back(boundaries[i][j].y);
		}

		boundary_xs.push_back(xs);
		boundary_ys.push_back(ys);
	}

	std::vector<std::vector<Index>>().swap(boundaries);
}



void CGAL_Image_Grid_Decomposition(std::vector<std::vector<int>> &image, Vector2d2 &boundaries)
{
	std::vector<std::vector<double>> xs;
	std::vector<std::vector<double>> ys;

	CGAL_Image_Grid_Decomposition(image, xs, ys);

	for (int i = 0; i < xs.size(); i++)
	{
		Vector2d1 one_boundary;
		for (int j = 0; j < xs[i].size(); j++)
		{
			one_boundary.push_back(Vector2d(xs[i][j], ys[i][j]));
		}
		boundaries.push_back(one_boundary);
	}
}

void CGAL_Image_Grid_Decomposition_Conservative(std::vector<std::vector<int>> &image, Vector2d2 &boundaries)
{
	std::vector<std::vector<double>> xs;
	std::vector<std::vector<double>> ys;

	CGAL_Image_Grid_Decomposition_Conservative(image, xs, ys);

	for (int i = 0; i < xs.size(); i++)
	{
		Vector2d1 one_boundary;
		for (int j = 0; j < xs[i].size(); j++)
		{
			one_boundary.push_back(Vector2d(xs[i][j], ys[i][j]));
		}
		boundaries.push_back(one_boundary);
	}
}

/***************************************************************************************************/
void CGAL_Surface_Decomposition(std::string path, std::vector<double> &face_sdf, int &regions_nb, std::vector<int> &face_segments)
{
	Polyhedron_3 polyhedron;

	Construct_Polyhedron(polyhedron, path);

	// create a property-map for SDF values
	typedef std::map<Polyhedron_3::Facet_const_handle, double> Facet_double_map;
	Facet_double_map internal_sdf_map;
	boost::associative_property_map<Facet_double_map> sdf_property_map(internal_sdf_map);

	// compute SDF values using default parameters for number of rays, and cone angle
	//CGAL::sdf_values(polyhedron, sdf_property_map);

	int index = 0;
	for (Polyhedron_3::Facet_const_iterator facet_it = polyhedron.facets_begin();
		facet_it != polyhedron.facets_end(); ++facet_it) {
		sdf_property_map[facet_it] = face_sdf[index];
		index++;
	}

	// create a property-map for segment-ids
	typedef std::map<Polyhedron_3::Facet_const_handle, std::size_t> Facet_int_map;
	Facet_int_map internal_segment_map;
	boost::associative_property_map<Facet_int_map> segment_property_map(internal_segment_map);

	// segment the mesh using default parameters for number of levels, and smoothing lambda
	// Any other scalar values can be used instead of using SDF values computed using the CGAL function
	//std::size_t number_of_segments = CGAL::segmentation_from_sdf_values(polyhedron, sdf_property_map, segment_property_map);
	
	//regions_nb = number_of_segments;

	const std::size_t number_of_clusters = 100;       // use 4 clusters in soft clustering
	const double smoothing_lambda = 0.08;  // importance of surface features, suggested to be in-between [0,1]
	regions_nb=CGAL::segmentation_from_sdf_values(
		polyhedron, sdf_property_map, segment_property_map, number_of_clusters, smoothing_lambda);

	for (Polyhedron_3::Facet_const_iterator facet_it = polyhedron.facets_begin();
		facet_it != polyhedron.facets_end(); ++facet_it) {
		face_segments.push_back(segment_property_map[facet_it]);
	}
}

bool InterpolationPoint(Vector3d v0, Vector3d v1, double psd0, double  psd1, double d, Vector3d &inter)
{
	if (Math::IsAlmostZero(d - psd0)) psd0 = d;
	if (Math::IsAlmostZero(psd1 - d)) psd1 = d;

	if (!Math::IsAlmostZero(psd1 - psd0) && ((d >= psd0&&d <= psd1) || (d >= psd1&&d <= psd0)))
	{
		if (d >= psd1&&d <= psd0)
		{
			Vector3d v = v0;
			v0 = v1;
			v1 = v;

			double d = psd0;
			psd0 = psd1;
			psd1 = d;
		}
		inter = (float)((d - psd0) / (psd1 - psd0))*v1 + (float)((psd1 - d) / (psd1 - psd0))*v0;
		return true;
	}
	else
	{
		return false;
	}
}

void CGAL_3D_Connecting_Segments(Vector2d2 &segments, Vector2d2 &lines)
{
	//save connecting relations
	std::vector<bool> used(segments.size(), false);

	std::vector<int> relations;
#pragma region get_relations
	for (int i = 0; i < segments.size(); i++)
	{
		for (int j = 0; j < segments.size(); j++)
		{
			if (i != j&&!used[i] && !used[j])
			{
				bool b_0_0 = Math::IsAlmostZero_Double(CGAL_2D_Distance_Point_Point(segments[i][0], segments[j][0]), 1.0E-09);
				bool b_0_1 = Math::IsAlmostZero_Double(CGAL_2D_Distance_Point_Point(segments[i][0], segments[j][1]), 1.0E-09);
				bool b_1_0 = Math::IsAlmostZero_Double(CGAL_2D_Distance_Point_Point(segments[i][1], segments[j][0]), 1.0E-09);
				bool b_1_1 = Math::IsAlmostZero_Double(CGAL_2D_Distance_Point_Point(segments[i][1], segments[j][1]), 1.0E-09);

				if ((b_0_0&&b_1_1) || (b_0_1&&b_1_0))
				{
					used[j] = true;
					continue;
				}

				if (b_0_0)
				{
					relations.push_back(i);
					relations.push_back(0);
					relations.push_back(j);
					relations.push_back(0);
					continue;
				}
				if (b_0_1)
				{
					relations.push_back(i);
					relations.push_back(0);
					relations.push_back(j);
					relations.push_back(1);
					continue;
				}
				if (b_1_0)
				{
					relations.push_back(i);
					relations.push_back(1);
					relations.push_back(j);
					relations.push_back(0);
					continue;
				}
				if (b_1_1)
				{
					relations.push_back(i);
					relations.push_back(1);
					relations.push_back(j);
					relations.push_back(1);
					continue;
				}
			}
		}
	}
#pragma endregion

	std::vector<std::vector<int>> ones;


	while (true)
	{
		int index = -1;
		int end = -1;

		for (int i = 0; i < segments.size(); i++)
		{
			if (!used[i]){
				index = i;
				end = 0;
				used[i] = true;
				break;
			}
		}

		if (index < 0)break;

		Vector2d1 line(1, segments[index][end]);

		std::vector<int> one(1, index);

		while (true)
		{
			end = 1 - end;
			bool search = false;
			for (int i = 0; i < relations.size(); i = i + 4)
			{
				if (relations[i] == index&&relations[i + 1] == end && !used[relations[i + 2]])
				{
					line.push_back(segments[relations[i + 2]][relations[i + 3]]);
					one.push_back(relations[i + 2]);
					index = relations[i + 2];
					end = relations[i + 3];
					used[index] = true;
					search = true;
					break;
				}
				if (relations[i + 2] == index&&relations[i + 3] == end && !used[relations[i]])
				{
					line.push_back(segments[relations[i]][relations[i + 1]]);
					one.push_back(relations[i]);
					index = relations[i];
					end = relations[i + 1];
					used[index] = true;
					search = true;
					break;
				}
			}
			if (!search){ break; }
		}

		ones.push_back(one);
		lines.push_back(line);
	}
}

void CGAL_3D_Connecting_Segments(Vector3d2 &segments, Vector3d2 &lines)
{
	//save connecting relations
	std::vector<bool> used(segments.size(), false);

	std::vector<int> relations;
#pragma region get_relations
	for (int i = 0; i < segments.size(); i++)
	{
		for (int j = 0; j < segments.size(); j++)
		{
			if (i != j&&!used[i] && !used[j])
			{
				bool b_0_0 = Math::IsAlmostZero_Double(CGAL_3D_Distance_Point_Point(segments[i][0], segments[j][0]), 1.0E-09);
				bool b_0_1 = Math::IsAlmostZero_Double(CGAL_3D_Distance_Point_Point(segments[i][0], segments[j][1]), 1.0E-09);
				bool b_1_0 = Math::IsAlmostZero_Double(CGAL_3D_Distance_Point_Point(segments[i][1], segments[j][0]), 1.0E-09);
				bool b_1_1 = Math::IsAlmostZero_Double(CGAL_3D_Distance_Point_Point(segments[i][1], segments[j][1]), 1.0E-09);

				if ((b_0_0&&b_1_1) || (b_0_1&&b_1_0))
				{
					used[j] = true;
					continue;
				}

				if (b_0_0)
				{
					relations.push_back(i);
					relations.push_back(0);
					relations.push_back(j);
					relations.push_back(0);
					continue;
				}
				if (b_0_1)
				{
					relations.push_back(i);
					relations.push_back(0);
					relations.push_back(j);
					relations.push_back(1);
					continue;
				}
				if (b_1_0)
				{
					relations.push_back(i);
					relations.push_back(1);
					relations.push_back(j);
					relations.push_back(0);
					continue;
				}
				if (b_1_1)
				{
					relations.push_back(i);
					relations.push_back(1);
					relations.push_back(j);
					relations.push_back(1);
					continue;
				}
			}
		}
	}
#pragma endregion

	std::vector<std::vector<int>> ones;


	while (true)
	{
		int index = -1;
		int end = -1;

		for (int i = 0; i < segments.size(); i++)
		{
			if (!used[i]){
				index = i;
				end = 0;
				used[i] = true;
				break;
			}
		}

		if (index < 0)break;

		Vector3d1 line(1, segments[index][end]);

		std::vector<int> one(1, index);

		while (true)
		{
			end = 1 - end;
			bool search = false;
			for (int i = 0; i < relations.size(); i = i + 4)
			{
				if (relations[i] == index&&relations[i + 1] == end && !used[relations[i + 2]])
				{
					line.push_back(segments[relations[i + 2]][relations[i + 3]]);
					one.push_back(relations[i + 2]);
					index = relations[i + 2];
					end = relations[i + 3];
					used[index] = true;
					search = true;
					break;
				}
				if (relations[i + 2] == index&&relations[i + 3] == end && !used[relations[i]])
				{
					line.push_back(segments[relations[i]][relations[i + 1]]);
					one.push_back(relations[i]);
					index = relations[i];
					end = relations[i + 1];
					used[index] = true;
					search = true;
					break;
				}
			}
			if (!search){ break; }
		}

		ones.push_back(one);
		lines.push_back(line);
	}
}

bool CheckingEndEdge(Index &edge_0, Index &edge_1)
{
	return (edge_0.x == edge_1.x&&edge_0.y == edge_1.y)||
		   (edge_0.y == edge_1.x&&edge_0.x == edge_1.y);
}

void CGAL_3D_Connecting_Segments(Vector3d2 &segments, std::vector<std::vector<Index>> &edges, Vector3d2 &lines)
{
	//save connecting relations
	std::vector<bool> used(segments.size(), false);

	std::vector<int> relations;
#pragma region get_relations
	for (int i = 0; i < segments.size(); i++)
	{
		for (int j = 0; j < segments.size(); j++)
		{
			if (i != j&&!used[i] && !used[j])
			{
				bool b_0_0 = CheckingEndEdge(edges[i][0], edges[j][0]);
				bool b_0_1 = CheckingEndEdge(edges[i][0], edges[j][1]);
				bool b_1_0 = CheckingEndEdge(edges[i][1], edges[j][0]);
				bool b_1_1 = CheckingEndEdge(edges[i][1], edges[j][1]);

				//bool b_0_0 = Math::IsAlmostZero(CGAL_3D_Distance_Point_Point(segments[i][0], segments[j][0]));
				//bool b_0_1 = Math::IsAlmostZero(CGAL_3D_Distance_Point_Point(segments[i][0], segments[j][1]));
				//bool b_1_0 = Math::IsAlmostZero(CGAL_3D_Distance_Point_Point(segments[i][1], segments[j][0]));
				//bool b_1_1 = Math::IsAlmostZero(CGAL_3D_Distance_Point_Point(segments[i][1], segments[j][1]));

				if ((b_0_0&&b_1_1) || (b_0_1&&b_1_0))
				{
					used[j] = true;
					continue;
				}

				if (b_0_0)
				{
					relations.push_back(i);
					relations.push_back(0);
					relations.push_back(j);
					relations.push_back(0);
					continue;
				}
				if (b_0_1)
				{
					relations.push_back(i);
					relations.push_back(0);
					relations.push_back(j);
					relations.push_back(1);
					continue;
				}
				if (b_1_0)
				{
					relations.push_back(i);
					relations.push_back(1);
					relations.push_back(j);
					relations.push_back(0);
					continue;
				}
				if (b_1_1)
				{
					relations.push_back(i);
					relations.push_back(1);
					relations.push_back(j);
					relations.push_back(1);
					continue;
				}
			}
		}
	}
#pragma endregion




	while (true)
	{
		int index = -1;
		int end = -1;

		for (int i = 0; i < segments.size(); i++)
		{
			if (!used[i]){
				index = i;
				end = 0;
				used[i] = true;
				break;
			}
		}

		if (index < 0)break;

		Vector3d1 line(1, segments[index][end]);

	

		while (true)
		{
			end = 1 - end;
			bool search = false;
			for (int i = 0; i < relations.size(); i = i + 4)
			{
				if (relations[i] == index&&relations[i + 1] == end && !used[relations[i + 2]])
				{
					line.push_back(segments[relations[i + 2]][relations[i + 3]]);
					index = relations[i + 2];
					end = relations[i + 3];
					used[index] = true;
					search = true;
					break;
				}
				if (relations[i + 2] == index&&relations[i + 3] == end && !used[relations[i]])
				{
					line.push_back(segments[relations[i]][relations[i + 1]]);
					index = relations[i];
					end = relations[i + 1];
					used[index] = true;
					search = true;
					break;
				}
			}
			if (!search){ break; }
		}

		lines.push_back(line);
	}
}


/***************************************************************************************************/
bool CGAL_3D_Mesh_Extract_Isoline(Vector3d1 &vecs, std::vector<int>& face_id_0, std::vector<int>& face_id_1, std::vector<int>& face_id_2, std::vector<double> &psd,
	double d, Vector3d2 &isolines)
{
	//compute gredient
	Vector3d1 vecs_gradients;
	Vector3d1 faces_gradients;
	CGAL_3D_Mesh_Gradient(vecs, face_id_0, face_id_1, face_id_2, psd, vecs_gradients, faces_gradients);

	Vector3d2 segments;
	std::vector<std::vector<Index>> edges;

	for (int i = 0; i < face_id_0.size(); i++)
	{
		int index_0 = face_id_0[i];
		int index_1 = face_id_1[i];
		int index_2 = face_id_2[i];

		Vector3d v_0 = vecs[index_0];
		Vector3d v_1 = vecs[index_1];
		Vector3d v_2 = vecs[index_2];

		Vector3d  n = Math::GetCrossproduct(v_0 - v_2, v_1 - v_2);
		Vector3d gradient = faces_gradients[i];

		Vector3d inter_0_1;
		Vector3d inter_1_2;
		Vector3d inter_2_0;

		bool b_0 = InterpolationPoint(v_0, v_1, psd[index_0], psd[index_1], d, inter_0_1);
		bool b_1 = InterpolationPoint(v_1, v_2, psd[index_1], psd[index_2], d, inter_1_2);
		bool b_2 = InterpolationPoint(v_2, v_0, psd[index_2], psd[index_0], d, inter_2_0);

		if (b_0&&b_1)
		{
			Vector3d1 segment;
			segment.push_back(inter_0_1);
			segment.push_back(inter_1_2);

			std::vector<Index> endedge;
			endedge.push_back(Index(index_0, index_1));
			endedge.push_back(Index(index_1, index_2));

			if (!Math::IsAlmostZero(CGAL_3D_Distance_Point_Point(segment[0], segment[1])))
			{
				Vector3d  gradient_n = Math::GetCrossproduct(segment[0] - segment[1], gradient);
				double angle = Math::GetAngleBetween(gradient_n, n);
				if (angle > Math::Math_PI / 2.0)
					std::reverse(segment.begin(),segment.end());

				segments.push_back(segment);
				edges.push_back(endedge);
			}
		}
		if (b_1&&b_2)
		{
			Vector3d1 segment;
			segment.push_back(inter_1_2);
			segment.push_back(inter_2_0);

			std::vector<Index> endedge;
			endedge.push_back(Index(index_1, index_2));
			endedge.push_back(Index(index_2, index_0));

			if (!Math::IsAlmostZero(CGAL_3D_Distance_Point_Point(segment[0], segment[1])))
			{
				Vector3d  gradient_n = Math::GetCrossproduct(segment[0] - segment[1], gradient);
				double angle = Math::GetAngleBetween(gradient_n, n);
				if (angle > Math::Math_PI / 2.0)
					std::reverse(segment.begin(), segment.end());

				segments.push_back(segment);
				edges.push_back(endedge);
			}
		}
		if (b_2&&b_0)
		{
			Vector3d1 segment;
			segment.push_back(inter_2_0);
			segment.push_back(inter_0_1);

			std::vector<Index> endedge;
			endedge.push_back(Index(index_2, index_0));
			endedge.push_back(Index(index_0, index_1));

			if (!Math::IsAlmostZero(CGAL_3D_Distance_Point_Point(segment[0], segment[1])))
			{
				Vector3d  gradient_n = Math::GetCrossproduct(segment[0] - segment[1], gradient);
				double angle = Math::GetAngleBetween(gradient_n, n);
				if (angle > Math::Math_PI / 2.0)
					std::reverse(segment.begin(), segment.end());

				segments.push_back(segment);
				edges.push_back(endedge);
			}
		}
	}

	std::cout << "CGAL_3D_Mesh_Extract_Isoline: " << segments.size() << std::endl;


	//std::ofstream export_fie("G:\\segments.obj");
	//int export_int = 1;

	//for (int i = 0; i < segments.size(); i++)
	//{
	//	CGAL_Export_Path_Segment(export_fie, export_int, "segment_" + Math::IntString(i), 1.0, 0.0, 0.0, segments[i][0], segments[i][1], 0.002);
	//}

	//export_fie.clear();
	//export_fie.close();


	if (segments.size() == 0) return false;

	CGAL_3D_Connecting_Segments(segments, isolines);

	//std::ofstream export_fie_0("G:\\isolines.obj");
	//export_int = 1;

	//for (int i = 0; i < isolines.size(); i++)
	//{
	//	for (int j = 0; j < isolines.size()-1; j++)
	//	{
	//		CGAL_Export_Path_Segment(export_fie_0, export_int, "isolines_" + Math::IntString(i) + "_" + Math::IntString(j), 1.0, 0.0, 0.0, isolines[i][j], isolines[i][(j + 1)], 0.002);
	//	}

	//}

	//export_fie_0.clear();
	//export_fie_0.close();


	return true;
}
/***************************************************************************************************/

//shortest geodesic distance
/***************************************************************************************************/

Vector3d NearestPoint(Polyhedron_3 &polyhedron, Tree &tree,Vector3d &source)
{
	Poly_point_3 query(source[0], source[1], source[2]);
	Point_and_primitive_id pp = tree.closest_point_and_primitive(query);
	return Vector3d(pp.first.x(), pp.first.y(), pp.first.z());
}

Vector3d1 CGAL_Project_Points_Onto_Surface(Vector3d1& vecs, std::vector<int>& face_id_0, std::vector<int>& face_id_1, std::vector<int>& face_id_2, Vector3d1& points)
{
	//construct polyhedron
	Polyhedron_3 polyhedron;
	Construct_Polyhedron(polyhedron, vecs, face_id_0, face_id_1, face_id_2);

	Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);
	tree.accelerate_distance_queries();

	Vector3d1 temp ;
	for (int i = 0; i < points.size(); i++)
	{
		temp.push_back(NearestPoint(polyhedron, tree, points[i]));
	}
	return temp;
}

Vector3d1 CGAL_Project_Points_Onto_Surface(std::string path, Vector3d1& points)
{
	Polyhedron_3 polyhedron;
	Construct_Polyhedron(polyhedron,path);

	Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);
	tree.accelerate_distance_queries();

	Vector3d1 temp;
	for (int i = 0; i < points.size(); i++)
	{
		temp.push_back(NearestPoint(polyhedron, tree, points[i]));
	}
	return temp;
}

void CGAL_BSplineCurveFit(Vector3d1& samples, Vector3d1& output)
{
	/********************************/
	Vector3d1().swap(samples);
	//float mult = 2.0f / (1000 - 1), t;
	//int i;
	//for (i = 0; i < 1000; ++i)
	//{
	//	t = -1.0f + mult*i;
	//	float angle = 2.0f*Wm5::Mathf::TWO_PI*t;
	//	float amplitude = 1.0f - t*t;
	//	Vector3d v(amplitude*Wm5::Mathf::Cos(angle), amplitude*Wm5::Mathf::Sin(angle),t);
	//	samples.push_back(v);
	//}
	/********************************/

	int dimension = 3;
	int numSamples = samples.size();

	Wm5::Vector3d* mSamples;
	mSamples = new1<Wm5::Vector3d>(numSamples);
	for (int i = 0; i < numSamples; ++i)
	{
		mSamples[i].X() = samples[i][0];
		mSamples[i].Y() = samples[i][1];
		mSamples[i].Z() = samples[i][2];
	}

	int  mNumCtrlPoints = numSamples/ 2;
	int mDegree = 6;

	Wm5::BSplineCurveFitd* mSpline;

	mSpline = new0 Wm5::BSplineCurveFitd(dimension, numSamples, (const double*)mSamples, mDegree, mNumCtrlPoints);

	// Sample it the same number of times as the original data.

	float mult = 1.0 / (numSamples*2 - 1);
	for (int i = 0; i < numSamples*2; ++i)
	{
		double *pos = new double[3];
		mSpline->GetPosition(mult*i, pos);
		output.push_back(Vector3d(pos[0],pos[1],pos[2]));
	}
}

//generate data
Vector2d1 GenerateConcavityData(int nb, int period, double max_r, double min_r)
{
	Vector2d1 cnc_path;
	for (int i = 0; i < nb; i++)
	{
		double angle = (double)i / (double)nb*Wm5::Mathf::PI*2.0;
		double r = (cos(angle*period) + 1.0) / 2.0*(max_r - min_r) + min_r;

		double x, y;
		x = r*cos(angle);
		y = r*sin(angle);
		cnc_path.push_back(Vector2d(x, y));
	}
	return cnc_path;
}



//
//void CGAL_Output_Boundary(std::string path, Vector3d1 &vecs, std::vector<std::vector<int>> &face_ids, std::vector<int> &triangles_lables)
//{
//	std::vector<int> face_id_0;
//	std::vector<int> face_id_1;
//	std::vector<int> face_id_2;
//
//	std::vector<int> lables_0(vecs.size(), -1);
//	std::vector<int> lables_1(vecs.size(), -1);
//	for (int i = 0; i < face_ids.size(); i++)
//	{
//		face_id_0.push_back(face_ids[i][0]);
//		face_id_1.push_back(face_ids[i][1]);
//		face_id_2.push_back(face_ids[i][2]);
//		if (triangles_lables[i] == 0)
//		{
//			lables_0[face_ids[i][0]]=0;
//			lables_0[face_ids[i][1]] = 0;
//			lables_0[face_ids[i][2]] = 0;
//		}
//		if (triangles_lables[i] == 1)
//		{
//			lables_1[face_ids[i][0]] = 0;
//			lables_1[face_ids[i][1]] = 0;
//			lables_1[face_ids[i][2]] = 0;
//		}
//	}
//
//	int nb = 0;
//
//	for (int i = 0; i < vecs.size(); i++)
//	{
//		if (lables_0[i] == 0 && lables_1[i] == 0)
//		{
//			nb++;
//		}
//	}
//
//	std::ofstream file("G:\\0_heat.source");
//	file << nb << std::endl;
//
//	for (int i = 0; i < vecs.size(); i++)
//	{
//		if (lables_0[i] == 0 && lables_1[i] == 0)
//		{
//			file << "1 " << i << std::endl;
//		}
//	}
//
//	file.clear();
//	file.close();
//}
//


double ComputeGapFromScallop(double surface_curvature, double R_cutter, double scallop)
{
	if (Math::IsAlmostZero(surface_curvature))
	{
		return 2.0*sqrt(2.0*scallop*R_cutter / (1.0 + R_cutter*surface_curvature));
	}
	else
	{
		if (surface_curvature <0.0&&-1.0 / surface_curvature < R_cutter)
			return  100000.0;
		else
			return 2.0*sqrt(2.0*scallop*R_cutter / (1.0 + R_cutter*surface_curvature));
	}
}

void Laplace_Mesh_Vertex_Values(Vector3d1 &vecs, std::vector<int> &face_id_0, std::vector<int> &face_id_1, std::vector<int> &face_id_2, std::vector<double> &values,int iterations)
{
	std::vector<std::vector<int>> neighs;
	CGAL_3D_Triangle_Mesh_Vecs_Neighbors(vecs, face_id_0, face_id_1, face_id_2, neighs);

	std::vector<double> iter_values(vecs.size(),0.0);

	for (int i = 0; i < iterations; i++)
	{
		for (int j = 0; j < vecs.size(); j++)
		{
			for (int k = 0; k < neighs[j].size(); k++)
			{
				iter_values[j] += values[neighs[j][k]];
			}
			iter_values[j] = iter_values[j]/(double)neighs[j].size();

			//iter_values[j] = 0.5*iter_values[j] + 0.5*values[i];
		}

		std::vector<double>().swap(values);
		values = iter_values;
		for (int j = 0; j < vecs.size(); j++)
			iter_values[j] = 0.0;
	}

}


void CGAL_3D_Mesh_Gradient(Vector3d1 &vecs, std::vector<int> &face_id_0, std::vector<int> &face_id_1, std::vector<int> &face_id_2, std::vector<double> &psd,
	Vector3d1 &vecs_gradients, Vector3d1 &faces_gradients)
{
	Vector3d1().swap(vecs_gradients);
	Vector3d1().swap(faces_gradients);

	for (int i = 0; i < vecs.size(); i++)
		vecs_gradients.push_back(Vector3d(0.0, 0.0, 0.0));

	std::vector<double> areas(vecs.size(), 0.0);

	for (int i = 0; i < face_id_0.size(); i++)
	{
		int index_0 = face_id_0[i];
		int index_1 = face_id_1[i];
		int index_2 = face_id_2[i];

		Vector3d v_0 = vecs[index_0];
		Vector3d v_1 = vecs[index_1];
		Vector3d v_2 = vecs[index_2];

		Vector3d   n = Math::GetCrossproduct(v_0 - v_2, v_1 - v_2);
		Math::SetVectorLength(n, 1.0);

		double area = CGAL_3D_One_Triangle_Area(v_0, v_1, v_2);

		Vector3d gradient(0.0, 0.0, 0.0);

		double d_0 = (float)psd[index_0];
		double d_1 = (float)psd[index_1];
		double d_2 = (float)psd[index_2];

		gradient += (float)psd[index_0] * Math::GetCrossproduct(n, v_2 - v_1);
		gradient += (float)psd[index_1] * Math::GetCrossproduct(n, v_0 - v_2);
		gradient += (float)psd[index_2] * Math::GetCrossproduct(n, v_1 - v_0);

		Vector3d face_gradient = gradient / (float)(2.0 * area);

		faces_gradients.push_back(face_gradient);

		vecs_gradients[index_0] += (float)area*face_gradient;
		areas[index_0] += area;

		vecs_gradients[index_1] += (float)area*face_gradient;
		areas[index_1] += area;

		vecs_gradients[index_2] += (float)area*face_gradient;
		areas[index_2] += area;

	}

	for (int i = 0; i < vecs.size(); i++)
	{
		vecs_gradients[i] = vecs_gradients[i] / (float)areas[i];
	}
}

//Save only four decimals for input point "p"
Point_3 FourDecimals(Point_3 p)
{
	double d_x = floor(p.x() * 10000.000f + 0.5) / 10000.000f;
	double d_y = floor(p.y() * 10000.000f + 0.5) / 10000.000f;
	double d_z = floor(p.z() * 10000.000f + 0.5) / 10000.000f;
	return Point_3(d_x,d_y,d_z);
}

//Project p onto the planar surface of 3d triangle
//Checking the position relationship between the p and 3d triangle
//face: 3d triangle
//p: 3d point
//return true: inside
//return false: outside
bool OutsidePointInsideTriangle(Poly_facet_iterator &face, Vector3d p)
{
	bool bbbb0;
	bool bbbb1;

	if (false)
	{
		Point_3 p0 = face->halfedge()->next()->next()->vertex()->point();
		Point_3 p1 = face->halfedge()->vertex()->point();
		Point_3 p2 = face->halfedge()->next()->vertex()->point();
		
		Plane_3 plane(p1, CGAL::cross_product(p2 - p1, p0 - p1));
		Vector2d p_2d = PointVector2d(point_to_2d(VectorPoint3d(p),plane));
		Vector2d p0_2d = PointVector2d(point_to_2d(p0, plane));
		Vector2d p1_2d = PointVector2d(point_to_2d(p1, plane));
		Vector2d p2_2d = PointVector2d(point_to_2d(p2, plane));

		Vector2d1 py;
		py.push_back(p0_2d);
		py.push_back(p1_2d);
		py.push_back(p2_2d);

		bbbb0= CGAL_2D_Location_Point_Polygon(p_2d, py);

		//return bbbb0;
	}

	if (true)
	{
		Point_3 p0 = face->halfedge()->next()->next()->vertex()->point();
		Point_3 p1 = face->halfedge()->vertex()->point();
		Point_3 p2 = face->halfedge()->next()->vertex()->point();
		Plane_3 plane(p1, CGAL::cross_product(p2 - p1, p0 - p1));
		Point_3 project = plane.projection(VectorPoint3d(p));

		Vector3d v0 = PointVector3d(p0);
		Vector3d v1 = PointVector3d(p1);
		Vector3d v2 = PointVector3d(p2);
		Vector3d vp = PointVector3d(project);

		double u, v, w;
		CGAL_Barycentric(vp, v0, v1, v2, u, v, w);

		if ((u >= 0.0&&u <= 1.0) && (v >= 0.0&&v <= 1.0) && (w >= 0.0&&w <= 1.0))
		{
			bbbb1 = true;
		}
		else
		{
			bbbb1 = false;
		}

		if (false)
		if (bbbb1 != bbbb0)
		{
			//output the projecting points
			if (DebugInformation())
			{
				std::ofstream  export_file_output("Z:\\Documents\\Windows\\SmartSFC\\workspace\\CFS\\debugdebug.obj");
				int export_index_0 = 1;
				CGAL_Export_Path_Segment(export_file_output, export_index_0, "p", 1.0, 0.0, 0.0, v0, p, 0.1);
				CGAL_Export_Path_Segment(export_file_output, export_index_0, "vp", 1.0, 0.0, 0.0, v0, vp, 0.1);
				CGAL_Export_Path_Segment(export_file_output, export_index_0, "debug", 1.0, 0.0, 0.0, v0, v1, 0.1);
				CGAL_Export_Path_Segment(export_file_output, export_index_0, "debug", 1.0, 0.0, 0.0, v1, v2, 0.1);
				CGAL_Export_Path_Segment(export_file_output, export_index_0, "debug", 1.0, 0.0, 0.0, v2, v0, 0.1);

				export_file_output.clear();
				export_file_output.close();
			}

			int dsads = 10;
		}
		return bbbb1;
	}


}




void OneIterationSmoothBoundary(Vector3d1 &vecs, std::vector<std::vector<int>> &edges, Vector3d1 &intersections)
{
	for (int i = 0; i < intersections.size(); i++)
	{
		Vector3d pre_scp = intersections[(i - 1 + intersections.size()) % intersections.size()];
		Vector3d next_scp = intersections[(i + 1 + intersections.size()) % intersections.size()];

		Vector3d scp_s = vecs[edges[i][0]];
		Vector3d scp_e = vecs[edges[i][1]];

		Vector3d pre_map = CGAL_3D_Projection_Point_Segment(pre_scp, scp_s, scp_e);
		Vector3d next_map = CGAL_3D_Projection_Point_Segment(next_scp, scp_s, scp_e);

		double d_pre_map = CGAL_3D_Distance_Point_Point(pre_scp, pre_map);
		double d_next_map = CGAL_3D_Distance_Point_Point(next_scp, next_map);

		if (Math::IsAlmostZero(d_pre_map + d_next_map))
			intersections[i] = (float)0.5*pre_map + (float)0.5*next_map;
		else
		{
			double w_pre = d_next_map / (d_pre_map + d_next_map);
			double w_next = d_pre_map / (d_pre_map + d_next_map);
			intersections[i] = (float)w_pre*pre_map + (float)w_next*next_map;
		}
	}
}



void ComputeRemeshTriangles(const Vector3d1 &vecs,const std::vector<int> &face_id_0, const std::vector<int> &face_id_1, const std::vector<int> &face_id_2,
	const std::vector<std::pair<int,int>> &edges, const Vector3d1 &cutting_points,const std::vector<int> &multi_cps, const std::vector<int> &two_class_lables, const std::string output_path)
{
	auto HandleBoundaryTriangle = [](const Vector3d1 &vecs, const std::vector<int> &remesh_lables, const std::vector<std::pair<int, int>> &edges, const Vector3d1 &cutting_points,
		std::vector<std::vector<int>> &remesh_triangles, std::vector<int> &remesh_triangles_lables, Vector3d1 &new_points, std::vector<int> &cutting_point_ids,
		const int index_0, const int index_1, const int index_2)
	{
		auto GetEdgeMiddlePoint = [](const int index_0, const int index_1, const std::vector<std::pair<int, int>> &edges)
		{
			for (int i = 0; i < edges.size(); i++)
				if ((edges[i].first == index_0&&edges[i].second == index_1) || (edges[i].second == index_0&&edges[i].first == index_1))
					return i;
			std::cerr << "GetEdgeMiddlePoint Error..." << std::endl;
			system("pause");
			return -1;
		};

		Vector3d v_0 = vecs[index_0];
		Vector3d v_1 = vecs[index_1];
		Vector3d v_2 = vecs[index_2];

		int lable_0 = remesh_lables[index_0];
		int lable_1 = remesh_lables[index_1];
		int lable_2 = remesh_lables[index_2];

		int cut_0_1 = GetEdgeMiddlePoint(index_0, index_1, edges);
		int cut_0_2 = GetEdgeMiddlePoint(index_0, index_2, edges);

		Vector3d v_0_1 = cutting_points[cut_0_1];
		Vector3d v_0_2 = cutting_points[cut_0_2];

		bool b_0 = CGAL_3D_Distance_Point_Point(v_0, v_0_1) <= 0.001;
		bool b_1 = CGAL_3D_Distance_Point_Point(v_1, v_0_1) <= 0.001;
		bool b_2 = CGAL_3D_Distance_Point_Point(v_0, v_0_2) <= 0.001;
		bool b_3 = CGAL_3D_Distance_Point_Point(v_2, v_0_2) <= 0.001;

		std::vector<int> triangle;
		triangle.push_back(index_0);
		triangle.push_back(index_1);
		triangle.push_back(index_2);

		//A: 0.0-0.0
		if (b_1&&b_3)
		{
			remesh_triangles.push_back(triangle);
			remesh_triangles_lables.push_back(lable_0);

			cutting_point_ids[cut_0_1] = index_1;
			cutting_point_ids[cut_0_2] = index_2;

			return;
		}

		//B: 0.0-0.5
		if (b_1&&!b_2&&!b_3)
		{
			new_points.push_back(v_0_2);
			int index_0_2 = new_points.size() - 1 + vecs.size();

			std::vector<int> triangle_0 = { index_0, index_1, index_0_2 };
			std::vector<int> triangle_1 = { index_2, index_0_2, index_1 };

			remesh_triangles.push_back(triangle_0);
			remesh_triangles_lables.push_back(lable_0);

			remesh_triangles.push_back(triangle_1);
			remesh_triangles_lables.push_back(lable_1);


			cutting_point_ids[cut_0_1] = index_1;
			cutting_point_ids[cut_0_2] = index_0_2;

			return;
		}

		//C: 0.5 0.0
		if (!b_0&&!b_1&&b_3)
		{
			new_points.push_back(v_0_1);
			int index_0_1 = new_points.size() - 1 + vecs.size();

			std::vector<int> triangle_0 = { index_0, index_0_1, index_2 };
			std::vector<int> triangle_1 = { index_1, index_2, index_0_1 };

			remesh_triangles.push_back(triangle_0);
			remesh_triangles_lables.push_back(lable_0);

			remesh_triangles.push_back(triangle_1);
			remesh_triangles_lables.push_back(lable_1);


			cutting_point_ids[cut_0_1] = index_0_1;
			cutting_point_ids[cut_0_2] = index_2;

			return;
		}

		//D: 0.5 0.5
		if (!b_0&&!b_1&&!b_2&&!b_3)
		{
			new_points.push_back(v_0_1);
			new_points.push_back(v_0_2);

			int index_0_1 = new_points.size() - 2 + vecs.size();
			int index_0_2 = new_points.size() - 1 + vecs.size();

			std::vector<int> triangle_0 = { index_0, index_0_1, index_0_2 };
			std::vector<int> triangle_1 = { index_1, index_0_2, index_0_1 };
			std::vector<int> triangle_2 = { index_2, index_0_2, index_1 };

			remesh_triangles.push_back(triangle_0);
			remesh_triangles_lables.push_back(lable_0);

			remesh_triangles.push_back(triangle_1);
			remesh_triangles_lables.push_back(lable_1);

			remesh_triangles.push_back(triangle_2);
			remesh_triangles_lables.push_back(lable_2);

			cutting_point_ids[cut_0_1] = index_0_1;
			cutting_point_ids[cut_0_2] = index_0_2;
			return;
		}

		//E: 0.5 1.0
		if (!b_0&&!b_1&&b_2&&!b_3)
		{
			new_points.push_back(v_0_1);
			int index_0_1 = new_points.size() - 1 + vecs.size();

			std::vector<int> triangle_0 = { index_0, index_0_1, index_2 };
			std::vector<int> triangle_1 = { index_1, index_2, index_0_1 };

			remesh_triangles.push_back(triangle_0);
			remesh_triangles_lables.push_back(lable_1);

			remesh_triangles.push_back(triangle_1);
			remesh_triangles_lables.push_back(lable_1);

			cutting_point_ids[cut_0_1] = index_0_1;
			cutting_point_ids[cut_0_2] = index_0;

			return;
		}

		//F: 1.0 0.5 
		if (b_0&&!b_1&&!b_2&&!b_3)
		{

			new_points.push_back(v_0_2);
			int index_0_2 = new_points.size() - 1 + vecs.size();

			std::vector<int> triangle_0 = { index_0, index_1, index_0_2 };
			std::vector<int> triangle_1 = { index_2, index_0_2, index_1 };

			remesh_triangles.push_back(triangle_0);
			remesh_triangles_lables.push_back(lable_1);

			remesh_triangles.push_back(triangle_1);
			remesh_triangles_lables.push_back(lable_1);

			cutting_point_ids[cut_0_1] = index_0;
			cutting_point_ids[cut_0_2] = index_0_2;

			return;
		}

		//G: 1.0 1.0
		if (b_0&&b_2)
		{
			remesh_triangles.push_back(triangle);
			remesh_triangles_lables.push_back(lable_1);


			cutting_point_ids[cut_0_1] = index_0;
			cutting_point_ids[cut_0_2] = index_0;

		}

		//H: 1.0 0.0
		if (b_0&&b_3)
		{
			remesh_triangles.push_back(triangle);
			remesh_triangles_lables.push_back(lable_1);

			cutting_point_ids[cut_0_1] = index_0;
			cutting_point_ids[cut_0_2] = index_2;
		}

		//I: 0.0 1.0
		if (b_1&&b_2)
		{
			remesh_triangles.push_back(triangle);
			remesh_triangles_lables.push_back(lable_1);

			cutting_point_ids[cut_0_1] = index_1;
			cutting_point_ids[cut_0_2] = index_0;
		}

		remesh_triangles.push_back(triangle);
		remesh_triangles_lables.push_back(lable_1);
		return;
	};

	Vector3d1 	remesh_points = vecs;
	Vector3d1 add_points;
	std::vector<std::vector<int>> remesh_triangles;
	std::vector<int> remesh_triangles_lables;

	std::vector<int> cutting_point_ids(cutting_points.size(),-1);

	for (int i = 0; i < face_id_0.size(); i++)
	{
		int index_0 = face_id_0[i];
		int index_1 = face_id_1[i];
		int index_2 = face_id_2[i];

		int lable_0 = two_class_lables[index_0];
		int lable_1 = two_class_lables[index_1];
		int lable_2 = two_class_lables[index_2];

		if (lable_0 == lable_1 && lable_0 == lable_2)
		{
			std::vector<int> triangle;
			triangle.push_back(index_0);
			triangle.push_back(index_1);
			triangle.push_back(index_2);

			remesh_triangles.push_back(triangle);
			remesh_triangles_lables.push_back(lable_0);
			continue;
		}
		if (lable_1 == lable_2) HandleBoundaryTriangle(vecs, two_class_lables, edges, cutting_points, remesh_triangles, remesh_triangles_lables, add_points, cutting_point_ids, index_0, index_1, index_2);
		if (lable_0 == lable_1) HandleBoundaryTriangle(vecs, two_class_lables, edges, cutting_points, remesh_triangles, remesh_triangles_lables, add_points, cutting_point_ids, index_2, index_0, index_1);
		if (lable_0 == lable_2) HandleBoundaryTriangle(vecs, two_class_lables, edges, cutting_points, remesh_triangles, remesh_triangles_lables, add_points, cutting_point_ids, index_1, index_2, index_0);
	}


	auto GetAddPointIndex=[](Vector3d1 &add_points, Vector3d v)
	{
		for (int i = 0; i < add_points.size(); i++)
		{
			double distance = CGAL_3D_Distance_Point_Point(add_points[i], v);
			if (Math::IsAlmostZero(distance))
			{
				return i;
			}
		}
		add_points.push_back(v);
		return (int)(add_points.size() - 1);
	};


	std::vector<int> new_index;
	Vector3d1 new_points;
	for (int i = 0; i < add_points.size(); i++)
		new_index.push_back(GetAddPointIndex(new_points, add_points[i]));

	for (int i = 0; i < new_points.size(); i++)
		remesh_points.push_back(new_points[i]);

	for (int i = 0; i < remesh_triangles.size(); i++)
	{
		for (int j = 0; j < remesh_triangles[i].size(); j++)
		{
			if (remesh_triangles[i][j] >= vecs.size())
				remesh_triangles[i][j] = vecs.size() + new_index[remesh_triangles[i][j] - vecs.size()];
		}
	}

	///////////////////////////////////////////////////////
	int nb = multi_cps[0];
	std::vector<std::vector<int>>  cutting_id_points;
	std::vector<int> points;
	for (int i = 0; i < cutting_point_ids.size(); i++)
	{
		//if (cutting_point_ids[i] < 0)
		//{
		//	std::cerr << "cutting_point_ids[i] < 0" << std::endl;
		//	system("pause");
		//}

		if (cutting_point_ids[i] >= 0)
		{
			if (cutting_point_ids[i] >= vecs.size())
				cutting_point_ids[i] = vecs.size() + new_index[cutting_point_ids[i] - vecs.size()];

			if (i >= nb)
			{
				//remove duplicate points
				auto RmoveDuplicatePoints = [](const std::vector<int> &points)
				{
					std::vector<int> ndps;
					for (int i = 0; i < points.size(); i++)
					{
						if (i == 0)
							ndps.emplace_back(points[i]);
						else
						{
							if (points[i] != ndps.back())
								ndps.emplace_back(points[i]);
						}
					}

					if (ndps.front() == ndps.back()) ndps.erase(ndps.begin());
					return ndps;
				};

				auto ndps = RmoveDuplicatePoints(points);
				if (ndps.size() < 3)
				{
					std::cerr << "if (ndps.size() < 3)" << std::endl;
					system("pause");
				}

				cutting_id_points.emplace_back(points);
				nb = multi_cps[cutting_id_points.size()];
				points.clear();
			}

			points.emplace_back(cutting_point_ids[i]);
		}
	}
	cutting_id_points.emplace_back(points);
	///////////////////////////////////////////////////////


	if (true)
	{
		for (int i = 0; i < cutting_id_points.size(); i++)
		{
			std::ofstream  export_file_output("Z:\\Documents\\Windows\\SmartSFC\\workspace\\CFS\\cutting_point_ids_"+std::to_string(i)+".obj");
			int export_index = 1;
			for (int j = 0; j < cutting_id_points[i].size(); j++)
			{
				CGAL_Export_Path_Segment(export_file_output, export_index, "intersection", 1,0,0,
					remesh_points[cutting_id_points[i][j]], remesh_points[cutting_id_points[i][(j + 1) % cutting_id_points[i].size()]], 0.05);
			}
			export_file_output.clear();
			export_file_output.close();
		}
	}

	auto TriangleMesh = [](const std::vector<int> &ids, const Vector3d1 &remesh_points,
		std::vector<std::vector<int>> &remesh_triangles, std::vector<int> &remesh_triangles_lables)
	{
		std::vector<int> nids;
		for (int i = 0; i < ids.size(); i++)
		{
			if (i == 0)
				nids.emplace_back(ids[i]);
			else
				if (ids[i] != nids.back())
					nids.emplace_back(ids[i]);
		}

		Vector3d1 points_3d;
		for (auto &id : nids) points_3d.emplace_back(remesh_points[id]);
		Vector2d1 points_2d = Math::Vector3d2d(points_3d);

		if (!CGAL_2D_Polygon_Simple(points_2d))
		{
			CGAL_2D_Polygon_Simple_0(points_2d);

			std::cerr << "if (!CGAL_2D_Polygon_Simple(points))" << std::endl;
			system("pause");
		}

		auto triangles = CGAL_2D_Polygon_Triangulation(points_2d);

		for (auto &triangle : triangles)
		{
			remesh_triangles.emplace_back(std::vector<int>{nids[triangle[2]], nids[triangle[1]], nids[triangle[0]]});
			remesh_triangles_lables.emplace_back(1);
		}
	};

	for (auto ids : cutting_id_points)
	{
		TriangleMesh(ids, remesh_points, remesh_triangles, remesh_triangles_lables);
	}


	//CGAL_Output_Boundary("G:\\abcd0.obj", remesh_points, remesh_triangles, remesh_triangles_lables);
	//CGAL_Output_Obj("G:\\remesh_full.obj", remesh_points, remesh_triangles);

	//CGAL_Output_Obj(path_0, remesh_points, remesh_triangles, remesh_triangles_lables, 0);
	CGAL_Output_Obj(output_path, remesh_points, remesh_triangles, remesh_triangles_lables, 1);
}
//std::vector<std::vector<int>> surface_vectices_to_vectices;

void ComputeEdgeLables(const int size_of_vertice, Halfedge_handle &start_hh, std::vector<std::pair<int,int>> &edges, std::vector<int> &lables)
{
	std::vector<int> two_class_lables(size_of_vertice, 0);
	std::vector<bool> edge_lables(size_of_vertice, false);
	std::vector<bool> over(size_of_vertice, false);
	for (int i = 0; i < edges.size(); i++)
	{
		edge_lables[edges[i].first] = true;
		edge_lables[edges[i].second] = true;
	}
	
	std::queue<Halfedge_handle> queue;
	std::queue<int> queue_1;

	queue.push(start_hh);
	queue_1.push(-1);
	over[start_hh->vertex()->id()] = true;

	std::vector<std::pair<Vector3d, int>> dsd;

	int iter = 0;
	while (queue.size() != 0)
	{
		Halfedge_handle hh = queue.front();
		int id = hh->vertex()->id();

		two_class_lables[hh->vertex()->id()] = 1;
		Vector3d v = PointVector3d(hh->vertex()->point());

		if (iter%100==0)
			std::cerr << iter << ": " << hh->vertex()->id() << " / " << size_of_vertice << std::endl;
		//////////
		
		dsd.emplace_back(std::pair<Vector3d, int>(v, queue_1.front()));

		if (false)
		{
			std::ofstream  export_file_output_0("Z:\\Documents\\Windows\\SmartSFC\\workspace\\CFS\\cutting\\" + std::to_string(iter) + "_" + std::to_string(hh->vertex()->id()) + "_" + std::to_string(queue_1.front()) + ".obj");
			int export_index_0 = 1;
			CGAL_Export_Path_Point(export_file_output_0, export_index_0,
				"inside" + std::to_string(iter) + "_" + std::to_string(hh->vertex()->id()) + "_" + std::to_string(queue_1.front()),
				1.0, 0.0, 0.0, v, 0.1);
			export_file_output_0.clear();
			export_file_output_0.close();
		}


		queue.pop();
		queue_1.pop();

		Halfedge_handle iter_hh = hh;
		do
		{
			int opposite_id = iter_hh->opposite()->vertex()->id();

			if (two_class_lables[opposite_id] == 0 && !over[opposite_id])
			{
				if (!(edge_lables[id] && edge_lables[opposite_id]))
				{
					queue.push(iter_hh->opposite());
					queue_1.push(iter);
					over[opposite_id] = true;
				}
				else
				{
					bool b = false;
					for (int i = 0; i < edges.size(); i++)
					{
						if ((edges[i].first == id&&edges[i].second == opposite_id) || (edges[i].second == id&&edges[i].first == opposite_id))
						{
							b = true;
							break;
						}
					}
					if (!b)
					{
						queue.push(iter_hh->opposite());
						queue_1.push(iter);
						over[opposite_id] = true;
					}
				}
			}

			if (!(edge_lables[id] && edge_lables[opposite_id]) && two_class_lables[opposite_id] == 0 && !over[opposite_id])
			{
				queue.push(iter_hh->opposite());
				queue_1.push(iter);
				over[opposite_id] = true;
			}

			iter_hh = iter_hh->opposite()->prev();
		} while (iter_hh != hh);

		iter++;
	}


	auto AAAA = [](std::vector<std::pair<Vector3d, int>> &dsd)
	{
		std::ofstream  export_file_output_0("Z:\\Documents\\Windows\\SmartSFC\\workspace\\CFS\\debug.obj");
		int export_index_0 = 1;
		int start_search = 22756;
		while (true)
		{
			if (start_search >= 0)
				CGAL_Export_Path_Point(export_file_output_0, export_index_0, "point", 1.0, 0.0, 0.0, dsd[start_search].first, 0.1);
			else
				break;

			start_search = dsd[start_search].second;
		}

		export_file_output_0.clear();
		export_file_output_0.close();
	};

	//AAAA(dsd);

	lables = two_class_lables;
}

//assumption: one projecting line can intersect with the three edges at most twice times
//            it's impossible to meet one edge for more than one times
void CGAL_Cut_Surface(Vector3d1 &boundary, Vector3d inside_point, std::string full_path, std::string output_path)
{
	Vector3d2 multi_boundary(1, boundary);
	CGAL_Cut_Surface_by_Multi_Boundaries(multi_boundary, inside_point, full_path, output_path);
	return;
}

//Cut a closed manifold mesh with boundaries
//multi_boundary: input boundaries
//inside_point: assign a point as the desired cutting surface
//full_path: file path of the input closed manifold mesh
//output_path: file path of output mesh
//assumption: one projecting line can intersect with the three edges of one triangle at most twice times
//            it's impossible to meet one edge for more than one times
void CGAL_Cut_Surface_by_Multi_Boundaries(Vector3d2 &multi_boundary, Vector3d inside_point, std::string full_path, std::string output_path)
{
	auto Intersection=[](const Halfedge_handle &hh, const int nb, const Vector3d inside, const Vector3d outside, Halfedge_handle &handle, Vector3d &intersection)
	{
		std::vector<Vector3d> result_vecs;
		std::vector<Halfedge_handle> result_handles;

		Halfedge_handle hh_3d = hh;
		for (int i = 0; i < nb; i++)
		{
			hh_3d = hh_3d->next();
			Vector3d edge_0 = PointVector3d(hh_3d->vertex()->point());
			Vector3d edge_1 = PointVector3d(hh_3d->opposite()->vertex()->point());

			double inside_d = CGAL_3D_Distance_Point_Segment(inside, edge_0, edge_1);
			if (Math::IsAlmostZero(inside_d))
			{
				result_vecs.emplace_back(inside);
				result_handles.emplace_back(hh_3d);
			}

			double outside_d = CGAL_3D_Distance_Point_Segment(outside, edge_0, edge_1);
			if (Math::IsAlmostZero(outside_d))
			{
				result_vecs.emplace_back(outside);
				result_handles.emplace_back(hh_3d);
			}
		}


		Point_3 p0 = hh->next()->next()->vertex()->point();
		Point_3 p1 = hh->vertex()->point();
		Point_3 p2 = hh->next()->vertex()->point();
		Vector_3 n = CGAL::cross_product(p2 - p1, p0 - p1);
		Vector3d nd(n.x(), n.y(), n.z());
		Math::SetVectorLength(nd, 1.0);
		Plane_3 plane(p1, Vector_3(nd[0], nd[1], nd[2]));

		//Point_2 point_to_2d(VectorPoint3d(inside), plane);

		Vector2d inside_2d = PointVector2d(point_to_2d(VectorPoint3d(inside), plane));
		Vector2d outside_2d = PointVector2d(point_to_2d(VectorPoint3d(outside), plane));

		hh_3d = hh;
		for (int i = 0; i < nb; i++)
		{
			hh_3d = hh_3d->next();
			Vector2d edge_0 = PointVector2d(point_to_2d(hh_3d->vertex()->point(), plane));
			Vector2d edge_1 = PointVector2d(point_to_2d(hh_3d->opposite()->vertex()->point(), plane));

			Vector3d edge_3d_0 = PointVector3d(hh_3d->vertex()->point());
			Vector3d edge_3d_1 = PointVector3d(hh_3d->opposite()->vertex()->point());

			Vector2d iter;

			if (CGAL_2D_Intersection_Segment_Segment(edge_0, edge_1, inside_2d, outside_2d, iter))
			{
				intersection = PointVector3d(point_to_3d(VectorPoint2d(iter), plane));
				intersection = CGAL_3D_Projection_Point_Segment(intersection, edge_3d_0, edge_3d_1);
				result_vecs.emplace_back(intersection);
				result_handles.emplace_back(hh_3d);
			}
		}

		if (result_vecs.empty())
		{
			return false;
		}
		else
		{
			double min_d = 1000000000000.0;
			for (int i = 0; i < result_vecs.size(); i++)
			{
				double dis = CGAL_3D_Distance_Point_Point(result_vecs[i],outside);
				if (min_d>dis)
				{
					min_d = dis;
					intersection = result_vecs[i];
					handle = result_handles[i];
				}
			}

			return true;
		}

	
	};

	//kd close query
	auto KD_Close_Query = [](kdtree *kd_tree, Vector3d query, const std::vector<Vector3d> &full_vecs,
		const std::vector<int> &full_face_id_0, const std::vector<int> &full_face_id_1, const std::vector<int> &full_face_id_2,
		const std::vector<std::vector<int>> &surface_vectices_to_face)
	{
		double *pos = new double[3];
		pos[0] = query[0];
		pos[1] = query[1];
		pos[2] = query[2];
		struct kdres *r = kd_nearest(kd_tree, pos);
		double position[3];
		int index = *(int*)kd_res_item(r, position);

		std::vector<int> faces = surface_vectices_to_face[index];
		std::vector<int> all_faces;
		int iter = 0;
		while (true)
		{
			for (int i = 0; i < faces.size(); i++) all_faces.emplace_back(faces[i]);

			std::vector<int> next_faces;
			for (int i = 0; i < faces.size(); i++)
			{
				auto vec_0 = full_face_id_0[faces[i]];
				auto vec_1 = full_face_id_1[faces[i]];
				auto vec_2 = full_face_id_2[faces[i]];
				auto d = CGAL_3D_Distance_Point_Triangle(query, full_vecs[vec_0], full_vecs[vec_1], full_vecs[vec_2]);
				if (Math::IsAlmostZero(d))return faces[i];

				for (int j = 0; j < surface_vectices_to_face[vec_0].size(); j++)
				{
					if (std::find(next_faces.begin(), next_faces.end(), surface_vectices_to_face[vec_0][j]) == next_faces.end() &&
						std::find(all_faces.begin(), all_faces.end(), surface_vectices_to_face[vec_0][j]) == all_faces.end())
					{
						next_faces.emplace_back(surface_vectices_to_face[vec_0][j]);
					}
				}

				for (int j = 0; j < surface_vectices_to_face[vec_1].size(); j++)
				{
					if (std::find(next_faces.begin(), next_faces.end(), surface_vectices_to_face[vec_1][j]) == next_faces.end() &&
						std::find(all_faces.begin(), all_faces.end(), surface_vectices_to_face[vec_1][j]) == all_faces.end())
					{
						next_faces.emplace_back(surface_vectices_to_face[vec_1][j]);
					}
				}

				for (int j = 0; j < surface_vectices_to_face[vec_2].size(); j++)
				{
					if (std::find(next_faces.begin(), next_faces.end(), surface_vectices_to_face[vec_2][j]) == next_faces.end() &&
						std::find(all_faces.begin(), all_faces.end(), surface_vectices_to_face[vec_2][j]) == all_faces.end())
					{
						next_faces.emplace_back(surface_vectices_to_face[vec_2][j]);
					}
				}
			}

			if (next_faces.empty()) break;

			if (iter > 5)break;


			faces = next_faces;
			iter++;
		}

		std::cerr << "KD_Close_Query" << std::endl;
		system("pause");

		return -1;
	};


	auto KD_Close_Query_0 = [](kdtree *kd_tree, Vector3d query, const std::vector<Vector3d> &full_vecs,
		const std::vector<int> &full_face_id_0, const std::vector<int> &full_face_id_1, const std::vector<int> &full_face_id_2,
		const std::vector<std::vector<int>> &surface_vectices_to_face)
	{
		double *pos = new double[3];
		pos[0] = query[0];
		pos[1] = query[1];
		pos[2] = query[2];
		struct kdres *r = kd_nearest(kd_tree, pos);
		double position[3];
		int index = *(int*)kd_res_item(r, position);

		std::vector<int> faces = surface_vectices_to_face[index];
		std::vector<int> all_faces;
		int iter = 0;
		int return_face_id = -1;
		double return_face_d = 1000000000000000.0;
		while (true)
		{
			for (int i = 0; i < faces.size(); i++) all_faces.emplace_back(faces[i]);

			std::vector<int> next_faces;
			for (int i = 0; i < faces.size(); i++)
			{
				auto vec_0 = full_face_id_0[faces[i]];
				auto vec_1 = full_face_id_1[faces[i]];
				auto vec_2 = full_face_id_2[faces[i]];
				auto d = CGAL_3D_Distance_Point_Triangle(query, full_vecs[vec_0], full_vecs[vec_1], full_vecs[vec_2]);
				//if (Math::IsAlmostZero(d))return faces[i];

				if (d < return_face_d)
				{
					return_face_d = d;
					return_face_id = faces[i];
				}

				for (int j = 0; j < surface_vectices_to_face[vec_0].size(); j++)
				{
					if (std::find(next_faces.begin(), next_faces.end(), surface_vectices_to_face[vec_0][j]) == next_faces.end() &&
						std::find(all_faces.begin(), all_faces.end(), surface_vectices_to_face[vec_0][j]) == all_faces.end())
					{
						next_faces.emplace_back(surface_vectices_to_face[vec_0][j]);
					}
				}

				for (int j = 0; j < surface_vectices_to_face[vec_1].size(); j++)
				{
					if (std::find(next_faces.begin(), next_faces.end(), surface_vectices_to_face[vec_1][j]) == next_faces.end() &&
						std::find(all_faces.begin(), all_faces.end(), surface_vectices_to_face[vec_1][j]) == all_faces.end())
					{
						next_faces.emplace_back(surface_vectices_to_face[vec_1][j]);
					}
				}

				for (int j = 0; j < surface_vectices_to_face[vec_2].size(); j++)
				{
					if (std::find(next_faces.begin(), next_faces.end(), surface_vectices_to_face[vec_2][j]) == next_faces.end() &&
						std::find(all_faces.begin(), all_faces.end(), surface_vectices_to_face[vec_2][j]) == all_faces.end())
					{
						next_faces.emplace_back(surface_vectices_to_face[vec_2][j]);
					}
				}
			}

			if (next_faces.empty()) break;

			if (iter > 5)break;


			faces = next_faces;
			iter++;
		}

		if (return_face_id >= 0)return return_face_id;

		std::cerr << "KD_Close_Query" << std::endl;
		system("pause");

		return -1;
	};


	auto AAA = [](std::string path,  Vector3d2 &multi_projects, std::vector<std::vector<Poly_facet_iterator>> &multi_project_faces,
		const std::vector<Vector3d> &full_vecs, const std::vector<int> &full_face_id_0, const std::vector<int> &full_face_id_1, const std::vector<int> &full_face_id_2)
	{
		std::ofstream  export_file_output("Z:\\Documents\\Windows\\SmartSFC\\workspace\\CFS\\projects.obj");
		int export_index_0 = 1;
		for (int i = 0; i < multi_projects.size(); i++)
		{
			for (int j = 0; j < multi_projects[i].size(); j++)
			{
				Vector3d end_0 = multi_projects[i][j];
				Vector3d end_1 = multi_projects[i][(j + 1) % multi_projects[i].size()];
				CGAL_Export_Path_Segment(export_file_output, export_index_0, "projection_" + Math::IntString(i), 1.0, 0.0, 0.0, end_0, end_1, 0.002);
			}
		}
		export_file_output.clear();
		export_file_output.close();

		std::vector<Vector3d> vecs;
		std::vector<int> face_id_0;
		std::vector<int> face_id_1;
		std::vector<int> face_id_2;
		for (int i = 0; i < multi_project_faces.size(); i++)
		{
			for (int j = 0; j < multi_project_faces[i].size(); j++)
			{
				Poly_facet_iterator cur_face = multi_project_faces[i][j];
				int face_id = cur_face->id();
				face_id_0.emplace_back(vecs.size());
				vecs.emplace_back(full_vecs[full_face_id_0[face_id]]);
				face_id_1.emplace_back(vecs.size());
				vecs.emplace_back(full_vecs[full_face_id_1[face_id]]);
				face_id_2.emplace_back(vecs.size());
				vecs.emplace_back(full_vecs[full_face_id_2[face_id]]);

			}
		}

		CGAL_Output_Obj("Z:\\Documents\\Windows\\SmartSFC\\workspace\\CFS\\triangles.obj", vecs, face_id_0, face_id_1, face_id_2);
	};
	
	auto BBBB = [](Halfedge_handle &cur_handle, Halfedge_handle &handle, Vector3d intersection,
		bool b, int iteration, int next_index, Vector3d1 &cutting_points, Vector3d2 &multi_projects, std::vector<std::vector<Poly_facet_iterator>> &multi_project_faces,  int i, Vector3d inside)
	{
		
		Point_3 p0 = cur_handle->next()->next()->vertex()->point();
		Point_3 p1 = cur_handle->vertex()->point();
		Point_3 p2 = cur_handle->next()->vertex()->point();
		Point_3 p3 = cur_handle->opposite()->vertex()->point();
		Vector3d v0 = PointVector3d(p0);
		Vector3d v1 = PointVector3d(p1);
		Vector3d v2 = PointVector3d(p2);
		Vector3d v3 = PointVector3d(p3);
		Vector3d center = (v0 + v1 + v2) / (float)3.0;

		//edge.push_back(handles[i]->vertex()->id());
		//edge.push_back(handles[i]->opposite()->vertex()->id());
		//handles
		std::string file_path;

		if (b)
			file_path = "Z:\\Documents\\Windows\\SmartSFC\\workspace\\CFS\\cutting\\cutting_points_" + Math::IntString(iteration) + "_" + Math::IntString(next_index) + "_1.obj";
		else
			file_path = "Z:\\Documents\\Windows\\SmartSFC\\workspace\\CFS\\cutting\\cutting_points_" + Math::IntString(iteration) + "_" + Math::IntString(next_index) + "_0.obj";


		std::ofstream  export_file_output(file_path);
		int export_index = 1;
		if (false)
		for (int j = 0; j < cutting_points.size() - 1; j++)
		{
			Vector3d end_0 = cutting_points[j];
			Vector3d end_1 = cutting_points[j + 1];
			CGAL_Export_Path_Segment(export_file_output, export_index, "cutting_points_" + Math::IntString(j), 1.0, 0.0, 0.0, end_0, end_1, 0.05);
		}
		//CGAL_Export_Path_Segment(export_file_output, export_index, "inside_outside_segment", 1.0, 0.0, 0.0, inside, multi_projects[i][next_index], 0.01);
		CGAL_Export_Path_Segment(export_file_output, export_index, "cur_tri_edge", 1.0, 0.0, 0.0, v1, v3, 0.001);
		CGAL_Export_Path_Point(export_file_output, export_index, "cur_tri_center_" + std::to_string(cur_handle->face()->id()), 1.0, 0.0, 0.0, center, 0.002);
		CGAL_Export_Path_Point(export_file_output, export_index, "cur_tri_inside", 1.0, 0.0, 0.0, inside, 0.002);
		CGAL_Export_Path_Point(export_file_output, export_index, "cur_tri_outside_" + std::to_string(multi_project_faces[i][next_index]->id()), 
			1.0, 0.0, 0.0, multi_projects[i][next_index], 0.002);

		if (b)
		{
			Vector3d vv0 = PointVector3d(handle->vertex()->point());
			Vector3d vv1 = PointVector3d(handle->opposite()->vertex()->point());
			CGAL_Export_Path_Segment(export_file_output, export_index, "next_edge", 1.0, 0.0, 0.0, vv0, vv1, 0.001);

			CGAL_Export_Path_Point(export_file_output, export_index, "inters_point", 1.0, 0.0, 0.0, intersection, 0.002);
		}

		export_file_output.clear();
		export_file_output.close();
	};

	auto CCCC = [](Vector3d1 &cutting_points)
	{
		std::ofstream  export_file_output("Z:\\Documents\\Windows\\SmartSFC\\workspace\\CFS\\cutting_points.obj");
		int export_index = 1;
		for (int i = 0; i < cutting_points.size() - 1; i++)
		{
			CGAL_Export_Path_Segment(export_file_output, export_index, "intersection_" + Math::IntString(i), 1.0, 0.0, 0.0, cutting_points[i], cutting_points[i + 1], 0.05);
		}
		export_file_output.clear();
		export_file_output.close();
	};

	auto DDDD = [](std::string path, std::vector<std::pair<int, int>> &edges,
		std::vector<Halfedge_handle> &handles, Vector3d1 &cutting_points, Vector3d1 &full_vecs)
	{
		std::ofstream  export_file_output_0(path);
		int export_index_0 = 1;

		for (int i = 0; i < handles.size(); i++)
		{
			auto name = "edge_" + std::to_string(i) + "_" + std::to_string(handles[i]->face()->id());
			CGAL_Export_Path_Segment(export_file_output_0, export_index_0, name, 1.0, 0.0, 0.0, full_vecs[edges[i].first], full_vecs[edges[i].second], 0.03);
			CGAL_Export_Path_Point(export_file_output_0, export_index_0, name, 1.0, 0.0, 0.0, cutting_points[i], 0.05);

			Point_3 p0 = handles[i]->next()->next()->vertex()->point();
			Point_3 p1 = handles[i]->vertex()->point();
			Point_3 p2 = handles[i]->next()->vertex()->point();
			Vector3d v0 = PointVector3d(p0);
			Vector3d v1 = PointVector3d(p1);
			Vector3d v2 = PointVector3d(p2);
			Vector3d center = (v0 + v1 + v2) / (float)3.0;
			CGAL_Export_Path_Point(export_file_output_0, export_index_0, name, 1.0, 0.0, 0.0, center, 0.05);

		}

		export_file_output_0.clear();
		export_file_output_0.close();
	};

	auto EEEE = [](Vector3d inside)
	{
		std::ofstream  export_file_output_0("Z:\\Documents\\Windows\\SmartSFC\\workspace\\CFS\\inside.obj");
		int export_index_0 = 1;
		CGAL_Export_Path_Point(export_file_output_0, export_index_0, "start_point", 1.0, 0.0, 0.0, inside, 0.05);
		export_file_output_0.clear();
		export_file_output_0.close();
	};

	auto FFFF = [](std::string path,Vector3d1 insides)
	{
		std::ofstream  export_file_output_0(path);
		int export_index_0 = 1;
		for (auto inside:insides)
			CGAL_Export_Path_Point(export_file_output_0, export_index_0, "start_point", 1.0, 0.0, 0.0, inside, 0.05);
		export_file_output_0.clear();
		export_file_output_0.close();
	};

	//build polyhedron from the full_path obj
	Polyhedron_3 polyhedron;
	Vector3d1 full_vecs;
	std::vector<int> full_face_id_0;
	std::vector<int> full_face_id_1;
	std::vector<int> full_face_id_2;
	Construct_Polyhedron(polyhedron, full_path, full_vecs, full_face_id_0, full_face_id_1, full_face_id_2);

	if (full_vecs.empty() || full_face_id_0.empty() || full_face_id_1.empty() || full_face_id_2.empty())
	{
		std::cerr << "if (full_vecs.empty() || full_face_id_0.empty() || full_face_id_1.empty() || full_face_id_2.empty())" << std::endl;
		system("pause");
	}

	//surface normals
	//full_face_iters
	Vector3d1 surface_normals;
	std::vector<Poly_facet_iterator> full_face_iters;
	for (Poly_facet_iterator iter = polyhedron.facets_begin(); iter != polyhedron.facets_end(); iter++)
	{
		full_face_iters.emplace_back(iter);
		
		Point_3 p0 = iter->halfedge()->next()->next()->vertex()->point();
		Point_3 p1 = iter->halfedge()->vertex()->point();
		Point_3 p2 = iter->halfedge()->next()->vertex()->point();
		Vector_3 n = CGAL::cross_product(p2 - p1, p0 - p1);
		Vector3d nd(n.x(), n.y(), n.z());
		Math::SetVectorLength(nd, 1.0);
		surface_normals.emplace_back(nd);
	}

	//surface_vectices_to_face and surface_vectices_to_vectices
	std::vector<std::vector<int>> surface_vectices_to_face;
	CGAL_3D_Triangle_Mesh_Vecs_Faces(full_vecs, full_face_id_0, full_face_id_1, full_face_id_2, surface_vectices_to_face);
	std::vector<std::vector<int>> surface_vectices_to_vectices;
	CGAL_3D_Triangle_Mesh_Vecs_Neighbors(full_vecs, full_face_id_0, full_face_id_1, full_face_id_2, surface_vectices_to_vectices);

	//kd tree
	std::vector<int> full_vec_ids(full_vecs.size(), 0);
	for (int i = 0; i < full_vecs.size(); i++) full_vec_ids[i] = i;
	kdtree *kd_tree = kd_create(3);
	for (int i = 0; i < full_vecs.size(); i++)
	{
		void* val = &full_vec_ids[i];
		kd_insert3(kd_tree, full_vecs[i][0], full_vecs[i][1], full_vecs[i][2], val);
	}

	//project points along input boundaries to the full_path obj
	Vector3d2 multi_projects;// projecting points
	std::vector<std::vector<Poly_facet_iterator>> multi_project_faces;//related faces of projecting points
	for (int i = 0; i < multi_boundary.size(); i++)
	{
		Vector3d1 one_projects;
		std::vector<Poly_facet_iterator> one_project_faces;
		for (int j = 0; j < multi_boundary[i].size(); j++)
		{
			std::cerr << i << " " << j << std::endl;
			Poly_point_3 query(multi_boundary[i][j][0], multi_boundary[i][j][1], multi_boundary[i][j][2]);
			auto fid = KD_Close_Query_0(kd_tree, multi_boundary[i][j], full_vecs, full_face_id_0, full_face_id_1, full_face_id_2, surface_vectices_to_face);
			one_projects.push_back(multi_boundary[i][j]);
			one_project_faces.push_back(full_face_iters[fid]);
		}
		multi_projects.push_back(one_projects);
		multi_project_faces.push_back(one_project_faces);
	}
	//output the projecting points
	if (DebugInformation())
		AAA("Z:\\Documents\\Windows\\SmartSFC\\workspace\\CFS\\projects.obj",multi_projects, multi_project_faces, full_vecs, full_face_id_0, full_face_id_1, full_face_id_2);

	Math::ClearFolder("Z:\\Documents\\Windows\\SmartSFC\\workspace\\CFS\\cutting\\");

	//searching for all of the cutting points on edges
	Vector3d1 cutting_points;
	std::vector<int> multi_cps;
	std::vector<Halfedge_handle> handles;
	std::vector<bool> face_used(full_face_id_0.size(),false);
	for (int i = 0; i < multi_projects.size(); i++)
	{
		std::cerr << "multi_projects: " << i <<" / "<< multi_projects.size()<< std::endl;

		int cur_face_id = multi_project_faces[i][0]->id();
		face_used[cur_face_id] = true;
		Poly_facet_iterator cur_face = multi_project_faces[i][0];
		Halfedge_handle cur_handle = cur_face->halfedge();
		Vector3d inside = multi_projects[i][0];

		int next_index = 1;
		int iteration = 0;
		while (true)
		{
			if (multi_projects.size()==1) std::cout << iteration << " / " << multi_projects[i].size() << " " << next_index << std::endl;

			//searching for the outside point of the current triangle
			bool goon = false;
			while (true)
			{
				if (cur_face_id == multi_project_faces[i][next_index]->id())
				{
					inside = multi_projects[i][next_index];
					if (next_index == 0) break;
					next_index = (next_index + 1) % multi_projects[i].size();
				}
				else
				{
					auto cur_face_normal = surface_normals[cur_face_id];
					auto next_face_normal = surface_normals[multi_project_faces[i][next_index]->id()];
					double angle = Math::GetAngleBetween(cur_face_normal, next_face_normal);

					if (angle<Math::Math_PI/2.0 && OutsidePointInsideTriangle(cur_face, multi_projects[i][next_index]))
					{
						inside = multi_projects[i][next_index];
						if (next_index == 0) break;
						next_index = (next_index + 1) % multi_projects[i].size();
					}
					else
					{
						//if (face_used[multi_project_faces[i][next_index]->id()])
						//{
						//	inside = multi_projects[i][next_index];
						//	if (next_index == 0) break;
						//	next_index = (next_index + 1) % multi_projects[i].size();
						//}
						goon = true;
						break;
					}
				}
			}


			if (!goon) break;
			Halfedge_handle handle;
			Vector3d intersection;
			bool b;
			if (iteration == 0)
				b = Intersection(cur_handle, 3, inside, multi_projects[i][next_index], handle, intersection);
			else
				b = Intersection(cur_handle, 2, inside, multi_projects[i][next_index], handle, intersection);
			
			//output the cutting points of each iteration
			if (true)
				BBBB(cur_handle, handle, intersection, b, iteration, next_index, cutting_points, multi_projects, multi_project_faces, i, inside);

			if (b)
			{
				cutting_points.push_back(intersection);
				handles.push_back(handle);

				//move next step
				inside = intersection;
				cur_handle = handle->opposite();
				cur_face_id = cur_handle->face()->id();
				face_used[cur_face_id] = true;
				cur_face = cur_handle->face();
				
				if (cur_face_id == multi_project_faces[i][0]->id())
				{
					break;
				}
			}
			else
			{
				std::cerr << "if (b)" << std::endl;
				system("pause");

				cutting_points.erase(cutting_points.begin() + cutting_points.size() - 1);
				handles.erase(handles.begin() + handles.size() - 1);

				next_index = (next_index + 1) % multi_projects[i].size();
				if (next_index == 0) 
					break;

				inside = cutting_points[cutting_points.size() - 1];
				cur_handle = handles[handles.size() - 1]->opposite();
				cur_face_id = cur_handle->face()->id();
				cur_face = cur_handle->face();
			}
			iteration++;
		}

		multi_cps.emplace_back(cutting_points.size());
	}

	//output the cutting points
	if (false)
	{
		CCCC(cutting_points);
	}

	/*******************************************/
	std::vector<std::pair<int,int>> edges;
	/*******************************************/
	for (int i = 0; i < handles.size(); i++)
		edges.emplace_back(handles[i]->vertex()->id(), handles[i]->opposite()->vertex()->id());

	//smoothing
	//for (int i = 0; i < 3; i++) OneIterationSmoothBoundary(full_vecs, edges, cutting_points);

	auto fid = KD_Close_Query(kd_tree, inside_point, full_vecs, full_face_id_0, full_face_id_1, full_face_id_2, surface_vectices_to_face);

	if (true)
	{
		DDDD("Z:\\Documents\\Windows\\SmartSFC\\workspace\\CFS\\edges.obj",edges, handles, cutting_points, full_vecs);
		EEEE(PointVector3d(full_face_iters[fid]->halfedge()->vertex()->point()));
	}

	//lables
	std::vector<int> two_class_lables;
	ComputeEdgeLables(full_vecs.size(), full_face_iters[fid]->halfedge(), edges, two_class_lables);

	Vector3d1 lable_0_vecs;
	Vector3d1 lable_1_vecs;
	for (int i = 0; i < two_class_lables.size(); i++)
	{
		if (two_class_lables[i] == 0) lable_0_vecs.emplace_back(full_vecs[i]);
		if (two_class_lables[i] == 1) lable_1_vecs.emplace_back(full_vecs[i]);
	}
	std::cerr << lable_0_vecs.size() << " " << lable_1_vecs.size() << std::endl;

	FFFF("Z:\\Documents\\Windows\\SmartSFC\\workspace\\CFS\\lable_0_vecs.obj", lable_0_vecs);
	FFFF("Z:\\Documents\\Windows\\SmartSFC\\workspace\\CFS\\lable_1_vecs.obj", lable_1_vecs);

	//remesh triangles
	ComputeRemeshTriangles(full_vecs, full_face_id_0, full_face_id_1, full_face_id_2, 
		edges, cutting_points, multi_cps, two_class_lables, output_path);

	std::cout << "over" << std::endl;
	return;
}


//
//
//void main()
//{
//
//	//Vector3d1 heat_surface_vertices;
//	//std::vector<int> heat_face_id_0;
//	//std::vector<int> heat_face_id_1;
//	//std::vector<int> heat_face_id_2;
//	//CGAL_3D_Read_Triangle_Mesh("D:\\CNCProduction\\NewResults_Pieces\\kitten_210\\kitten_210\\path\\0_heat_sub_split.obj", heat_surface_vertices, heat_face_id_0, heat_face_id_1, heat_face_id_2);
//	//std::vector<double> vertices_d;
//	//std::vector<int> vertices_index;
//	//std::ifstream file("D:\\CNCProduction\\NewResults_Pieces\\kitten_210\\kitten_210\\path\\heat.dist");
//	//for (int i = 0; i < heat_surface_vertices.size(); i++){
//	//	double d;
//	//	int index;
//	//	file >> d >> index;
//	//	vertices_d.push_back(d);
//	//	vertices_index.push_back(index);
//	//}
//	//file.clear();
//	//file.close();
//
//	//double distance  = ComputeGapFromScallop(0.0, 2.0, 0.02) / 2.0;
//
//	//Vector3d2 offsets;
//	//CGAL_3D_Mesh_Extract_Isoline(heat_surface_vertices, heat_face_id_0, heat_face_id_1, heat_face_id_2, vertices_d, distance, offsets);
//
//
//	/*std::string path = "D:\\CNCProduction\\NewResults\\maxplanck\\maxplanck";
//
//	Vector3d1 surface_vertices;
//	std::vector<int> face_id_0;
//	std::vector<int> face_id_1;
//	std::vector<int> face_id_2;
//
//	CGAL_3D_Read_Triangle_Mesh(path + ".off", surface_vertices, face_id_0, face_id_1, face_id_2);
//
//	std::vector<double> max_curvature;
//	std::vector<double> min_curvature;
//	CGAL_3D_Mesh_Curvature(surface_vertices, face_id_0, face_id_1, face_id_2, max_curvature, min_curvature);
//
//	Vector3d1 colors;
//	for (int i = 0; i < surface_vertices.size(); i++)
//	{
//		if (min_curvature[i] < -0.5)
//		{
//			Vector3d color(1.0, 0.0, 0.0);
//			colors.push_back(color);
//		}
//		else
//		{
//			Vector3d color(0.0, 1.0, 0.0);
//			colors.push_back(color);
//		}
//	}
//
//	CGAL_Output_Obj(path + "_0.5_curvature.obj", surface_vertices, colors, face_id_0, face_id_1, face_id_2);
//	CGAL_Mesh_Laplace_Smooth_by_Curvature(surface_vertices, face_id_0, face_id_1, face_id_2, -0.5);
//
//	std::vector<double>().swap(max_curvature);
//	std::vector<double>().swap(min_curvature);
//	Vector3d1().swap(colors);
//	CGAL_3D_Mesh_Curvature(surface_vertices, face_id_0, face_id_1, face_id_2, max_curvature, min_curvature);
//
//
//	for (int i = 0; i < surface_vertices.size(); i++)
//	{
//		if (min_curvature[i] < -0.5)
//		{
//			Vector3d color(1.0, 0.0, 0.0);
//			colors.push_back(color);
//		}
//		else
//		{
//			Vector3d color(0.0, 1.0, 0.0);
//			colors.push_back(color);
//		}
//	}
//
//
//	CGAL_Output_Obj(path + "_0.5_curvature_after.obj", surface_vertices, colors, face_id_0, face_id_1, face_id_2);*/
//
//
//	Vector3d inside;
//	Vector3d2 boundarys;
//	CGAL_3D_Triangle_Mesh_Boundary("D:\\CNCProduction\\Revision\\engineering\\cad2\\cad2\\0_309.obj", boundarys, inside);
//	CGAL_Cut_Surface_by_Multi_Boundaries(boundarys, inside, "D:\\CNCProduction\\Revision\\engineering\\cad2\\cad2.off", "G:\\cut.obj");
//
//	
//	//std::string path = "D:\\CNCProduction\\new_setup\\threepieces\\404\\404";
//	//std::string obj_path = path + "\\path\\0_heat_sub_split.obj";
//	//std::string obj_full_path = path + "\\0_full.off";
//	//std::string heat_path = path + "\\path\\heat.dist";
//	//std::string scale_path = path + "\\path\\scale.source";
//	//Heat_to_Gradient(obj_path, obj_full_path, heat_path, scale_path);
//	//Heat_to_Gradient("G:\\remesh_2_piece.obj", "G:\\0_404_full_0.off", "G:\\heat.dist", "G:\\0_heat_factor1.source");
//	//
//	system("pause");
//}
