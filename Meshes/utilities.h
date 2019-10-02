//
//  utilities.h
//  RIS
//
//  Created by Vitaliy Kurlin on 21/10/2016.
//  Copyright © 2016 Vitaliy Kurlin. All rights reserved.
//
#ifndef utilities_h
#define utilities_h
#pragma once

#include <fstream>
#include <iostream>
#include <deque>
#include <set>
#include <string>
#include <iterator>
#include <utility>
#include <algorithm>
#include <limits>

// OpenCV
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
using namespace cv;

//Boost
//#define BOOST_FILESYSTEM_NO_DEPRECATED
#include <boost/graph/graphviz.hpp>
#include "boost/graph/topological_sort.hpp"
#include <boost/graph/graph_traits.hpp>
#include "boost/filesystem/operations.hpp"
#include "boost/filesystem/path.hpp"
#include "boost/progress.hpp"
#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/prim_minimum_spanning_tree.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <boost/graph/filtered_graph.hpp>

// CGAL
#include <CGAL/enum.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/boost/graph/graph_traits_Delaunay_triangulation_2.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/Polygon_2.h>

#include "segmentation.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Delaunay_triangulation_2<K> Triangulation;

const double Min_Angle = M_PI / 6;
//
class Node1
{
public:
    int index; // index of triangle in faces
    int uplink;
    double birth;
    int height; // height of the tree rooted at Node
    int bar; // index of the region containing the triangle of this node
    std::vector<int> live; // indices of live triangles from Tree of this node
    std::vector<int> edges; // indices of edges of the corresponding triangle in DelEdges
    
    Node1 ()
    {
        index = 0;  // index in faces: 0 means external boundary
        uplink = 0; // index of the parent node in faces
        birth = 0; //
        height = 1; // no nodes below the new node
        bar = -1;
        live.clear();
        edges.resize(3);
    }
};

typedef CGAL::Triangulation_vertex_base_with_info_2<int, K> VB;
typedef CGAL::Triangulation_face_base_with_info_2<Node1, K> FB;
typedef CGAL::Triangulation_data_structure_2<VB,FB> TDS;
typedef CGAL::Delaunay_triangulation_2<K, TDS> DT;
typedef K::Point_2 P2;
typedef K::Segment_2 Segment;
typedef CGAL::Polygon_2<K> Polygon;
typedef DT::Edge DE;
typedef DT::Face DF;
typedef DT::Vertex_handle VH;
typedef DT::Face_handle FH;
typedef DT::Vertex_iterator VI;
typedef DT::Edge_iterator EI;
typedef DT::Finite_vertices_iterator FVI;
typedef DT::Finite_edges_iterator FEI;
typedef DT::All_faces_iterator AFI;
typedef DT::Finite_faces_iterator FFI;
//

struct Compare_Points
{
    bool operator()(const Point& a, const Point& b) const
    {
        if ( a.x < b.x ) return true;
        if ( a.x > b.x ) return false;
        return ( a.y < b.y );
    }
};

String Str (int n);
string Round (double v);
string Round1 (double v);
string Round2 (double v);
string Round3 (double v);
string Round4 (double v);
string Round5 (double v);

class Index_Value
{
public:
    int index;
    double value;
    Index_Value() {}
    Index_Value(int i, double v) { index = i; value = v; }
};
bool Decreasing_Values (Index_Value const& p1, Index_Value const& p2);
bool Increasing_Values (Index_Value const& p1, Index_Value const& p2);

class Benchmark
{
public:
    double num_faces = 0, num_edges = 0, num_vertices = 0, BR2 = 0, BR1 = 0, CUE = 1, USE = 1, ASA = 1, comp = 0, RMS = 0, time = 0;
    Benchmark (double v) { CUE = v; USE = v; ASA = v; }
    void Add (size_t v, size_t e, size_t f) { num_faces += f; num_edges += e; num_vertices += v;}
    void Add (Benchmark const& b)
    {
        Add( b.num_vertices, b.num_edges, b.num_faces );
        BR2 += b.BR2; BR1 += b.BR1; CUE += b.CUE; USE += b.USE; ASA += b.ASA; RMS += b.RMS; comp += b.comp; time += b.time;
    }
    void Divide_by_Pixels (int n) { CUE *= 1.0/n; USE *= 1.0/n; ASA *= 1.0/n; }
    void Convert_in_Percents() { BR2 *= 100; CUE *= 100; ASA *= 100; }
    void Divide_by_Images (int m)
    {
        num_faces *= 1.0 / m; num_edges *= 1.0 / m; num_vertices *= 1.0 / m;
        BR2 /= m; BR1 /= m; comp /= m; RMS /= m; time /= m;
        Divide_by_Pixels( m );
    }
    void Print (std::ostream& fout)
    {
        //fout<<" f="<<num_faces<<" e="<<num_edges; //<<" v="<<num_vertices;
        fout<<" BR2="<<BR2<<" CUE="<<CUE<<" ASA="<<ASA<<" Comp="<<Round1( comp )<<" RMS="<<RMS; //
    }
    void Print_All (int n, std::ostream& fout)
    {
        String s = "\ncurrent t=" + Round1( time / n )+ " f=" + Round1( num_faces / n ) + " e=" + Round1( num_edges / n ) + " BR2=" + Round3( BR2 / n ) + " CUE=" + Round3( CUE / n ) + " ASA=" + Round3( ASA / n ) + " Comp=" + Round1( comp / n ) + " RMS=" + Round3( RMS / n );
        fout << s;
    }
    void Print_All (int n, std::ostream& fout1, std::ostream& fout2 ) { Print_All( n, fout1 ); Print_All( n, fout2 ); }
    void Print () { Print( std::cout ); }
    void Print (std::ostream& fout1, std::ostream& fout2 ) { Print( fout1 ); Print( fout2 ); }
};

class iSegmentation
{
public:
    int num_faces, num_edges, num_vertices;
    std::vector<Point> boundary;
    Mat_<bool> boundary_mask;
    Mat_<int32_t> indices;
    std::vector< std::vector<Point> > superpixels; // list of pixels within every superpixel
};

void Write_Image (String const& name, Mat const& image);

void Print (String s);
void Print (String s, std::ostream& fout);
void Print (String s, std::ostream& file1, std::ostream& file2);
void Create_Directory (String name);

bool Find_Sine_Cosine (Point2d a, Point2d b, Point2d c, double& sine, double& cosine);
bool Collinear (Point2d a, Point2d b, Point2d c);
void Points_to_Line_Coef (Point2d p1, Point2d p2, double& a, double& b, double& c);
bool Distance_Point_to_Line (Point2d p, double a, double b, double c, double& d);
bool Is_Nonstrictly_Inside_Polygon (Polygon const& poly, P2 point);

void Scale_Image (const Mat& img, Mat& img_scaled, int scale);
void Transpose (const Mat_<int>& matrix, Mat_<int>& transpose);
bool Read_Image_Names_Sizes (std::string file_name, std::vector<std::string>& image_names, std::vector<Point>& image_sizes);

bool Read_BSD_Human (std::string path, iSegmentation& segmentation);
bool Read_BSD_Human (std::string path, Segmentation& segmentation);
void Boundary_Recall (iSegmentation const& s0, iSegmentation const& s1, Benchmark& b);
void Undersegmentation_Errors (iSegmentation const& s0, iSegmentation const& s1, Benchmark& benchmark);
bool Find_Benchmarks (PixelSegmentation& segmentation, String name, Point sizes, Benchmark& average, std::ostream& fout);


#endif /* utilities_h */
