//
//  mesh.h
//  RIS
//
//  Created by Vitaliy Kurlin on 21/10/2016.
//  Copyright © 2016 Vitaliy Kurlin. All rights reserved.
//
#ifndef mesh_h
#define mesh_h
#pragma once

//#include "colors.h"
#include "utilities.h"

// OpenMesh
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>

struct MyTraits : public OpenMesh::DefaultTraits
{
    typedef OpenMesh::Vec3f Normal;
    
    //couldn't get it to work with this double precision when calculating edge lengths
    //typedef OpenMesh::Vec3d Point;
    VertexTraits
    {
    public:
        //bool corner = false;
    };
    EdgeTraits
    {
    public:
        double length = -1;
        int index = -1;
        bool shrink = false, chain_found = false;
    };
    FaceTraits
    {
        int label = -1;
        double area = 0, perimeter = 0, mean = 0;
        Vec3b color = Vec3b( 0, 0, 0 );
        //bool passed = false;
    };
    VertexAttributes(OpenMesh::Attributes::Status);
    FaceAttributes(OpenMesh::Attributes::Status | OpenMesh::Attributes::Normal);
    EdgeAttributes(OpenMesh::Attributes::Status);
};

typedef OpenMesh::PolyMesh_ArrayKernelT<MyTraits> Mesh;
typedef Mesh::Point MP;
typedef Mesh::VertexHandle MVH;
typedef Mesh::EdgeHandle MEH;
typedef Mesh::HalfedgeHandle MHH;
typedef Mesh::FaceHandle MFH;

MP P2M (P2 p);
P2 MP2 (MP p);
Point2d P (MP p);
MP P (Point2d p);
P2 P2d (Point2d p);

bool Print_Face (Mesh& mesh, MFH fh);
void Print_Vertex (Mesh& mesh, MVH vh);
void Print_Vertices (Mesh& mesh);
void Print_Neighbors (Mesh& mesh, MVH vh);
void Print_Info (Mesh& mesh);

double Area_Face (Mesh& mesh, MFH fh);
double Area_Faces (Mesh& mesh);
void Read_Mesh (Mesh& mesh, String file);
void Write_Mesh (Mesh& mesh, String file);
void Scale_Mesh (Mesh& mesh, double scale);
bool Enlarge_Mesh (Mesh& mesh, Size size);
bool Draw_Mesh (Mesh& mesh, int scale, int shift, Scalar edge_color, Scalar vertex_color, int thickness, Mat& image);
bool Draw_Reconstruction (Mesh& mesh, int scale, Mat& image);
bool Type_Labels (Mesh& mesh, int scale, int shift, Scalar color, int thickness, Mat& image);

bool Remove_Trivial_Vertices (Mesh& mesh);
bool Remove_Small_Faces (Mesh& mesh, double min_area_face);
bool Straighten_Chains (Mesh& mesh, double color_offset, double min_color_dif);
bool Merge_Faces (Mesh& mesh, double min_color_dif);

void Update_Labels (Mesh& mesh, MFH fh, Mat_<int>& labels);
void Find_Compactness (Mesh& mesh, int total, double& comp);
void Find_Means (const Mat& image, const Mat_<int>& labels, Mesh& mesh);
bool Find_Benchmarks (Mesh& mesh, String name, Point sizes, Benchmark& average, std::ostream& fout);
bool Find_Benchmarks (const Mat& image, PixelSegmentation& s, Mesh& mesh, String path, Point sizes, Benchmark& benchmark_segm, Benchmark& benchmark_mesh);

#endif /* mesh_h */

