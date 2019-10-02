//
//  mesh.cpp
//  RIS
//
//  Created by Vitaliy Kurlin on 21/10/2016.
//  Copyright © 2016 Vitaliy Kurlin. All rights reserved.
//

#include "utilities.h"
#include "segmentation.h"
#include "mesh.h"
//#include "colors.h"

Point2d P (MP p) { return Point2d( p[0], p[1] ); }
MP P (Point2d p) { return MP( p.x, p.y, 0 ); }
MP P2M (P2 p) { return MP( p.x(), p.y(), 0 ); }
P2 MP2 (MP p) { return P2( p[0], p[1] ); }
P2 P2d (Point2d p) { return P2( p.x, p.y ); }

void Print_Vertex (Mesh& mesh, MVH vh) { std::cout<<" v"<<P( mesh.point( vh ) )<<"deg"<<mesh.valence( vh ); }

void Print_Edge (Mesh& mesh, MHH h)
{
    Print_Vertex( mesh, mesh.from_vertex_handle( h ) );
    Print_Vertex( mesh, mesh.to_vertex_handle( h ) );
}

void Print_Edge (Mesh& mesh, MEH eh) { Print_Edge( mesh, mesh.halfedge_handle( eh, 0 ) ); }

bool Print_Face (Mesh& mesh, MFH fh)
{
    std::cout<<"\nf"<<mesh.data( fh ).label<<" val="<<mesh.valence( fh )<<" area="<<mesh.data( fh ).area;
    for ( auto fv_it = mesh.fv_iter( fh ); fv_it.is_valid(); fv_it++ ) Print_Vertex( mesh, *fv_it );
    return true;
}

bool Print_Face_Vertices (Mesh& mesh, MFH fh)
{
    std::cout<<"\nf"<<mesh.data( fh ).label<<" val="<<mesh.valence( fh )<<" area="<<mesh.data( fh ).area;
    for ( auto fv_it = mesh.fv_iter( fh ); fv_it.is_valid(); fv_it++ ) std::cout<<" "<<P( mesh.point( *fv_it ) );
    return true;
}

void Print_Face (Mesh& mesh, int label)
{
    for ( auto f_it = mesh.faces_begin(); f_it != mesh.faces_end(); f_it++ )
        if ( mesh.data( *f_it ).label == label ) { Print_Face( mesh, *f_it ); break; }
}

void Print_Neighbors (Mesh& mesh, MVH vh)
{
    Print_Vertex( mesh, vh );
    std::cout<<" neighbors: ";
    for ( auto v_it = mesh.vv_iter( vh ); v_it.is_valid(); v_it++ )
        std::cout<<P( mesh.point( *v_it ) );
}

void Print_Vertices (Mesh& mesh)
{
    for ( auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++ )
        if ( mesh.valence( *v_it ) == 4 )
            Print_Vertex( mesh, *v_it );
}

void Print_Info (Mesh& mesh)
{
    std::cout<<" f="<<mesh.n_faces()<<" e="<<mesh.n_edges(); //<<" v="<<mesh.n_vertices()
}

void Read_Mesh (Mesh& mesh, String file)
{
    if ( ! OpenMesh::IO::read_mesh( mesh, file ) ) std::cout << "\nError reading "<<file;
}

void Write_Mesh (Mesh& mesh, String file)
{
    if ( ! OpenMesh::IO::write_mesh( mesh, file ) ) std::cout << "\nError writing "<<file;
}

void Find_Compactness (Mesh& mesh, int total, double& compactness)
{
    MP p0, p1;
    compactness = 0;
    double perimeter = 0, l = 0;
    for ( auto f_it = mesh.faces_begin(); f_it != mesh.faces_end(); f_it++ )
    {
        if ( ! mesh.is_valid_handle( *f_it ) or mesh.valence( *f_it ) < 3 )  continue;
        perimeter = 0;
        for ( auto fh_it = mesh.fh_iter( *f_it ); fh_it.is_valid(); fh_it++ )
        {
            p0 = mesh.point( mesh.from_vertex_handle( *fh_it ) );
            p1 = mesh.point( mesh.to_vertex_handle( *fh_it ) );
            l = norm( P( p0 ) - P( p1 ) );
            mesh.data( mesh.edge_handle( *fh_it ) ).length = l;
            perimeter += l;
        }
        mesh.data( *f_it ).area = -Area_Face( mesh, *f_it );
        double comp =  4 * M_PI * mesh.data( *f_it ).area / pow( perimeter, 2 );
        if ( comp >= 1 ) std::cout<<"\ne Error in Find_Compactness: comp="<<comp;
        compactness += comp;
    }
    compactness /= mesh.n_faces(); // average over all faces;
    compactness *= 100; // in percents
}

double Area_Polygon (const Polygon& poly)
{
    double area = 0;
    if ( poly.size() < 3 ) return 0;
    for ( int i = 0; i<poly.size(); i++)
    {
        double x = poly[i].x(), y = poly[i].y();
        int i_n = (i + 1) % poly.size();
        double x_n = poly[i_n].x(), y_n = poly[i_n].y();
        area += (y_n - y) * (x + x_n);
    }
    return area/2;
}

double Area_Polygon (const std::vector<Point2d>& polygon)
{
    double area = 0;
    if ( polygon.size() < 3 ) return 0;
    for ( int i = 0; i<polygon.size(); i++)
    {
        double x = polygon[i].x, y = polygon[i].y;
        int i_n = (i + 1) % polygon.size();
        double x_n = polygon[i_n].x, y_n = polygon[i_n].y;
        area += (y_n - y) * (x + x_n);
    }
    return area/2;
}

bool Convert (Mesh& mesh, MFH fh, Polygon& poly)
{
    poly.clear();
    for ( auto fv_it = mesh.fv_iter( fh ); fv_it.is_valid(); fv_it++ )
        poly.push_back( P2( mesh.point( *fv_it )[0]+0.1, mesh.point( *fv_it )[1]+0.1 ) );
    return true;
}

bool Convert (Mesh& mesh, MFH fh, std::vector<Point2d>& polygon)
{
    polygon.clear();
    polygon.reserve( mesh.valence( fh ) );
    for ( auto fv_it = mesh.fv_iter( fh ); fv_it.is_valid(); ++fv_it )
        polygon.push_back( P( mesh.point( *fv_it ) ) );
    return true;
}

void Scale_Mesh (Mesh& mesh, double scale)
{
    for ( auto it = mesh.vertices_begin(); it != mesh.vertices_end(); ++it ) mesh.point( *it ) *= scale;
}

bool Enlarge_Mesh (Mesh& mesh, Size size)
{
    for ( auto it = mesh.vertices_begin(); it != mesh.vertices_end(); ++it )
    {
        if ( mesh.point( *it )[0] + 1 >= size.width ) mesh.point( *it )[0]++;
        if ( mesh.point( *it )[1] + 1 >= size.height ) mesh.point( *it )[1]++;
    }
    return true;
}

void Draw_Edges (Mesh& mesh, Scalar color, Mat& image)
{
    MP p0, p1;
    for ( auto it = mesh.edges_begin(); it != mesh.edges_end(); it++ )
    {
        p0 = mesh.point( mesh.from_vertex_handle( mesh.halfedge_handle( *it, 0 ) ) );
        p1 = mesh.point( mesh.from_vertex_handle( mesh.halfedge_handle( *it, 1 ) ) );
        line( image, P( p0 ), P( p1 ), color, 1 ); //CV_AA );
    }
}

void Find_Bounding_Box (Mesh& mesh, MFH fh, std::pair<Point2d,Point2d>& box)
{
    Point2d p;
    box.first = Point2d( 1e+8, 1e+8 );
    box.second = Point2d( -1e+8, -1e+8 );
    for ( auto fv_it = mesh.fv_iter( fh ); fv_it.is_valid(); ++fv_it )
    {
        p = P( mesh.point( *fv_it ) );
        if ( box.first.x > p.x ) box.first.x = p.x;
        if ( box.first.y > p.y ) box.first.y = p.y;
        if ( box.second.x < p.x ) box.second.x = p.x;
        if ( box.second.y < p.y ) box.second.y = p.y;
    }
}

double Area_Face (Mesh& mesh, MFH fh)
{
    //Polygon poly;
    std::vector<Point2d> polygon;
    polygon.reserve( mesh.valence( fh ) );
    Convert( mesh, fh, polygon );
    //if ( ! poly.is_simple() ) std::cout<<"\nNon-simple face"<<mesh.data( fh ).label<<" "<<poly;
    return Area_Polygon( polygon );
}

double Area_Faces (Mesh& mesh)
{
    double area = 0;
    for ( auto f_it = mesh.faces_begin(); f_it != mesh.faces_end(); f_it++ )
    {
        mesh.data( *f_it ).area = -Area_Face( mesh, *f_it );
        area += mesh.data( *f_it ).area;
    }
    return area;
}

//
bool Type_Labels (Mesh& mesh, int scale, int shift, Scalar color, int thickness, Mat& image)
{
    Point2d c( 0, 0 ), s( -shift, shift );
    for ( auto f_it = mesh.faces_begin(); f_it != mesh.faces_end(); f_it++ )
    {
        c = Point2d( 0, 0 );
        for ( auto fv_it = mesh.fv_iter( *f_it ); fv_it.is_valid(); fv_it++ ) c += P( mesh.point( *fv_it ) );
        c *= 1.0 / mesh.valence( *f_it );
        putText( image, std::to_string( mesh.data( *f_it ).label ), (s + c) * scale, FONT_HERSHEY_PLAIN, thickness, color );
    }
    return true;
}//

bool Draw_Mesh (Mesh& mesh, int scale, int shift, Scalar edge_color, Scalar vertex_color, int thickness, Mat& image)
{
    MHH hh;
    Point2d p0, p1, s( shift, shift );
    resize( image, image, image.size() * scale );
    for ( auto e_it = mesh.edges_begin(); e_it != mesh.edges_end(); e_it++ )
    {
        hh = mesh.halfedge_handle( *e_it, 0 );
        p0 = ( P( mesh.point( mesh.from_vertex_handle( hh ) ) ) + s ) * scale;
        p1 = ( P( mesh.point( mesh.to_vertex_handle( hh ) ) ) + s ) * scale;
        line( image, p0, p1, edge_color, thickness, LINE_AA );
    }
    for ( auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++ )
        if ( mesh.valence( *v_it ) > 2 )
            circle( image, ( P( mesh.point( *v_it ) ) + s ) * scale, 2 * thickness, vertex_color, -1 );
        else circle( image, ( P( mesh.point( *v_it ) ) + s ) * scale, 2 * thickness, edge_color, -1 );
    //Type_Labels( mesh, scale, shift, vertex_color, thickness, image );
    return true;
}

bool Draw_Reconstruction (Mesh& mesh, int scale, Mat& image)
{
    resize( image, image, image.size() * scale );
    for ( auto f_it = mesh.faces_begin(); f_it != mesh.faces_end(); f_it++ )
    {
        int ind = 0;
        std::vector<Point> contour( mesh.valence( *f_it ) );
        int num_vertices = mesh.valence( *f_it );
        for ( auto fv_it = mesh.fv_iter( *f_it ); fv_it.is_valid(); fv_it++ )
            contour[ ind++ ] = Point2i( int( mesh.point( *fv_it )[0] * scale ), int( mesh.point( *fv_it )[1] * scale ) ) ;
        const Point* points[1] = { &contour[0] }; // vector of 1 pointer to the vector contour containing num_vertices
        Scalar c = Scalar( int( mesh.data( *f_it ).mean ) );
        if ( image.channels() == 3 ) c = Scalar(  mesh.data( *f_it ).color );
        fillPoly( image, points, &num_vertices, 1, c );
    }
    return true;
}

bool Angles_Ok (Mesh& mesh, MHH h, Point2d p0, Point2d p1)
{
    MHH h_prev = mesh.prev_halfedge_handle( h );
    double sine, cosine;
    Find_Sine_Cosine( P( mesh.point( mesh.from_vertex_handle( h_prev ) ) ), p0, p1, sine, cosine );
    if ( abs( sine ) < sin( Min_Angle ) && cosine > cos( Min_Angle ) ) return false;
    MHH h_next = mesh.next_halfedge_handle( mesh.opposite_halfedge_handle( h ) );
    Find_Sine_Cosine( P( mesh.point( mesh.to_vertex_handle( h_next ) ) ), p0, p1, sine, cosine );
    if ( abs( sine ) < sin( Min_Angle ) && cosine > cos( Min_Angle ) ) return false;
    return true;
}

bool Intersections (Mesh& mesh, MHH h_begin, MHH h_end, double a, double b, double c, Point2d start_point, Point2d end_point, bool debug)
{
    if ( h_begin == h_end ) return true; // single edge complements a chain
    MHH h = mesh.next_halfedge_handle( h_begin );
    Point2d p0 = P( mesh.point( mesh.from_vertex_handle( h ) ) ), p1;
    Segment straight( P2d( start_point ), P2d( end_point ) );
    while ( h != h_end )
    {
        p1 = P( mesh.point( mesh.to_vertex_handle( h ) ) );
        double prod = ( a * p0.x + b * p0.y + c ) * ( a * p1.x + b * p1.y + c );
        if ( debug ) std::cout<<" p0="<<p0<<"p1="<<p1<<prod;
        if ( prod <= 0 ) // on different sides
        {
            Segment edge( P2d( p0 ), P2d( p1 ) );
            if ( do_intersect( straight, edge ) ) return true; // intersection exists
            if ( debug ) std::cout<<"%";
        }
        h = mesh.next_halfedge_handle( h );
        p0 = p1;
    }
    return false;
}

bool Find_Elbow (Mesh& mesh, const std::vector<MHH>& chain, double color_offset, int& elbow)
{
    bool debug = false;
    int check_face = -1;
    Point2d check_point( 259, 209 );
    double a,b,c,d, max_offset = 0;
    MHH h0 = chain[0];
    MHH h1 = chain[ (int) chain.size()-1 ];
    MVH vh0 = mesh.from_vertex_handle( h0 );
    MVH vh1 = mesh.to_vertex_handle( h1 );
    Point2d p0 = P( mesh.point( vh0 ) );
    Point2d p1 = P( mesh.point( vh1 ) );
    Points_to_Line_Coef( p0, p1, a, b, c );
    //if ( p0 == check_point ) std::cout<<"\np0="<<p0<<" p1="<<p1<<" a="<<a<<" b="<<b<<" c="<<c;
    double s = chain.size();
    elbow = -1;
    for ( int i = 0; i+1 < s; i++ )
    {
        Distance_Point_to_Line( P( mesh.point( mesh.to_vertex_handle( chain[i] ) ) ), a, b, c, d );
        if ( max_offset < d or ( max_offset == d and fabs( i - 0.5 * s ) < abs( elbow - 0.5 * s ) ) ) { max_offset = d; elbow = i; }
    }
    if ( elbow < 0 ) return false; // no elbows found at all
    
    // Check that no small angles appear
    MHH h0back = mesh.opposite_halfedge_handle( h0 );
    MHH h1back = mesh.opposite_halfedge_handle( h1 );
    bool angles_ok = Angles_Ok( mesh, h0, p0, p1 ) and Angles_Ok( mesh, h1back, p1, p0 );
    if ( max_offset < 1 and angles_ok ) return false; // no big elbows found
    
    // Check if the colors of adjacent faces have a large difference
    MFH fh0 = mesh.face_handle( h0 );
    MFH fh1 = mesh.opposite_face_handle( h0 );
    double delta_color = color_offset / max_offset; // main formula relating distance and color
    double dif = fabs( mesh.data( fh0 ).mean - mesh.data( fh1 ).mean );
    if ( dif >= delta_color ) return true; // too large difference in color
    
    // Check if potential interesections can forbid straightening
    if ( ! angles_ok ) return true; // large enough angles are needed further
    MVH vh = mesh.to_vertex_handle( chain[ elbow ] );
    if ( mesh.data( fh0 ).label == check_face or mesh.data( fh1 ).label == check_face )
    {
        std::cout<<"\n f0="<<mesh.data( fh0 ).label<<" f1="<<mesh.data( fh1 ).label<<" a="<<a<<" b="<<b<<" c="<<c;
        Print_Vertex( mesh, vh0 ); Print_Vertex( mesh, vh ); Print_Vertex( mesh, vh1 );
    }
    if ( Intersections( mesh, mesh.next_halfedge_handle( h1 ), mesh.prev_halfedge_handle( h0 ), a, b, c, p0, p1, debug ) )
        return true;
    if ( Intersections( mesh, mesh.next_halfedge_handle( h0back ), mesh.prev_halfedge_handle( h1back ), a, b, c, p0, p1, debug ) )
        return true;
    return false; // small angle found
}

bool Collapse_Chain (Mesh& mesh, const std::vector<MHH>& chain)
{
    //std::cout<<"\nCollapse: "; Print_Vertex( mesh, mesh.from_vertex_handle( chain[0] ) );
    //for ( auto h : chain ) Print_Vertex( mesh, mesh.to_vertex_handle( h ) );
    mesh.request_vertex_status();
    mesh.request_edge_status();
    mesh.request_face_status();
    for ( int i = 1; i < chain.size(); i++ ) mesh.collapse( chain[i] );
    return true;
}

bool Straighten_Chain (Mesh& mesh, double color_offset, const std::vector<MHH>& chain)
{
    if ( chain.size() == 1 ) return true;
    int elbow = -1;
    if ( ! Find_Elbow( mesh, chain, color_offset, elbow ) ) { Collapse_Chain( mesh, chain ); return true; }
    std::vector<MHH> subchain;
    subchain.assign( chain.begin(), chain.begin() + elbow + 1 );
    if ( chain.size() <= subchain.size() ) std::cout<<"\nError in Straighten_Chain: subchain not shorter"; //else std::cout<<".";
    Straighten_Chain( mesh, color_offset, subchain );
    subchain.assign( chain.begin() + elbow + 1, chain.end() );
    if ( chain.size() <= subchain.size() ) std::cout<<"\nError in Straighten_Chain: subchain not shorter"; //else std::cout<<".";
    Straighten_Chain( mesh, color_offset, subchain );
    return true;
}

bool Analyse_Chain (Mesh& mesh, double color_offset, const std::vector<MHH>& chain, std::vector< std::vector<MHH> >& chains)
{
    if ( chain.size() == 1 ) return true;
    int elbow = -1;
    if ( ! Find_Elbow( mesh, chain, color_offset, elbow ) ) { Collapse_Chain( mesh, chain ); return true; }
    if ( elbow < 0 ) { chains.push_back( chain ); return true; } // remove this chain
    std::vector<MHH> subchain;
    subchain.assign( chain.begin(), chain.begin() + elbow + 1 );
    if ( chain.size() <= subchain.size() ) std::cout<<"\nError in Straighten_Chain: subchain not shorter"; //else std::cout<<".";
    Analyse_Chain( mesh, color_offset, subchain, chains );
    subchain.assign( chain.begin() + elbow + 1, chain.end() );
    if ( chain.size() <= subchain.size() ) std::cout<<"\nError in Straighten_Chain: subchain not shorter"; //else std::cout<<".";
    Analyse_Chain( mesh, color_offset, subchain, chains );
    return true;
}

bool Find_Chain (Mesh& mesh, MHH hh, std::vector<MHH>& chain)
{
    chain.assign( 1, hh );
    MHH h = mesh.next_halfedge_handle( hh );
    while ( mesh.valence( mesh.from_vertex_handle( h ) ) == 2 )
    {
        chain.push_back( h );
        h = mesh.next_halfedge_handle( h );
    }
    // Check a specific face for debugging
    int check_face = -1;
    if ( check_face >= 0 )
    {
        MFH fh0 = mesh.face_handle( hh );
        MFH fh1 = mesh.opposite_face_handle( hh );
        if ( mesh.data( fh0 ).label == check_face or mesh.data( fh1 ).label == check_face )
        {
            std::cout<<"\n f0="<<mesh.data( fh0 ).label<<" f1="<<mesh.data( fh1 ).label;
            MVH vh = mesh.from_vertex_handle( chain[0] );
            std::cout<<"\nChain from"; Print_Vertex( mesh, vh );
            for ( auto h : chain ) Print_Vertex( mesh, mesh.to_vertex_handle( h ) );
        }
    }
    return true;
}

bool Find_Edge_in_Chain (MHH hh, const std::vector<MHH>& chain)
{
    for ( auto h : chain ) if ( hh == h ) return true;
    return false;
}

bool Chain_Removable (Mesh& mesh, const std::vector<MHH>& chain)
{
    MFH fh0 = mesh.face_handle( chain[0] );
    MFH fh1 = mesh.opposite_face_handle( chain[0] );
    for ( auto h_it = mesh.fh_iter( fh0 ); h_it.is_valid(); h_it++ )
        if ( mesh.opposite_face_handle( *h_it ) == fh1 and ! Find_Edge_in_Chain( *h_it, chain ) )
            return false;
    //std::cout<<"\nFaces will be merged: f0="<<mesh.data( fh0 ).label<<" f1="<<mesh.data( fh1 ).label;
    return true;
}

bool Remove_Trivial_Vertex (Mesh& mesh, MVH vh)
{
    if ( mesh.valence( vh ) != 2 ) return false;
    auto h0 = mesh.voh_iter( vh ), h1 = h0;
    h1++;
    if ( ! Collinear( P( mesh.point( mesh.to_vertex_handle( *h0 ) ) ), P( mesh.point( vh) ), P( mesh.point( mesh.to_vertex_handle( *h1 ) ) ) ) ) return false;
    mesh.request_vertex_status();
    mesh.request_edge_status();
    mesh.request_face_status();
    mesh.collapse( *h0 );
    return true;
}

bool Remove_Trivial_Vertices (Mesh& mesh)
{
    for ( auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++ )
        if ( ! mesh.is_valid_handle( *v_it ) or mesh.valence( *v_it ) != 2 ) continue;
        else Remove_Trivial_Vertex( mesh, *v_it );
    mesh.garbage_collection( true, true, true );
    return true;
}

bool Remove_Chain (Mesh& mesh, double min_color_dif, const std::vector<MHH>& chain, std::vector< std::vector<MHH> >& chains)
{
    MFH fh0 = mesh.face_handle( chain[0] );
    MFH fh1 = mesh.opposite_face_handle( chain[0] );
    double dif = fabs( mesh.data( fh0 ).mean - mesh.data( fh1 ).mean );
    std::cout<<"\nd="<<dif;
    if ( dif >= min_color_dif and mesh.data( fh0 ).area > 9 and mesh.data( fh1 ).area > 9 ) return false; // too large difference
    if ( ! Chain_Removable( mesh, chain ) ) return false;
    if ( chain.size() + 1 != mesh.valence( fh0 ) and chain.size() + 1 != mesh.valence( fh1 ) ) return false;
    std::cout<<" f0="<<mesh.data( fh0 ).label<<" f1="<<mesh.data( fh1 ).label<<" a0="<<mesh.data( fh0 ).area<<" a1="<<mesh.data( fh1 ).area;
    /* MVH vh0 = mesh.from_vertex_handle( chain[ 0 ] );
    MVH vh1 = mesh.to_vertex_handle( chain[ (int)chain.size() - 1 ] );
    Collapse_Chain( mesh, chain );
    if ( mesh.valence( vh0 ) == 2 and ! Remove_Trivial_Vertex( mesh, vh0 ) );
    if ( mesh.valence( vh1 ) == 2 and ! Remove_Trivial_Vertex( mesh, vh1 ) );
    */
    chains.push_back( chain );
    return true;
}

void Remove_Edge (Mesh& mesh, MHH h)
{
    if ( ! mesh.is_valid_handle( h ) ) return;
    MEH eh = mesh.edge_handle( h );
    if ( mesh.status( eh ).deleted() ) return;
    if ( ! mesh.is_simple_link(eh) ) return;
    mesh.request_edge_status();
    mesh.request_face_status();
    mesh.remove_edge( eh );
}

bool Straighten_Chains (Mesh& mesh, double color_offset, double min_color_dif)
{
    MEH eh;
    std::vector<MHH> chain, subchain;
    std::vector< std::vector<MHH> > chains;
    for ( auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++ )
    {
        if ( ! mesh.is_valid_handle( *v_it ) ) continue;
        if ( mesh.valence( *v_it ) == 2 ) continue; // start chains only from vertices of deg>2
        for ( auto h_it = mesh.voh_iter( *v_it ); h_it.is_valid(); h_it++ )
        {
            if ( ! mesh.is_valid_handle( * h_it ) ) continue;
            eh = mesh.edge_handle( *h_it );
            if ( mesh.is_boundary( eh ) or mesh.data( eh ).chain_found ) continue;
            mesh.data( eh ).chain_found = true;
            //Print_Edge( mesh, eh );
            Find_Chain( mesh, *h_it, chain );
            mesh.data( mesh.edge_handle( chain[ (int)chain.size()-1 ] ) ).chain_found = true;
            //if ( ! Remove_Chain( mesh, min_color_dif, chain, chains ) )
                Straighten_Chain( mesh, color_offset, chain );
            //Analyse_Chain( mesh, color_offset, chain, chains );
            //std::cout<<"*";
        }
    }
    mesh.garbage_collection( true, true, true );
    //std::cout<<"\nChains were straightened";
    /*
    for ( auto chain : chains ) Collapse_Chain( mesh, chain );
    mesh.garbage_collection( true, true, true );
    std::cout<<"\nChains were removed";
    */
    return true;
}

bool Find_Single_Edges (Mesh& mesh, double min_color_dif, std::vector<MEH>& edges)
{
    MHH hh;
    MFH fh0, fh1;
    MVH vh0, vh1;
    int ind = 0;
    double dif;
    std::map<int, MEH> edge_handles;
    std::vector<Index_Value> sorted_edges;
    for ( auto e_it = mesh.edges_begin(); e_it != mesh.edges_end(); e_it++ )
    {
        if ( mesh.is_boundary( *e_it ) ) continue;
        hh = mesh.halfedge_handle( *e_it, 0 );
        vh0 = mesh.from_vertex_handle( hh );
        if ( mesh.valence( vh0 ) < 3 ) continue;
        vh1 = mesh.to_vertex_handle( hh );
        if ( mesh.valence( vh1 ) < 3 ) continue;
        if ( ! mesh.is_valid_handle( *e_it ) ) continue;
        if ( ! mesh.is_simple_link( *e_it ) ) continue;
        fh0 = mesh.face_handle( hh );
        fh1 = mesh.opposite_face_handle( hh );
        dif = fabs( mesh.data( fh0 ).mean - mesh.data( fh1 ).mean );
        if ( dif >= min_color_dif and mesh.data( fh0 ).area > 9 and mesh.data( fh1 ).area > 9 ) continue; // too large difference
        //std::cout<<"\nd="<<dif<<" f0="<<mesh.data( fh0 ).label<<" f1="<<mesh.data( fh1 ).label<<" a0="<<mesh.data( fh0 ).area<<" a1="<<mesh.data( fh1 ).area<<" erase: "; Print_Vertex( mesh, vh0 ); Print_Vertex( mesh, vh1 );
        edge_handles.insert( make_pair( ind, *e_it ) );
        sorted_edges.push_back( Index_Value( ind, dif) );
        ind++;
    }
    sort( sorted_edges.begin(), sorted_edges.end(), Increasing_Values );
    for ( size_t e = 0; e < sorted_edges.size(); e++ )
    {
        //std::cout<<"\nd="<<sorted_edges[e].value; Print_Edge( mesh, edge_handles[ sorted_edges[e].index ] );
        edges.push_back( edge_handles[ sorted_edges[e].index ] );
    }
    return true;
}

bool Find_Chain (Mesh& mesh, MFH fh, std::vector<MHH>& chain)
{
    auto h_it = mesh.fh_iter( fh );
    if ( mesh.is_boundary( *h_it ) ) h_it++;
    if ( mesh.is_boundary( *h_it ) ) h_it++;
    if ( mesh.is_boundary( *h_it ) ) { std::cout<<"\nError in Find_Chain: "; Print_Face( mesh, fh ); }
    chain.assign( 1, *h_it );
    MFH f_opposite = mesh.opposite_face_handle( *h_it );
    MHH h = mesh.prev_halfedge_handle( *h_it );
    while ( mesh.opposite_face_handle( h ) == f_opposite )
    {
        chain.insert( chain.begin(), h );
        h = mesh.prev_halfedge_handle( h );
    }
    h = mesh.next_halfedge_handle( *h_it );
    while ( mesh.opposite_face_handle( h ) == f_opposite )
    {
        chain.push_back( h );
        h = mesh.next_halfedge_handle( h );
    }
    return true;
}

bool Remove_Small_Faces (Mesh& mesh, double min_area_face)
{
    std::vector<MHH> chain;
    for ( auto f_it = mesh.faces_begin(); f_it != mesh.faces_end(); f_it++ )
        if ( mesh.data( *f_it ).area < min_area_face )
        {
            //Print_Face( mesh, *f_it );
            Find_Chain( mesh, *f_it, chain );
            if ( chain.size() + 1 != mesh.valence(  *f_it ) ) continue; // removing this chain won't collapse the face
            //std::cout<<"\nChain:"; for ( auto h : chain ) Print_Edge( mesh, h );
            Collapse_Chain( mesh, chain );
        }
    mesh.garbage_collection( true, true, true );
    Remove_Trivial_Vertices( mesh );
    return true;
}

bool Merge_Faces (Mesh& mesh, double min_color_dif)
{
    MFH fh, fh0, fh1;
    std::vector<MEH> edges;
    double weight, area;
    Find_Single_Edges( mesh, min_color_dif, edges );
    for ( auto edge : edges )
    {
        if ( ! mesh.is_valid_handle( edge ) ) continue;
        //std::cout<<"\nErase:"; Print_Edge( mesh, edge );
        mesh.request_edge_status();
        mesh.request_face_status();
        if ( ! mesh.is_simple_link( edge ) ) continue;
        fh0 = mesh.face_handle( mesh.halfedge_handle( edge, 0 ) );
        fh1 = mesh.opposite_face_handle( mesh.halfedge_handle( edge, 0 ) );
        area = mesh.data( fh0 ).area + mesh.data( fh1 ).area;
        weight = mesh.data( fh0 ).mean * mesh.data( fh0 ).area + mesh.data( fh1 ).mean * mesh.data( fh1 ).area;
        fh = mesh.remove_edge( edge );
        Polygon polygon;
        Convert( mesh, fh, polygon );
        if ( ! polygon.is_simple() ) { mesh.reinsert_edge( edge ); continue; }
        mesh.data( fh ).area = area;
        mesh.data( fh ).mean = 1.0 * weight / area;
    }
    mesh.garbage_collection( true, true, true );
    //Remove_Trivial_Vertices( mesh );
    return true;
}

void Update_Labels (Mesh& mesh, MFH fh, Mat_<int>& labels)
{
    std::pair<Point2d,Point2d> box;
    Find_Bounding_Box( mesh, fh, box );
    Polygon polygon;
    Convert( mesh, fh, polygon );
    //std::cout<<"\nl"<<mesh.data( fh ).label;
    //for ( int i = 0; i < polygon.size(); i++ ) std::cout<<" ["<<polygon[i].x()<<","<<polygon[i].y()<<"]";
    for ( int j = int( box.first.x ); j <= int( box.second.x ); j++ ) // columns
        for ( int i = int( box.first.y ); i <= int( box.second.y ); i++ ) // rows
            if ( Is_Nonstrictly_Inside_Polygon( polygon, P2( j+0.5, i+0.5 ) ) )
                labels( i, j ) = mesh.data( fh ).label;
}

bool Mesh_to_Superpixels (Mesh& mesh, Point sizes, std::vector< std::vector<Point> >& superpixels, Mat_<int>& indices, int& index)
{
    index = 0;
    std::vector<Point> pixels;
    pixels.reserve( sizes.x * sizes.y / 10 );
    std::pair<Point2d,Point2d> box;
    indices = Mat_<int>( sizes.y, sizes.x ); //  rows x cols
    for ( int i = 0; i < indices.rows; i++ )
        for ( int j = 0; j < indices.cols; j++ ) indices( i, j ) = -1;
    //std::cout<<" sizes="<<sizes<<" indices.rows="<<indices.rows<<" indices.cols="<<indices.cols;
    for ( auto f_it = mesh.faces_begin(); f_it != mesh.faces_end(); f_it++ )
    {
        Find_Bounding_Box( mesh, *f_it, box );
        Polygon polygon;
        Convert( mesh, *f_it, polygon );
        if ( ! polygon.is_simple() ) { std::cout<<"\nNon-simple polygon:"; Print_Face_Vertices( mesh, *f_it ); }
        pixels.clear();
        for ( int j = int( box.first.x ); j <= int( box.second.x ); j++ ) // columns
            for ( int i = int( box.first.y ); i <= int( box.second.y ); i++ ) // rows
                if ( Is_Nonstrictly_Inside_Polygon( polygon, P2( j+0.5, i+0.5 ) ) )
                    pixels.push_back( Point( j, i ) );
        if ( pixels.size() == 0 ) // a new superpixel was found
        { std::cout<<"\nEmpty superpixel: "<<polygon; continue; }
        for ( int k = 0; k < pixels.size(); k++ ) indices( pixels[k] ) = index;
        mesh.data( *f_it ).label = index; // labels of faces are updated
        index++;
    }
    if ( index != mesh.n_faces() ) std::cout<<"\nError in Mesh_to_Superpixels: index="<<index<<" f="<<mesh.n_faces();
    //
    std::vector<Point> empty;
    empty.reserve( sizes.x * sizes.y / 10 );
    superpixels.assign( index, empty );
    for ( int j = 0; j < indices.cols; j++ )
        for ( int i = 0; i < indices.rows; i++ )
        {
            if ( indices( i, j ) < 0 || indices( i, j ) >= index ) { std::cout<<" i="<<indices( i, j )<<Point( j, i ); return false; }
            else superpixels[ indices( i, j ) ].push_back( Point( j, i ) );
        }
    return true;
}

bool Mesh_to_Segmentation (Mesh& mesh, Point sizes, iSegmentation& segmentation)
{
    Mat bmat( sizes, CV_8UC1, Scalar( 0 ) );
    Draw_Edges( mesh, Scalar( 255 ), bmat );
    adaptiveThreshold( bmat, bmat, 255, ADAPTIVE_THRESH_MEAN_C, THRESH_BINARY, 3, 0 );
    //Write_Image( "/Users/kurlin/C++/RIS/Output/bmat.png", bmat );
    segmentation.num_faces = (int)mesh.n_faces();
    segmentation.num_edges = (int)mesh.n_edges();
    segmentation.num_vertices = (int)mesh.n_vertices();
    segmentation.boundary_mask = Mat_<bool>( sizes.y, sizes.x );
    segmentation.boundary.clear();
    for ( int i = 0; i < bmat.rows; i++ )
        for ( int j = 0; j < bmat.cols; j++ )
            if ( bmat.at<uchar>( i, j ) == 255 ) segmentation.boundary_mask( i, j ) = true;
            else segmentation.boundary_mask( i, j ) = false;
    if ( ! Mesh_to_Superpixels( mesh, sizes, segmentation.superpixels, segmentation.indices, segmentation.num_faces ) ) return false;
    return true;
}

bool Off_to_Segmentation (String path, String name, Point sizes, iSegmentation& segmentation)
{
    Mesh mesh;
    String file_name = path + name + ".off";
    if ( ! OpenMesh::IO::read_mesh( mesh, file_name ) ) { std::cout<<"\nError reading "<<file_name; return false; }
    segmentation.num_edges = (int)mesh.n_edges();
    Mesh_to_Segmentation( mesh, sizes, segmentation );
    return true;
}

bool Pixel_Based_Segmentation (PixelSegmentation& s, Point sizes, iSegmentation& segmentation )
{
    std::vector<Point> empty;
    empty.reserve( sizes.x * sizes.y / 100 );
    segmentation.superpixels.assign( s.segments.size(), empty );
    segmentation.indices = Mat_<int>( sizes.y, sizes.x );
    for ( int i = 0; i < s.segments.size(); i++ )
    {
        for ( int j = 0; j < s.segments[i].xs.size(); j++ )
        {
            Point p( s.segments[i].xs[j], s.segments[i].ys[j] );
            segmentation.superpixels[ s.segments[i].label - 1 ].push_back( p );
            segmentation.indices( p ) = s.segments[i].label - 1;
        }
    }
    segmentation.boundary_mask = Mat_<bool>( sizes.y, sizes.x, false );
    Mat b = s.get_boundary_pixels();
    Mat boundaries( sizes, CV_8UC1, Scalar( 0 ) );
    segmentation.boundary.reserve( sizes.x * sizes.y / 10 );
    for ( int i = 0; i < sizes.y; i++ )
        for ( int j = 0; j < sizes.x; j++ )
            if ( b.at<uchar>( i, j ) == 1 )
            {
                segmentation.boundary.push_back( Point( j, i ) );
                segmentation.boundary_mask( Point( j, i ) ) = true;
                boundaries.at<uchar>( i, j ) = 255;
            }
    //Write_Image( "/Users/kurlin/C++/RIS/Output/boundaries.png", boundaries );
    return true;
}

bool Find_Benchmarks (Mesh& mesh, String name, Point sizes, Benchmark& average, std::ostream& fout)
{
    int image_index;
    iSegmentation segmentation, ground_truth;
    Mesh_to_Segmentation( mesh, sizes, segmentation );
    Benchmark benchmark(1);
    benchmark.num_faces = segmentation.num_faces;
    benchmark.num_edges = segmentation.num_edges;
    benchmark.num_vertices = segmentation.num_vertices;
    image_index = 0;
    while ( Read_BSD_Human( name + "_" + Str( ++image_index ) + ".dat", ground_truth ) )
    {
        Boundary_Recall( ground_truth, segmentation, benchmark );
        Undersegmentation_Errors( ground_truth, segmentation, benchmark );
    }//
    average.Add( benchmark );
    benchmark.Print( std::cout, fout );
    return true;
}

bool Find_Means (const Mat& image, const Mat_<int>& labels, Mesh& mesh, std::vector<double>& means)
{
    if ( image.channels() != 1 ) return false;
    //
    int max_label = 0;
    for ( auto f_it = mesh.faces_begin(); f_it != mesh.faces_end(); f_it++ )
        if ( max_label < mesh.data( *f_it ).label ) max_label = mesh.data( *f_it ).label;
    means.assign( max_label + 1, 0 ); // starting from label 0
    std::vector<int> areas( max_label + 1, 0 );
    // means.assign( mesh.n_faces(), 0 );
    //std::vector<int> areas( mesh.n_faces(), 0 );
    for ( int i = 0; i < labels.rows; i++ )
        for ( int j = 0; j < labels.cols; j++ )
        {
            means[ labels( i, j ) ] += image.at<uchar>( i, j );
            areas[ labels( i, j ) ]++;
            //if ( labels( i, j ) == 123 ) std::cout<<" c123"<<Point(j,i)<<"="<<image.at<uchar>( i, j );
        }
    for ( auto f_it = mesh.faces_begin(); f_it != mesh.faces_end(); f_it++ )
    {
        if ( mesh.data( *f_it ).label >= means.size() ) std::cout<<"\nError in Find_Means: l="<<mesh.data( *f_it ).label;
        if ( areas[ mesh.data( *f_it ).label ] > 0 ) means[ mesh.data( *f_it ).label ] /= areas[ mesh.data( *f_it ).label ];
        //std::cout<<" m"<< mesh.data( *f_it ).label<<"="<< means[ mesh.data( *f_it ).label ]<<" a="<<areas[ mesh.data( *f_it ).label ];
    }
    return true;
}

bool Find_Colors (const Mat& image, const Mat_<int>& labels, Mesh& mesh, std::vector<Vec3b>& colors)
{
    if ( image.channels() != 3 ) return false;
    colors.assign( mesh.n_faces(), 0 ); // starting from label 0
    std::vector<int> areas( mesh.n_faces(), 0 );
    for ( int i = 0; i < labels.rows; i++ )
        for ( int j = 0; j < labels.cols; j++ )
        {
            colors[ labels( i, j ) ] += image.at<Vec3b>( i, j );
            areas[ labels( i, j ) ]++;
        }
    for ( auto f_it = mesh.faces_begin(); f_it != mesh.faces_end(); f_it++ )
    {
        if ( mesh.data( *f_it ).label >= colors.size() ) std::cout<<"\nError in Find_Means: ="<<mesh.data( *f_it ).label;
        if ( areas[ mesh.data( *f_it ).label ] > 0 ) colors[ mesh.data( *f_it ).label ] /= areas[ mesh.data( *f_it ).label ];
    //std::cout<<" c"<< mesh.data( *f_it ).label<<"="<< colors[ mesh.data( *f_it ).label ]<<" a="<<areas[ mesh.data( *f_it ).label ];
    }
    return true;
}

void Save_Means (const std::vector<double>& means, Mesh& mesh)
{
    for ( auto f_it = mesh.faces_begin(); f_it != mesh.faces_end(); f_it++ )
        mesh.data( *f_it ).mean = means[ mesh.data( *f_it ).label ];
}

void Save_Colors (const std::vector<Vec3b>& colors, Mesh& mesh)
{
    for ( auto f_it = mesh.faces_begin(); f_it != mesh.faces_end(); f_it++ )
        mesh.data( *f_it ).color = colors[ mesh.data( *f_it ).label ];
}

void Find_Means (const Mat& image, const Mat_<int>& labels, Mesh& mesh)
{
    std::vector<double> means;
    Find_Means( image, labels, mesh, means );
    Save_Means( means, mesh );
}

void Find_Mesh_Colors (const Mat& image, const Mat_<int>& labels, Mesh& mesh)
{
    if ( image.channels() == 1 )
    {
        std::vector<double> means;
        Find_Means( image, labels, mesh, means );
        Save_Means( means, mesh );
    }
    if ( image.channels() == 3 )
    {
        std::vector<Vec3b> colors;
        Find_Colors( image, labels, mesh, colors );
        Save_Colors( colors, mesh );
    }
}

double Reconstruction_Error (const Mat& image, const Mat_<int>& labels, size_t num_superpixels)
{
    double error = 0;
    std::vector<double> means( num_superpixels, 0 ); // starting from label 0
    std::vector<int> areas( num_superpixels, 0 );
    for ( int i = 0; i < labels.rows; i++ )
        for ( int j = 0; j < labels.cols; j++ )
        {
            //if ( image.channels() == 1 )
            means[ labels( i, j ) ] += image.at<uchar>( i, j );
            areas[ labels( i, j ) ]++;
        }
    for ( int i = 0; i < means.size(); i++ ) means[i] /= areas[i];
    for ( int i = 0; i < labels.rows; i++ )
        for ( int j = 0; j < labels.cols; j++ )
        {
            if ( means[ labels( i, j ) ] > 255 ) std::cout<<"\nError in Reconstruction_Error: m"<<labels( i, j )<<"="<<means[ labels( i, j ) ];
            error += pow( image.at<uchar>( i, j ) - means[ labels( i, j ) ] , 2 );
        }
    error /= image.rows * image.cols;
    error = sqrt( error ) * 100.0 / 255;
    return error;
}

void Read_Means (Mesh& mesh, std::vector<double>& means)
{
    means.resize( mesh.n_faces() );
    for ( auto f_it = mesh.faces_begin(); f_it != mesh.faces_end(); f_it++ )
        means[ mesh.data( *f_it ).label ] = mesh.data( *f_it ).mean;
}

double Reconstruction_Error (const Mat& image, const Mat_<int>& labels, Mesh& mesh)
{
    double error = 0;
    std::vector<double> means;
    Read_Means( mesh, means );
    for ( int i = 0; i < labels.rows; i++ )
        for ( int j = 0; j < labels.cols; j++ )
        {
            error += pow( image.at<uchar>( i, j ) - means[ labels( i, j ) ] , 2 );
        }
    error /= image.rows * image.cols;
    error = sqrt( error ) * 100.0 / 255;
    return error;
}

bool Find_Benchmarks (const Mat& image, PixelSegmentation& s, Mesh& mesh, String path, Point sizes, Benchmark& benchmark_segm, Benchmark& benchmark_mesh)
{
    int image_index;
    bool discrete = ( s.segments.size() > 0 ); // non-empty segmentation into pixel-based superpixels
    iSegmentation segmentation, mesh_superpixels, ground_truth;
    if ( discrete)
    {
        Pixel_Based_Segmentation( s, sizes, segmentation );
        benchmark_segm.RMS = Reconstruction_Error( image, segmentation.indices, s.segments.size() );
    }
    if ( ! Mesh_to_Segmentation( mesh, sizes, mesh_superpixels ) ) return false;
    // labels of faces were updated in the previous line
    Find_Means( image, mesh_superpixels.indices, mesh );
    //Find_Mesh_Colors( image, mesh_superpixels.indices, mesh );
    benchmark_mesh.RMS = Reconstruction_Error( image, mesh_superpixels.indices, mesh.n_faces() );
    image_index = 0;
    while ( Read_BSD_Human( path + "_" + Str( ++image_index ) + ".dat", ground_truth ) )
    {
        if ( discrete)
        {
            Boundary_Recall( ground_truth, segmentation, benchmark_segm );
            Undersegmentation_Errors( ground_truth, segmentation, benchmark_segm );
        }
        Boundary_Recall( ground_truth, mesh_superpixels, benchmark_mesh );
        Undersegmentation_Errors( ground_truth, mesh_superpixels, benchmark_mesh );
    }//
    benchmark_segm.Convert_in_Percents();
    benchmark_mesh.Convert_in_Percents();
    return true;
}
