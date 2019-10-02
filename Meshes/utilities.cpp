//
//  utulities.cpp
//  RIS
//
//  Created by Vitaliy Kurlin on 21/10/2016.
//  Copyright © 2016 Vitaliy Kurlin. All rights reserved.
//
#include "utilities.h"

String Str (int n) { return std::to_string( n ); }

string Round (double v) { char b[999]; sprintf( b, "%.0f", v ); return string(b); }
string Round1 (double v) { char b[999]; sprintf( b, "%.1f", v ); return string(b); }
string Round2 (double v) { char b[999]; sprintf( b, "%.2f", v ); return string(b); }
string Round3 (double v) { char b[999]; sprintf( b, "%.3f", v ); return string(b); }
string Round4 (double v) { char b[999]; sprintf( b, "%.4f", v ); return string(b); }
string Round5 (double v) { char b[999]; sprintf( b, "%.5f", v ); return string(b); }

bool Decreasing_Values (Index_Value const& p1, Index_Value const& p2) { return p1.value > p2.value; }
bool Increasing_Values (Index_Value const& p1, Index_Value const& p2) { return p1.value < p2.value; }

void Write_Image (String const& name, Mat const& image)
{
    try{ imwrite( name, image ); }
    catch (std::runtime_error& ex) { fprintf(stderr, "Exception converting to PNG: %s\n", ex.what()); }
}

void Print (String s) { std::cout<<s; }
void Print (String s, std::ostream& fout) { fout<<s; std::cout<<s; }
void Print (String s, std::ostream& file1, std::ostream& file2 ) { std::cout<<s; file1<<s; file2<<s; }

void Create_Directory (String name)
{
    boost::filesystem::path dir( name );
    boost::filesystem::create_directory( dir );
}

double Det (Point2d a, Point2d b) { return a.x * b.y - a.y * b.x; }

bool Find_Sine_Cosine (Point2d a, Point2d b, Point2d c, double& sine, double& cosine)
{
    Point2d v1 = a-b, v2 = c-b;
    double l1 = norm( v1 ), l2 = norm( v2 ), prod = l1 * l2;
    if ( prod < 1e-8  ) return false ;
    double det = Det( v1, v2 );
    sine = det / prod;
    cosine = ( v1.x * v2.x + v1.y * v2.y ) / prod;
    return true;
}

bool Collinear (Point2d a, Point2d b, Point2d c)
{
    if ( Det( a-b, b-c ) != 0 ) return false;
    else return true;
}

void Points_to_Line_Coef (Point2d p1, Point2d p2, double& a, double& b, double& c)
{ // (p1.y - p2.y) x + (p2.x - p1.x) y + (p1.x p2.y - p1.y p2.x)
    a = p1.y - p2.y;
    b = p2.x - p1.x;
    c = Det( p1, p2 );
}

bool Distance_Point_to_Line (Point2d p, double a, double b, double c, double& d)
{
    double l = a*a + b*b;
    if (l <= 1e-32) return false;
    d = fabs(a * p.x + b * p.y + c) / std::sqrt(l);
    return true;
}

bool Read_Image_Names_Sizes (std::string file_name, std::vector<std::string>& image_names, std::vector<Point>& image_sizes)
{
    std::ifstream file;
    file.open( file_name );
    image_names.clear();
    image_sizes.clear();
    if ( ! file.is_open() ) { std::cout << "\nError opening " << file_name; return false; }
    int x, y;
    std::string line, name;
    while( getline( file, line ) )  // read one line from ifs
    {
        std::istringstream iss( line ); // access line as a stream
        iss >> name >> x >> y;
        image_names.push_back( name );
        image_sizes.push_back( Point2i( x, y ) );
    }
    return true;
}

void Scale_Image (const Mat& img, Mat& img_scaled, int scale)
{
    resize( img, img_scaled, img.size() * scale );
    for ( int i = 0; i < img.rows; i++ )
        for ( int j = 0; j < img.cols; j++ )
            for ( int ii = 0; ii < scale; ii++ )
                for ( int jj = 0; jj < scale; jj++ )
                    img_scaled.at<char>( i * scale + ii, j * scale + jj ) = img.at<char>( i, j );
}

bool Scale_Image (const Mat& image_old, int border, int scale_down, Mat& image_new)
{
    //image_new = Rect( image_old, );
    
    return true;
}

void Transpose (const Mat_<int>& matrix, Mat_<int>& transpose)
{
    transpose = Mat_<int>( matrix.size().width, matrix.size().height );
    for ( int i = 0; i < transpose.rows; i++ )
        for ( int j = 0; j < transpose.cols; j++ )
            transpose( i, j ) = matrix( j, i );
}

bool Is_Nonstrictly_Inside_Polygon (Polygon const& poly, P2 point)
{
    switch( CGAL::bounded_side_2( poly.vertices_begin(), poly.vertices_end(), point, K() ) )
    {
        case CGAL::ON_BOUNDED_SIDE : return true;
        case CGAL::ON_BOUNDARY: return true;
        case CGAL::ON_UNBOUNDED_SIDE: return false;
    }
}

bool Read_BSD_Human (std::string path, iSegmentation& segmentation)
{
    std::ifstream file;
    //std::cout<<"\nFile to BSD human: "<<path;
    if ( ! boost::filesystem::exists( path ) ) { return false; } //cout<<"\nFile "<<path<<" not found";
    file.open( path );
    int32_t* dims = new int32_t[2];
    file.read( (char*)dims, sizeof(int32_t)*2 );
    int32_t width = dims[0];
    int32_t height = dims[1];
    //std::cout<<"\nFile "<<path<<" found"<<" w="<<width<<" h="<<height;
    segmentation.boundary.clear();
    int num_pixels = width * height;
    int32_t* data = new int32_t[ num_pixels ];
    uchar* bdata = new uchar[ num_pixels ];
    file.read( (char*)data, num_pixels * sizeof(int32_t) ) ;
    file.read( (char*)bdata, num_pixels * sizeof(uchar) );
    segmentation.indices = Mat_<int32_t>( height, width, data ).clone();
    segmentation.num_faces = 0;
    for ( int i = 0; i < height; i++ )
        for ( int j = 0; j < width; j++ )
        {
            if ( segmentation.num_faces < segmentation.indices( i, j ) ) segmentation.num_faces = segmentation.indices( i, j );
            segmentation.indices( i, j )--;
        }
    //std::cout<<"\nnum_faces="<<segmentation.num_faces; //<<" min="<<*min_element( data, data + num_pixels );;
    Mat_<uchar> boundary_mesh = Mat_<uchar>( height, width, bdata ).clone();
    for ( int i = 0; i < height; i++ )
        for ( int j = 0; j < width; j++ )
            if ( boundary_mesh( i, j ) == 1 ) segmentation.boundary.push_back( Point( j, i ) );
    delete data;
    delete bdata;
    delete dims;
    file.close();
    return true;
}

bool Read_BSD_Human (std::string path, PixelSegmentation& segmentation)
{
    std::ifstream file;
    //std::cout<<"\nFile to BSD human: "<<path;
    if ( ! boost::filesystem::exists( path ) ) { cout<<"\nFile "<<path<<" not found"; return false; } //
    PixelSegmentation::load_from_file( file, segmentation );
    return true;
}

bool Pixel_Near_Boundary (Point p, Mat_<bool> boundary_mask, int radius)
{
    for ( int j = max( 0, p.x - radius ); j < min( boundary_mask.cols - 1, p.x + radius ); j++ )
        for ( int i = max( 0, p.y - radius ); i < min( boundary_mask.rows - 1, p.y + radius ); i++ )
            if ( boundary_mask( i, j ) ) return true;
    return false;
}

void Boundary_Recall (iSegmentation const& s0, iSegmentation const& s1, Benchmark& b) // s0 ground, s1 experiment
{
    int bpixels2 = 0;
    for ( size_t k = 0; k < s0.boundary.size(); k++ )
        if ( Pixel_Near_Boundary ( s0.boundary[k], s1.boundary_mask, 2) ) bpixels2++;
    double br2 = 1.0 * bpixels2 / s0.boundary.size(); //std::cout<<" BR2="<< br2;
    if ( b.BR2 < br2 ) b.BR2 = br2;
    int bpixels1 = 0;
    for ( size_t k = 0; k < s0.boundary.size(); k++ )
        if ( Pixel_Near_Boundary ( s0.boundary[k], s1.boundary_mask, 1) ) bpixels1++;
    double br1 = 1.0 * bpixels1 / s0.boundary.size(); //cout<<" BR2="<< br2;
    if ( b.BR1 < br1 ) b.BR1 = br1;
}

void Undersegmentation_Errors (iSegmentation const& s0, iSegmentation const& s1, Benchmark& benchmark) // s0 ground, s1 exp
{
    Benchmark b(0);
    double max_area_in = 0;
    int cur_area_in; //, current_use;
    std::vector<int> num_pixels( s0.num_faces, 0 ); // pixels contained in ground truth segments
    for ( int k = 0; k < s1.superpixels.size(); k++ ) // for each fixed superpixel
    {
        max_area_in = 0;
        num_pixels.assign( s0.num_faces, 0 );
        for ( size_t l = 0; l < s1.superpixels[k].size(); l++ )
        {
            num_pixels[ s0.indices( s1.superpixels[k][l] ) ]++;
            cur_area_in = num_pixels[ s0.indices( s1.superpixels[k][l] ) ];
            if ( max_area_in < cur_area_in ) max_area_in = cur_area_in;
        }
        b.ASA += max_area_in;
        if ( max_area_in == (int)s1.superpixels[k].size() ) continue; // the whole superpixel is inside a ground truth segment
        //cout<<"\ns"<<k;
        b.CUE += (int)s1.superpixels[k].size() - max_area_in;
        /*
         for ( int l = 0; l < num_pixels.size(); l++ )
         if ( num_pixels[l] > 0 ) // non-empty intersection with the l-th segment
         b.USE += min( num_pixels[l], (int)s1.superpixels[k].size() - num_pixels[l] );
         */
    }//
    b.Divide_by_Pixels( s1.indices.cols * s1.indices.rows );
    if ( benchmark.CUE > b.CUE ) benchmark.CUE = b.CUE;
    //if ( benchmark.USE > b.USE ) benchmark.USE = b.USE;
    if ( benchmark.ASA > b.ASA ) benchmark.ASA = b.ASA;
    //cout<<" final USE="<<benchmark.USE;
}
