#define OM_STATIC_BUILD

#include "utilities.h"
#include "colors.h"
#include "mesh.h"
#include "seeds.h"
#include "slic.h"

bool Boundary_Point (Mat_<int>const& labels, Point p)
{
    if ( p.x > 0 || p.x < labels.cols || p.y > 0 || p.y < labels.rows ) return false;
    else return true;
}

bool Boundary_Pixel (Mat_<int>const& labels, Point p, Point& neighbor)
{
    neighbor = Point( p.x-1, p.y );
    if ( p.x > 0 and labels( neighbor ) != labels( p ) ) return true;
    neighbor = Point( p.x, p.y-1 );
    if ( p.y > 0 and labels( neighbor ) != labels( p ) ) return true;
    neighbor = Point( p.x+1, p.y );
    if ( p.x+1 < labels.cols and labels( neighbor  ) != labels( p ) ) return true;
    neighbor = Point( p.x, p.y+1 );
    if ( p.y+1 < labels.rows and labels( neighbor ) != labels( p ) ) return true;
    return false;
}

bool External_Point (Mat_<int>const& labels, Point p)
{
    if ( p.x < 0 || p.x > labels.cols || p.y < 0 || p.y > labels.rows ) return true;
    else return false;
}

int Label (const Mat_<int>& m, Point2d p)
{
    if ( p.x < 0 ) { return -2; } // std::cout<<" p="<<p;
    if ( p.y > m.rows ) { return -3; }
    if ( p.x > m.cols ) {return -4; }
    if ( p.y < 0 ) { return -5; }
    return m.at<int>( Point( int( p.x ), int( p.y ) ) );
}

void Mark_Boundary (Mat_<bool>& b, Point2d p)
{
    if ( p.x > 0 and p.x < b.cols and p.y > 0 and p.y < b.rows )
        b.at<bool>( Point( int( p.x ), int( p.y ) ) ) = true;
}

int Label_on_Left (const Mat_<int>& labels, Point2d current, Point2d next)
{
    Point2d forward = next - current;
    Point2d left( forward.y, -forward.x );
    Point2d center = current + 0.5 * forward + 0.5 * left;
    return Label( labels, center );
}

int Label_on_Right (const Mat_<int>& labels, Point2d current, Point2d next)
{
    Point2d forward = next - current;
    Point2d left( forward.y, -forward.x );
    Point2d center = current + 0.5 * forward - 0.5 * left;
    return Label( labels, center );
}

Point Next_Corner (Mat_<int>const& labels, Point past, Point current, int& label_new, int& label_neighb, Mat_<bool>& passed_boundary, bool print_info)
{
    // Analyse the pixel on the left
    Point forward = current - past;
    Point2d left( forward.y, -forward.x );
    Point2d pixel_left = Point2d( past ) + 0.5 * Point2d( forward ) + 0.5 * left;
    int label_our = Label( labels, pixel_left ); //Label_on_Left( labels, past, current );
    if ( print_info ) std::cout<<" pl="<<pixel_left<<" lo="<<label_our;
    Mark_Boundary( passed_boundary, pixel_left );
   
    // Analyse the pixel on the right
    Point2d pixel_right = pixel_left - left;
    label_neighb = Label( labels, pixel_right ); // Label_on_Right( labels, past, current );
    if ( print_info ) std::cout<<" pr="<<pixel_right<<" ln="<<label_neighb;
    if ( label_our == label_neighb )
    { std::cout<<"\nError in Next_Corner:"<<" label at "<<past<<current<<" is "<<label_our<<"=neighb="<<label_neighb; exit(0); }
    
    // Boundary cases
    label_new = -1; // default value
    int x_max = labels.cols, y_max = labels.rows;
    if ( current.x == 0 && past.x > 0 ) return Point( 0, current.y + 1 ); // turn down at the left boundary
    if ( current.y == y_max && past.y < y_max ) return Point( current.x + 1, y_max ); //right at bottom
    if ( current.x == x_max && past.x < x_max ) return Point( x_max, current.y - 1 ); // up at right
    if ( current.y == 0 && past.y > 0 ) { return Point( current.x - 1, 0 ); }// turn left at the top boundary
    // Internal cases
    Point2d dif = current - past;
    dif *= 1.0 / norm( dif );
    forward = Point( int( dif.x ), int( dif.y ) ); // normalised vector of direction
    Point right( -forward.y, +forward.x ); //std::cout<<" forward="<<forward<<" right="<<right;
    int label_right = label_new; // default
    if ( ! External_Point( labels, current + right ) )
    {
        label_right = Label_on_Right( labels, current, current + forward );
        if ( label_right == label_our ) return current + right; // turn to the right
    }
    //std::cout<<" not to the right";
    int label_left = Label_on_Left( labels, current, current + forward ); //std::cout<<" ll="<<label_left;
    if ( label_left == label_our )  // continue forward
    {
        if ( label_neighb != label_right ) label_new = label_right; // a 3rd different label
        return current + forward;
    }
    if ( label_neighb != label_left ) label_new = label_left; // a 3rd different label
    return current - right; // turn to the left
}

class Arrow
{
public:
    Point tail = origin;
    Point head = origin;
    bool passed = false;
    bool empty = true;
    Arrow (Point t, Point h, bool p, bool e) { tail = t, head = h, passed = p; empty = e; }
};

bool Find_Arrow (const std::vector<Arrow>& arrows, int& i)
{
    for ( i = 0; i < arrows.size(); i++ )
        if ( !arrows[i].passed && !arrows[i].empty ) break;
    if ( i < arrows.size() ) return true; else return false;
}

bool Find_Face (const Mat_<int>& labels, Mat_<bool>& passed_boundary, Point past, Point current, std::map<Point,MVH,Compare_Points>& all_vhandles, std::vector<Arrow>& arrows, Mesh& mesh, MFH& fh, bool print_info)
{
    MVH vh;
    int label_new, label_neighb;
    Point initial = current, next;
    std::vector<MVH> face_vhandles; // for a current face
    face_vhandles.reserve( labels.rows * labels.cols / 100 );
    while ( true )
    {
        if ( print_info ) std::cout<<"\npast="<<past<<" cur="<<current;
        next = Next_Corner( labels, past, current, label_new, label_neighb, passed_boundary, print_info );
        if ( External_Point( labels, next ) ) std::cout<<"\nError: external point "<<next;
        if ( print_info) std::cout<<" next="<<next;
        auto it = all_vhandles.find( current );
        if ( it !=  all_vhandles.end() or // past vertex
            ! Collinear( past, current, next ) or // new corner at the current point
            ( label_new >=0 ) ) // future vertex at the current point
        {
            //if ( print_info) std::cout<<"\nc="<<current<<" l_our="<<label_our<<" l_new="<<label_new;
            if ( it !=  all_vhandles.end() ) vh = it->second;
            else
            {
                vh = mesh.add_vertex( MP( current.x, current.y, 0 ) );
                all_vhandles.insert( make_pair( current, vh ) );
            }
            face_vhandles.push_back( vh );
            // Check if the new label to the right of the next arrow represents a new superpixel
            if ( label_new >= 0 and !arrows[ label_new ].passed and arrows[ label_new ].empty ) // new superpixel
            {
                arrows[ label_new ].empty = false;
                if ( print_info) std::cout<<"\na"<<current<<next<<"l_new="<<label_new<<" l_neighb="<<label_neighb;
                arrows[ label_new ].tail = next;
                arrows[ label_new ].head = current;
            }
        }
        if ( next == initial ) break; // closed boundary
        past = current;
        current = next;
    }
    fh = mesh.add_face( face_vhandles );
    mesh.data( fh ).area = -Area_Face( mesh, fh ); // negative orientation
    return true;
}

bool Find_Faces (Mat_<int>& labels, Mat_<bool>& passed_boundary, std::map<Point,MVH,Compare_Points>& all_vhandles, std::vector<PixelSegment>& segments, std::vector<Arrow>& arrows, Mesh& mesh, bool print_info)
{
    MVH vh;
    MFH fh, fh_extra;
    double area_found;
    int check_face = -1;
    Point initial, past, current, next, p, neighbor;
    int label, label_new;
    int extra_faces = 0;
    //std::cout<<"\nsup="<<(int)segments.size()<<" labels.size="<<labels.size()<<" rows="<<labels.rows<<" cols="<<labels.cols;
    while ( Find_Arrow( arrows, label ) )
    {
        past = arrows[ label ].tail;
        current = arrows[ label ].head;
        if ( label == check_face ) print_info = true; else print_info = false;
        if ( print_info ) std::cout<<"\ns"<<label<<" start: "<<past<<"->"<<current;
        arrows[ label ].passed = true; // this flag means that the superpixel was considered
        
        // Start the loop along the boundary from the arrow (past,current)
        Find_Face( labels, passed_boundary, past, current, all_vhandles, arrows, mesh, fh, print_info );
        mesh.data( fh ).label = label;
        //std::cout<<" l"<<label<<"#f"<<mesh.n_faces();
        
        // Check completeness
        area_found = mesh.data( fh ).area;
        while ( area_found < segments[ label ].area() )
        {
            if ( print_info ) std::cout<<"\nFace "<<label<<" incomplete: "<<area_found<<" orig="<<segments[ label ].area();
            for ( int i = 0; i < segments[ label ].xs.size(); i++ )
            {
                p = Point( segments[ label ].xs[i], segments[ label ].ys[i] );
                if ( passed_boundary( p ) or ! Boundary_Pixel( labels, p, neighbor ) ) continue;
                Point dif = p - neighbor, tail = p, head = p;
                if ( dif.y < 0 ) head.x += 1;
                else if ( dif.y > 0 ) tail.x += 1;
                else if ( dif.x < 0 ) tail.y += 1;
                else if ( dif.x > 0 ) head.y += 1;
                if ( print_info ) std::cout<<" p="<<p<<" nb="<<neighbor<<" "<<tail<<"->"<<head;
                label_new = (int)segments.size() + extra_faces; // label of new face
                extra_faces++; // new face added
                past = tail; //Point( tail.y, tail.x );
                current = head; //Point( head.y, head.x );
                //arrows.push_back( Arrow( past, current, true, true ) );
                Find_Face( labels, passed_boundary, past, current, all_vhandles, arrows, mesh, fh_extra, print_info );
                mesh.data( fh_extra ).label = label_new;
                Update_Labels( mesh, fh_extra, labels );
                area_found += mesh.data( fh_extra ).area;
                if ( print_info )
                    std::cout<<"\nl"<<mesh.data( fh_extra ).label<<" extra_area="<<mesh.data( fh_extra ).area<<" total="<<area_found;
                break;
            }
        }//
        if ( label == check_face ) Print_Face( mesh, fh );
    }
    /*
    if ( segments.size() + extra_faces != mesh.n_faces() )
    {
        std::cout<<"\nError in Find_Faces:";
        std::cout<<" total_area="<<total_area;
        int l = 83; std::cout<<"\ns83="<<segments[l].xs.size();
        for ( int i = 0; i < segments[l].xs.size(); i++ ) std::cout<<Point( segments[ l ].xs[i], segments[ l ].ys[i] );
    }*/
    return true;
}

bool Labels_to_Mesh (const Mat& image, PixelSegmentation& segmentation, Mesh& mesh)
{
    Mat_<int> labels = segmentation.segmentation_data - 1;
    //Transpose( segmentation.segmentation_data - 1, labels );
    int num_superpixels = (int)segmentation.number_segments();
    Mat_<bool> passed_boundary( labels.size(), false );
    //for ( int j = 0; j < labels.cols; j++ ) std::cout<<" l"<<j<<"="<<labels( 0, j );
    mesh.clear();
    Arrow a( Point(1,0), Point(0,0), false, true );
    std::vector<Arrow> arrows( num_superpixels, a ); // default arrows
    std::map<Point,MVH,Compare_Points> all_vhandles; // vertex_handles associated with pixel positions
    //Print_Matrix
    //for ( int i = 319; i < 325; i++ ) { std::cout<<"\n";
    //    for ( int j = 0; j < 9; j++ ) std::cout<<" l"<<Point( i, j )<<"="<<labels.at<int>( Point( i, j ) ); }
    std::vector<double> areas;
    areas.reserve( segmentation.segments.size() );
    for ( auto s : segmentation.segments ) { areas.push_back( s.area() ); } //std::cout<<" a"<<areas.size()<<"="<<s.area(); }
    
    arrows[ labels.at<int>(0,0) ].empty = false; // initialization of the first arrow
    //std::cout<<" arrows="<<arrows.size();
    if ( ! Find_Faces( labels, passed_boundary, all_vhandles, segmentation.segments, arrows, mesh, false ) ) return false;
    // Check the total area
    double area = 0;
    for ( auto f_it = mesh.faces_begin(); f_it != mesh.faces_end(); f_it++ ) area += mesh.data( *f_it ).area;
    if ( area != image.rows * image.cols ) std::cout<<"\nError in Labels_to_Mesh: area="<<area;
    Find_Means( image, labels, mesh );
    return true;
}

void Transpose (const Mat& source, Mat& result)
{
    result = Mat( source.cols, source.rows, CV_8U );
    for ( int i = 0; i < result.rows; i++ )
        for ( int j = 0; j < result.cols; j++ )
            result.at<uchar>(i,j) = source.at<uchar>(j,i);
}

int main ()
{
    String data_path = "/Users/philsmith/Documents/Xcode Projects/Meshes/";
    String input_folder = data_path + "BSDS500/data/images/test/";
    String output_folder = data_path + "Output/";
    
    std::vector<std::string> image_names;
    std::vector<Point> image_sizes;
    
    Read_Image_Names_Sizes( output_folder + "BSD500sizes.txt", image_names, image_sizes );
    
    int num_images = (int)image_names.size();
    
    int num_superpixels = 100; // 60, 150, 260, 384,
    String method = "SLIC"; //"SEEDS"; // "Voronoi"
    
    String parameter = Str( num_superpixels );
    
    String type_method = method + parameter;
    bool user_image = true;
    String user_image_name = "36046";
    bool cameraman = false;
    
    Create_Directory( output_folder );
    
    bool find_benchmarks = false; //false; //
    bool run_segmentation = true;
    bool save_images = true;
    bool input_gray = false;
    bool means_gray = true;
    
    int scale_up = 4;
    int scale_down = 1;
    int border = 1;
    
    double color_offset = 30; // Was 30
    double min_color_dif = 2;
    double min_area_face = 10;
    int eps_radius = 12;
    
    String type_images = "BSD"; //"user";
    
    //if ( num_images > 1 ) save_images = false;
    
    std::ofstream out_segm, out_mesh;
    
    if ( find_benchmarks )
    {
        out_segm.open( output_folder + type_method + "_segmentation.txt" );
        out_mesh.open( output_folder + type_method + "_mesh.txt" );
    }
    
    // Loop over images
    
    Mesh mesh;
    std::vector<int> failed_images; //( 473, -1 );
    String name, path;
    Benchmark average_segm( 0 ), average_mesh( 0 );
    Point sizes;
    Mat image, image_means, image_mesh, image_segm, image_input, image_scaled, image_gray, image_color, image_reconstruction;
    
    for (int i = 0; i < 1; i++ )
    {
        if ( type_images == "BSD" ) name = image_names[i];
        
        if (user_image) name = user_image_name;
        
        if (cameraman) name = "Cameraman256";
        
        Print( Str(i) + " " + name + "\n", out_segm, out_mesh );
        
        output_folder = data_path + "Output/" + name + "/";
        
        Create_Directory( output_folder );
        
        output_folder += type_method + "/";
        
        Create_Directory( output_folder );
        
        path = input_folder + name + ".jpg";
        
        if (cameraman) path = "/Users/philsmith/Documents/Xcode Projects/Superpixel Meshes/Cameraman/Cameraman256.jpg";
        
        if ( input_gray ) image_input = cv::imread( path, IMREAD_GRAYSCALE );
        else image_input = cv::imread( path, IMREAD_COLOR );
        
        if ( !image_input.data ) { std::cout<<"\nNo file "<<path; exit(0); }
        
        if ( scale_down == 1 ) image_scaled = image_input.clone();
        else resize( image_input, image_scaled, Size( image_input.size() / scale_down ) );
        
        if ( input_gray ) { image_gray = image_scaled.clone(); cvtColor( image_scaled, image_color, COLOR_GRAY2BGR ); }
        else { cvtColor( image_scaled, image_gray, COLOR_BGR2GRAY ); image_color = image_scaled.clone(); }
        
        image = image_gray.clone();
        
        if ( method == "SLIC" ) image = image_color.clone();
        
        //Write_Image( output_folder + "Image.png", image );
        
        // image is the input for the method below, image_input has the original size and color
        
        Benchmark benchmark_segm( 1 ), benchmark_mesh( 1 );
        PixelSegmentation segmentation;
        
        if ( run_segmentation )// Run pixel segmentation
        {
            auto start = chrono::high_resolution_clock::now();
            
            if ( method == "SLIC") segmentation = run_slic( image, num_superpixels, 1 );
            
            else segmentation = run_seeds( image, num_superpixels );
            
            benchmark_segm.time = chrono::duration_cast<chrono::milliseconds>( chrono::high_resolution_clock::now() - start ).count();
            
            cout<<"Number of segments = " <<segmentation.number_segments() << endl;
            
            if ( method == "SLIC" ) image = image_gray.clone();
            
            if (save_images)
            {
                segmentation.compute_mean( image, image_means );
                
                if ( scale_down > 1 ) Scale_Image( image_means, image_means, scale_down );
                
                //Write_Image( output_folder + "Means.png", image_means );
            }
            
            // Build a boundary mesh from pixel-based superpixels
            
            if (!Labels_to_Mesh( image, segmentation, mesh )) { failed_images.push_back( i ); continue; }
            
            if ( scale_down > 1 )
            {
                Scale_Mesh( mesh, scale_down ); // points are multiplied by scale_down to get the input size
                
                Enlarge_Mesh( mesh, image_input.size() ); //
            }
            
            if ( save_images )
            {
                image_segm = image_input.clone(); // big input
                
                if ( input_gray ) cvtColor( image_segm, image_segm, COLOR_GRAY2BGR ); // for drawing colored boundaries
                copyMakeBorder( image_segm, image_segm, border, border, border, border, BORDER_CONSTANT, White );
                image_mesh = image_segm.clone();
                Write_Mesh( mesh, output_folder + type_method + ".off" );
                Draw_Mesh( mesh, scale_up, border, Red, Red, 2, image_segm );
                Write_Image( output_folder + type_method + ".png", image_segm );
            }
            
            sizes = image_sizes[i];
            
            if ( find_benchmarks )
            {
                benchmark_segm.num_faces = (int)mesh.n_faces();
                benchmark_segm.num_edges = (int)mesh.n_edges();
                benchmark_segm.num_vertices = (int)mesh.n_vertices();
                //benchmark_segm.RMS = segmentation.reconstruction_error( image );
                Find_Compactness( mesh, sizes.x * sizes.y, benchmark_segm.comp );
            }
            
            // Build a resolution-independent mesh from pixel-based superpixels
            
            start = chrono::high_resolution_clock::now();
            Remove_Small_Faces( mesh, min_area_face );
            Straighten_Chains( mesh, color_offset * scale_down, min_color_dif );
            Merge_Faces( mesh, min_color_dif );
            Write_Mesh( mesh, output_folder + type_method + "_RIMe.off" );
            benchmark_mesh.time = chrono::duration_cast<chrono::milliseconds>( chrono::high_resolution_clock::now() - start ).count();
        }
        
        else Read_Mesh( mesh, output_folder + name + ".off" );
        
        if ( save_images )
        {
            Draw_Mesh( mesh, scale_up, border, Red, Red, 2, image_mesh ); // labels of faces are the same as in the segmentation
            Write_Image( output_folder + type_method + "_RIMe.png", image_mesh );
        }
        
        // Find benchmarks for a resolution-independent mesh
        
        if ( find_benchmarks )
        {
            benchmark_mesh.num_faces = (int)mesh.n_faces();
            benchmark_mesh.num_edges = (int)mesh.n_edges();
            benchmark_mesh.num_vertices = (int)mesh.n_vertices();
            Find_Compactness( mesh, sizes.x * sizes.y, benchmark_mesh.comp );
        }
        
        // Find BR2, CUE, ASA, discrete RMS for segm and mesh from the same ground truth
        
        if ( find_benchmarks and type_images == "BSD" )
        {   // labels of faces are updated in the next line for computing benchmarks and means
            if (!Find_Benchmarks( image, segmentation, mesh, data_path + "BSD500/data/" + name, sizes, benchmark_segm, benchmark_mesh ))
            { failed_images.push_back( i ); continue; }
            Print( " t=" + Round1( benchmark_segm.time ) + " f=" + Str( benchmark_segm.num_faces ) + " e=" + Str( benchmark_segm.num_edges ), out_segm );
            benchmark_segm.Print( std::cout, out_segm );
            average_segm.Add( benchmark_segm );
            int act_images = i + 1 - (int)failed_images.size();
            //int act_images = num_images - i - (int)failed_images.size();
            std::cout<<"\nmesh:";
            Print( " t=" + Round1( benchmark_mesh.time ) + " f=" + Str( benchmark_mesh.num_faces ) + " e=" + Str( benchmark_mesh.num_edges ), out_mesh );
            benchmark_mesh.Print( std::cout, out_mesh );
            average_mesh.Add( benchmark_mesh );
            average_segm.Print_All( act_images, cout, out_segm );
            average_mesh.Print_All( act_images, cout, out_mesh );
        }
        
        if ( save_images )
        {
            image_reconstruction = image_gray.clone();
            //Draw_Reconstruction( mesh, scale_up, image_reconstruction );
            //Write_Image( output_folder + "Reconstruction.png", image_reconstruction  );
        }
    }
    
    if ( find_benchmarks )
    {
        int act_images = num_images - (int)failed_images.size();
        std::cout<<"\n" + type_method + " averages for "<<act_images<<" images, color_offset="<<color_offset<<", min_color_dif="<<min_color_dif;
        average_segm.Divide_by_Images( act_images );
        Print( "\ngrid t=" + Round1( average_segm.time )+ " f=" + Round1( average_segm.num_faces ) + " e=" + Round1( average_segm.num_edges ), out_segm );
        average_segm.Print( std::cout, out_segm ); out_segm.close();
        average_mesh.Divide_by_Images( act_images );
        Print( "\nRIMe t=" + Round1( average_mesh.time ) + " f=" + Round1( average_mesh.num_faces ) + " e=" + Round1( average_mesh.num_edges ), out_mesh );
        average_mesh.Print( std::cout, out_mesh ); out_mesh.close();
    }
    
    if ( failed_images.size() > 0 ) { std::cout<<"\nFailed images:"; for ( auto i : failed_images ) std::cout<<" "<<i; }
    
    cout << endl;
    
    return 0;
}
