//#pragma once

//#include <opencv2/core/core.hpp>
//#include <cstdint>
#include "utilities.h"

using namespace std;

#ifndef SEGMENTATION_H 
#define SEGMENTATION_H

class Segmentation {
    public:
        double boundary_recall(Segmentation& ground_truth);
        double boundary_recall(Segmentation& ground_truth, int epsilon);
        double undersegmentation_error(Segmentation& ground_truth);
        double achievable_segmentation_accuracy(Segmentation& ground_truth);
        double reconstruction_error(const cv::Mat& image);
        /// Determine the area of the intersection of the face with <label> with the pixel at (x,y).
        double intersection_area(int label, int x, int y);
        // /// Get the number of regions in this segmentation.
        unsigned long number_segments();
        // /// The pixel width of the represented image.
        // int size();

        int width, height;
};

struct intersection_result {
    int area_in;
    int area_out;
};

class PixelSegmentation;

class PixelSegment {
    public:
        PixelSegmentation& segmentation;
        int32_t label;
        int x1;
        int x2;
        int y1;
        int y2;
        vector<int> xs;
        vector<int> ys;
        double perimeter();
        double area();
        bool bbox_intersect(PixelSegment& other);
        intersection_result intersection(PixelSegment& other);
        PixelSegment(PixelSegmentation& seg);
};

class PixelSegmentation : public Segmentation {
    public:
        static PixelSegmentation& load_from_file(ifstream& file);
        static void load_from_file (std::ifstream& file, PixelSegmentation& s);
        void output_to_file(ofstream file);
        double reconstruction_error(const cv::Mat& image);
        double compactness();
        double boundary_recall(PixelSegmentation& ground_truth);
        double boundary_recall(PixelSegmentation& ground_truth, int epsilon);
        double undersegmentation_error(PixelSegmentation& ground_truth);
        double achievable_segmentation_accuracy(PixelSegmentation& ground_truth);
        cv::Mat get_boundary_pixels() const;
        void initialise_segments();
        void compute_mean(cv::Mat& image, cv::Mat& output);

        PixelSegmentation(int32_t* data, uchar* b_data, int width, int height);
        PixelSegmentation(cv::Mat& data, cv::Mat& boundary_data);
        PixelSegmentation();
    
        unsigned long number_segments();
    
        cv::Mat segmentation_data;
        cv::Mat boundary_data;
        vector<PixelSegment> segments;

};

class MeshSegmentation : Segmentation {
    public:
        void load_from_file(ifstream file);
        void output_to_file(ofstream file);

        MeshSegmentation();

};
#endif