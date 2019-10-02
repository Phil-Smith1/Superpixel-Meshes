#ifndef slic_h
#define slic_h
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include "segmentation.h"
#include <fstream>
#include <iostream>

using namespace cv;
#define NR_ITERATIONS 10

PixelSegmentation run_slic(Mat& image, int target_superpixel_number, int m);

#endif /* donald_slic_h */