#include <iostream>
using namespace std;
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
using namespace cv;

#include "gvfsnake.h"

#if 0
int cursorx = 0;
int cursory = 0;

void on_mouse
    (
    int event,
    int x,
    int y,
    int flags,
    void* userdata
    )
{
    if (event == CV_EVENT_MOUSEMOVE && (flags & CV_EVENT_FLAG_LBUTTON))
        {
        cout << "x: " << x << ", y: " << y << endl;
        }

}
#endif

double contour_pts_x[100] = { 250.3502, 246.33792, 241.7719, 236.60872, 231.03105, 225.13277, 218.97568, 212.54716, 205.96244, 199.30511, 192.48044,
185.58844, 178.60824, 171.53826, 164.43261, 157.23314, 150.0169, 142.76033, 135.50221, 128.2431, 121.06032, 113.96166, 106.97162, 100.30941, 94.209824,
88.763155, 83.995415, 79.854739, 76.176088, 72.976062, 70.38184, 68.073567, 66.086524, 64.491298, 63.059189, 61.870738, 60.880399, 59.995726, 59.290111,
58.295225, 57.124785, 56.099996, 55.240014, 54.552689, 54.084658, 53.920722, 54.213592, 55.003703, 56.402774, 58.509822, 61.382511, 65.057925, 69.374158,
74.182682, 79.503169, 85.216432, 91.219317, 97.468918, 103.99454, 110.7219, 117.57963, 124.57452, 131.69826, 138.87638, 146.11103, 153.37194, 160.61904,
167.83743, 174.95768, 181.99897, 188.92386, 195.63975, 202.21461, 208.57259, 214.62692, 220.48295, 226.01711, 231.18304, 236.10608, 240.6117, 244.5881,
248.22511, 251.38049, 254.16078, 256.67276, 258.70408, 260.41658, 261.87783, 262.86972, 263.58069, 264.04099, 264.05068, 263.76181, 263.19287, 262.21149,
261.13227, 259.71121, 257.79795, 255.33481, 252.20847 };

double contour_pts_y[100] = { 238.2637, 244.31107, 249.96139, 255.07378, 259.71992, 263.94801, 267.79835, 271.18494, 274.25021, 277.15666, 279.65042,
281.93813, 283.956, 285.6162, 287.13427, 288.11413, 288.91764, 289.28789, 289.23744, 288.92075, 287.82622, 286.3419, 284.35841, 281.46861, 277.65401,
272.86128, 267.37826, 261.42269, 255.17077, 248.64737, 241.86318, 234.97953, 227.99049, 220.90336, 213.78181, 206.61365, 199.4165, 192.20493, 184.97327,
177.77982, 170.60898, 163.4156, 156.20123, 148.96915, 141.71955, 134.45539, 127.19586, 119.9768, 112.85444, 105.91222, 99.253068, 93.006061, 87.173127,
81.730296, 76.781829, 72.301558, 68.223001, 64.525001, 61.329517, 58.598541, 56.217025, 54.250683, 52.837055, 51.757156, 51.082856, 51.125813, 51.529528,
52.359978, 53.7891, 55.550043, 57.750228, 60.511059, 63.585732, 67.102972, 71.109727, 75.399006, 80.1073, 85.20729, 90.543198, 96.24357, 102.30977,
108.59488, 115.13999, 121.8465, 128.66169, 135.63797, 142.69423, 149.81042, 157.00841, 164.23535, 171.48562, 178.75163, 186.00796, 193.2491, 200.44853,
207.63375, 214.75254, 221.75485, 228.58717, 235.1462 };

int main
    (
    int argc,
    char* argv[]
    )
{
    Mat im = imread("testimage2.png");
    int debugidx = 0;
    for (debugidx = 0; debugidx < 100; ++debugidx)
        {
        circle(im, Point(contour_pts_x[debugidx], contour_pts_y[debugidx]), 2, Scalar(255, 0, 0), -1);
        }
    Mat gray;
    cvtColor(im, gray, CV_BGR2GRAY);
    Mat grayf;
    gray.convertTo(grayf, CV_32FC1);
    normalize(grayf, grayf, 0.0, 1.0, NORM_MINMAX, CV_32FC1);
    init_stack();
    const int32_t inwidth = grayf.cols;
    const int32_t inheight = grayf.rows;
    const float* in = (float*)grayf.data;
    float* Ex = (float*)alloc_from_stack(inwidth * inheight * sizeof(float));
    float* Ey = (float*)alloc_from_stack(inwidth * inheight * sizeof(float));
    calc_external_force(in, Ex, Ey, inwidth, inheight);
    double** internal_energy_coefficients = create_internal_force_pentadiagonal_matrix(100);
    int32_t i = 0;
    double** px = alloc_matrix(100, 1);
    double** py = alloc_matrix(100, 1);
    memcpy(px[1] + 1, contour_pts_x, 100 * sizeof(double));
    memcpy(py[1] + 1, contour_pts_y, 100 * sizeof(double));
    for (debugidx = 1; debugidx <= 100; ++debugidx)
        {
        circle(im, Point(px[debugidx][1], py[debugidx][1]), 2, Scalar(255, 0, 0), -1);
        }
    for (i = 0; i < ITER_TIMES_CONTOUR; ++i)
        {
        contour_update(px, py, Ex, Ey, internal_energy_coefficients, 100, inwidth, inheight);
        }
    for (debugidx = 1; debugidx <= 100; ++debugidx)
        {
        circle(im, Point(px[debugidx][1], py[debugidx][1]), 2, Scalar(0, 0, 255), -1);
        }
    imshow("demo", im);
    waitKey(0);
    free_stack();
    return EXIT_SUCCESS;
}