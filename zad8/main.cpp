#include <iostream>
#include <opencv2/highgui.hpp>
#include <opencv2/opencv.hpp>
#include <vector>

using namespace std;
using namespace cv;

bool compareContourAreas(vector<Point> contour1, vector<Point> contour2) {
    double i = fabs(contourArea(Mat(contour1)));
    double j = fabs(contourArea(Mat(contour2)));
    return (i < j);
}

int main(int argc, char **argv) {
    vector<int> lower = {0, 141, 87};
    vector<int> upper = {15, 255, 255};
    namedWindow("options", cv::WINDOW_AUTOSIZE);
    createTrackbar("lh", "options", &lower[0], 255);
    createTrackbar("ls", "options", &lower[1], 255);
    createTrackbar("lv", "options", &lower[2], 255);
    createTrackbar("hh", "options", &upper[0], 255);
    createTrackbar("hs", "options", &upper[1], 255);
    createTrackbar("hv", "options", &upper[2], 255);

    VideoCapture camera(0);
    if (!camera.isOpened())
        return -1;

    while (waitKey(1) != 27) {
        Mat src_img;

        camera.read(src_img);
        Mat bw_img;
        Mat red_img;
        Mat img_edges;

        cvtColor(src_img, bw_img, COLOR_BGR2HSV);
        inRange(bw_img, Scalar(lower[0], lower[1], lower[2]), Scalar(upper[0], upper[1], upper[2]), red_img);

        Mat k = getStructuringElement(MORPH_ELLIPSE, {10, 10});
        morphologyEx(red_img, red_img, MORPH_CLOSE, k);

        Canny(red_img, red_img, 100, 100);

        auto kern = getStructuringElement(MORPH_ELLIPSE, {50, 50});
        morphologyEx(red_img, red_img, MORPH_CLOSE, kern);

        vector<vector<Point>> contours;
        vector<vector<Point>> contours_biggest;
        findContours(red_img, contours, RETR_LIST, CHAIN_APPROX_SIMPLE);
        drawContours(src_img, contours, -1, {255, 255, 255});

        sort(contours.begin(), contours.end(), compareContourAreas);

        if (contours.size() > 1) {
            vector<Point> biggest = contours[contours.size() - 1];
            vector<Point> secondBiggest = contours[contours.size() - 2];

            auto m = moments(biggest, false);
            Point p = {(int) (m.m10 / m.m00), (int) (m.m01 / m.m00)};
            Point p_2 = {(int) ((m.m10 / m.m00) - 4), (int) ((m.m01 / m.m00) - 4)};
            int x = (int) (m.m10 / m.m00);
            int y = (int) (m.m01 / m.m00);
            ellipse(src_img, p, Size{5, 5}, 0.0, 0, 360, {0, 255, 255}, 3);

            auto m1 = moments(secondBiggest, false);
            Point p1 = {(int) (m1.m10 / m1.m00), (int) (m1.m01 / m1.m00)};
            Point p1_2 = {(int) ((m1.m10 / m1.m00) - 4), (int) ((m1.m01 / m1.m00) - 4)};
            int x1 = (int) (m1.m10 / m1.m00);
            int y1 = (int) (m1.m01 / m1.m00);
            ellipse(src_img, p1, Size{5, 5}, 0.0, 0, 360, {0, 255, 255}, 3);

            if (abs(y1 - y) < 20) {
                line(src_img, p, p1, cv::Scalar(255, 0, 0), 2, cv::LINE_8);
                line(src_img, p_2, p1_2, cv::Scalar(0, 255, 0), 2, cv::LINE_8);
            }
        }

        imshow("colour detection", src_img);
    }
    return 0;
}