//
// Created by pupa on 12/17/20.
//

#include "ImageSampler.h"

#define STB_IMAGE_IMPLEMENTATION
#include <stb_image.h>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>
#include <stdint.h>

#include <iostream>
using Eigen::Vector2d;

double inline clip(double x, double l, double r) {
    return (x > l)? (x < r ? x : r):l;
}

Eigen::MatrixXd read_density(std::string filename) {
    Eigen::MatrixXd image;

    int width_, height_, n_;
    if (filename.empty()) return image;
    std::uint8_t* img = stbi_load(filename.c_str(), &width_, &height_, &n_, 1);
    if ( img == nullptr) return image;

    image.resize(height_, width_);
    for (int j = 0; j < height_; j++)
        for (int i = 0; i < width_; i++)
            image(j, i) = img[j * width_ + i];

    stbi_image_free(img);
    return  image;
}


Eigen::MatrixX2d dither_sampling(Eigen::MatrixXd& img, int n_sample) {
    if(n_sample == -1) n_sample = 128*128;
    Eigen::MatrixXd f = img/img.sum()*n_sample;
    std::vector<Eigen::RowVector2i> samples_;

    for(auto r = 1; r < f.rows() -1; r++) {
        for(auto c = 1; c < f.cols() -1; c++) {
            if(std::round(f(r, c)) == 1)
                samples_.push_back(Eigen::RowVector2i(r, c));

            double error = f(r, c) - std::round(f(r, c));
            f(r, c+1) += error * 7/16.0;
            f(r+1, c-1) += error * 3/16.0;
            f(r+1, c) += error * 5/16.0;
            f(r+1, c+1) += error * 1/16.0;
        }
    }
    Eigen::MatrixX2d samples(samples_.size(), 2);
    for(size_t i = 0; i < samples_.size(); i++) {
        samples(i, 0) = double(samples_[i][0])/img.rows();
        samples(i, 1) = double(samples_[i][1])/img.cols();
    }
    return samples;
}

double ImageSampler::d(double y, double x) const {
    x = (cols() - 1) * clip(x, 0, 1);
    y = (rows() - 1) * clip(y, 0, 1);

    size_t floor_x = std::floor(x), floor_y = std::floor(y);
    size_t ceil_x = std::ceil(x), ceil_y = std::ceil(y);

    double _x = x - floor_x, _y = y - floor_y;

    double floor_pixel = pixel_(floor_y, floor_x) * (1 - _x) + pixel_(ceil_y, floor_x) * _x;
    double ceil_pixel = pixel_(floor_y, ceil_x) * (1 - _x) + pixel_(ceil_y, ceil_x) * _x;
    return (floor_pixel * (1 - _y) + ceil_pixel * _y);
}

double inline area_2(double x1, double y1, double x2, double y2, double x3, double y3) {
    return abs((x1 - x3) * (y2 - y3) - (x2 - x3) * (y1 - y3));
}


Eigen::Vector3d ImageSampler::centroid(double x1, double y1, double x2, double y2, double x3, double y3, size_t k) {
    auto [max_x, min_x] = std::minmax({x1, x2, x3});
    auto [max_y, min_y] = std::minmax({y1, y2, y3});
    double area = area_2(x1, y1, x2, y2, x3, y3);
    double dx = (max_x - min_x)/k;
    double dy = (max_y - min_y)/k;
    double integral_x = 0, integral_y = 0;
    double total_d = 0;
    for(size_t i = 0; i <= k; i++ ) {
        for(size_t j = 0; j <= k; j++) {
            double x = min_x + i * dx;
            double y = min_y + j * dy;
            double alpha = area_2(x, y, x2, y2, x3, y3)/area;
            double beta = area_2(x1, y1, x, y, x3, y3)/area;
            double gamma = (1-alpha -beta);
            if( 0 <= alpha && alpha <= 1 &&  beta >= 0 && beta <= 1.0 && gamma >= 0 && gamma <= 1.0 ){
                double d = std::max(0.001, this->d(y, x));
                integral_x += d * x;
                integral_y += d * y;
                total_d += d;
            }
        }
    }
    return {integral_x, integral_y, total_d};
}
