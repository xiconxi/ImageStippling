//
// Created by pupa on 12/17/20.
//

#pragma once

#include <string>
#include <vector>
#include <cstdint>
#include <iostream>

#include <Eigen/Core>

Eigen::MatrixXd read_density(std::string filename);

//! floyd steinberg dithering with point budget
Eigen::MatrixX2d dither_sampling(Eigen::MatrixXd& img, int n_sample);

class ImageSampler: public Eigen::MatrixXd
{
public:
    explicit ImageSampler(const Eigen::MatrixXd& d): Eigen::MatrixXd(d){}

    double d(double y, double x) const ;

    bool is_empty() const { return !this->size(); }


    Eigen::Vector3d centroid(double x1, double y1, double x2, double y2, double x3, double y3, size_t k = 3);

private:
    double pixel_(size_t i, size_t j) const {
        return this->operator()(i, j);
    }
};
