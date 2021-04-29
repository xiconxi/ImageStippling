#include "ImageSampler.h"

#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <random>

#define JCV_REAL_TYPE double
#define JC_VORONOI_IMPLEMENTATION
#include "jc_voronoi.h"

class CentroidVoronoiDiagram
{
public:
    explicit CentroidVoronoiDiagram(jcv_rect bbox, Eigen::MatrixXd& d)
        : density_(d), bbox_(bbox) {
    }

    ~CentroidVoronoiDiagram() { if(diagram_.internal) jcv_diagram_free(&diagram_); }

    void random_sampling(int n_sample) {
        if(n_sample == -1)
            n_sample = 128*128;
        points_.reserve(n_sample);
        std::random_device rd;  //Will be used to obtain a seed for the random number engine
        std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
        std::uniform_real_distribution<> dis(.0, 1.0);

        double d_scale = density_.maxCoeff();

        while(points_.size() < n_sample) {
            double y = dis(gen);
            double x = dis(gen);
            double d = dis(gen)*d_scale;
            if(density_.d(y, x) < d) continue;
            points_.push_back(jcv_point{x, y});
        }
    }

    jcv_point polygon_centroid(const jcv_graphedge* graph_edge);

    jcv_point varying_polygon_centroid(const jcv_graphedge* graph_edge, size_t n_sample = 5);

    double relax_points()
    {
        assert(!density_.is_empty());
        const jcv_site* sites = jcv_diagram_get_sites(&diagram_);
        double max_diff = 1e-9;
        for (int i = 0; i < diagram_.numsites; ++i)
        {
            const jcv_site* site = &sites[i];
            const jcv_point pre_p = points_[site->index];
            if (!density_.is_empty())
                points_[site->index] = varying_polygon_centroid(site->edges);
            else
                points_[site->index] = polygon_centroid(site->edges);

            max_diff = std::max(max_diff, jcv_point_dist_sq(&pre_p, &points_[site->index]) );
        }
        return max_diff;
    }

    double lloyd(){
        if(!points_.size()) random_sampling(-1);
        jcv_diagram_generate(points_.size(), points_.data(), &bbox_, 0, &diagram_);
        return relax_points();
    }

    std::string export_svg(std::string file);
private:
    jcv_diagram            diagram_{nullptr};
    std::vector<jcv_point> points_;
    ImageSampler           density_;
    jcv_rect               bbox_;
};

int main(int argc, char** argv)
{
    assert(argc == 3);
    jcv_rect bounding_box = {{0.0f, 0.0f}, {1.0f, 1.0f}};
    Eigen::MatrixXd img_density = read_density(argv[1]);
    img_density = (255 - img_density.array()).max(10).pow(1.2);
    std::cout << img_density.maxCoeff() << std::endl;
    CentroidVoronoiDiagram cvd( bounding_box, img_density);
//    cvd.random_sampling(5000);
    for(int i = 0; i < 500; i++) {
        double max_diff = cvd.lloyd();
        if((i+1)%4 == 0)
            cvd.export_svg(argv[2]+std::to_string(int(i/4))+".svg");
        if(max_diff < 3e-7 and i > 30) break;
        std::cout << "\r" << "lloyd: " << i << ' ' << max_diff << std::endl ;
    }
    cvd.export_svg(std::string(argv[2])+"-stippling.svg");
}

jcv_point CentroidVoronoiDiagram::polygon_centroid(const jcv_graphedge* graph_edge) {
    double total_det = 0;
    jcv_point center{0, 0};
    for(;graph_edge;graph_edge = graph_edge->next) {
        jcv_point p1 = graph_edge->pos[0], p2 = graph_edge->pos[1];
        double det = p1.x * p2.y - p2.x * p1.y;
        total_det += det;
        center.x += (p1.x + p2.x) * det;
        center.y += (p1.y + p2.y) * det;
    }
    center.x /= 3 * total_det;
    center.y /= 3 * total_det;
    return center;
}

jcv_point CentroidVoronoiDiagram::varying_polygon_centroid(const jcv_graphedge* graph_edge, size_t k) {
    jcv_point pc = polygon_centroid(graph_edge);
    jcv_point center{0, 0};
    double W = 0;
    for(;graph_edge;graph_edge = graph_edge->next) {
        jcv_point p1 = graph_edge->pos[0], p2 = graph_edge->pos[1];
        Eigen::RowVector3d centroid = density_.centroid(p1.x, p1.y, p2.x, p2.y, pc.x, pc.y, k);
        center.x += centroid.x();
        center.y += centroid.y();
        W += centroid.z();
    }
    center.x /= W;
    center.y /= W;
    return center;
}

std::string  CentroidVoronoiDiagram::export_svg(std::string output_file) {
    std::stringstream svg_str;
    double w = 512;
    double h = density_.is_empty() ? 512: density_.rows()/double(density_.cols())*512;
    svg_str << "<svg width=\"" << w << "\" height=\"" << h << "\" viewBox=\"0 0 " << w << ' ' << h
            << R"(" xmlns="http://www.w3.org/2000/svg" >)"
            << "<rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n";
    const jcv_site* sites = jcv_diagram_get_sites(&diagram_);
    for (size_t i = 0; i < diagram_.numsites && density_.is_empty(); i++) {
        svg_str << "<polygon points=\"";
        jcv_graphedge* graph_edge = sites[i].edges;
        while (graph_edge) {
            jcv_point p1 = graph_edge->pos[0], p2 = graph_edge->pos[1];
            svg_str << p1.x*w << ',' << p1.y*h << ' ' << p2.x*w << ',' << p2.y*h << ' ';
            graph_edge = graph_edge->next;
        }
        svg_str << R"(" fill="#1d1d9b" stroke="white" stroke-width="2" />)" << "\n";
    }

    for( int i = 0; i < diagram_.numsites; ++i ) {
        const jcv_site* site = &sites[i];
        jcv_point p = site->p, pc = polygon_centroid(site->edges);
        if(density_.is_empty()){
            svg_str << "<circle cx=\"" << pc.x * w  << "\" cy=\"" << pc.y*h << R"(" r="1" fill="red"/>)" << '\n';
            svg_str << "<circle cx=\"" << p.x * w  << "\" cy=\"" << p.y*h << R"(" r="1" fill="yellow"/>)" << '\n';
        }else
            svg_str << "<circle cx=\"" << p.x * w  << "\" cy=\"" << p.y*h << R"(" r="1" fill="black"/>)" << '\n';
    }

    svg_str << "</svg>";

    if(output_file.size() != 0) {
        std::ofstream svg_file(output_file);
        svg_file << svg_str.str();
        svg_file.close();
    }


    return svg_str.str();
}
