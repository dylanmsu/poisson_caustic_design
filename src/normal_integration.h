#include "mesh.h"
#include <vector>

#include "ceres/ceres.h"
#include "costFunctor.h"

using ceres::AutoDiffCostFunction;
using ceres::NumericDiffCostFunction;
using ceres::CENTRAL;
using ceres::CostFunction;
using ceres::Problem;
using ceres::Solver;
using ceres::Solve;

class normal_integration
{
private:
    void gatherVertexInformation(Mesh &mesh, uint vertexIndex, std::vector<int> &neighborList, std::vector<int> &neighborMap, std::vector<int> & gridNeighbors);
    void addResidualBlocks(Mesh &mesh, Problem *problem, uint vertexIndex, std::vector<int> &neighbors, std::vector<int> &neighborMap, std::vector<int> &gridNeighbors, double *vertices, std::vector<double> &trg_normal);

    std::vector<std::vector<int> > neighborsPerVertex;
    std::vector<std::vector<int> > neighborMapPerVertex;
    std::vector<std::vector<int> > eightNeighborsPerVertex;
    std::vector<std::vector<double>> x_sources;

    double* vertices;
public:
    normal_integration(/* args */);
    ~normal_integration();

    void initialize_data(Mesh &mesh);

    void perform_normal_integration(Mesh &mesh, std::vector<std::vector<double>> &desired_normals);
};