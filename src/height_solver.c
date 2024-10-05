#include "height_solver.h"

struct ParabolicError
{
    static constexpr bool ComputesJacobian = false; // Set to true to indicate that the functor computes the Jacobian.

    template<typename Scalar, int Inputs, int Outputs>
    void operator()(const Eigen::Matrix<Scalar, Inputs, 1> &xval, Eigen::Matrix<Scalar, Outputs, 1> &fval) const {
        //fval.resize(xval.size());
        fval.resize(1);

        // set the heights of the mesh based on the input vector
        for(lsqcpp::Index i = 0; i < xval.size(); ++i) {
            caustic_design.mesh->source_points[i][2] = xval(i);
        }

        Scalar E_int = 0.0f;
        for (lsqcpp::Index i = 0; i < xval.size(); ++i) {
            // Calculate the current normal
            std::vector<double> normal = caustic_design.calculate_vertex_normal(caustic_design.mesh->source_points, i);

            // Get the target normal
            std::vector<double> normal_trg = { caustic_design.normals[0][i], caustic_design.normals[1][i], caustic_design.normals[2][i]};

            normal_trg = normalize(normal_trg);
            normal = normalize(normal);

            // Calculate the difference
            std::vector<double> diff = vector_subtract(normal, normal_trg);

            //fval(i) = diff[0];

            // Sum up the squared components of the difference
            double energy = diff[0] * diff[0];// + diff[1] * diff[1] + diff[2] * diff[2];

            E_int += energy;
        }

        fval(0) = E_int;
    }
};