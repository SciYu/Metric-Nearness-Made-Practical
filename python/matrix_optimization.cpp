#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <algorithm>
#include <cmath>
#include <vector>

namespace py = pybind11;

using Matrix = py::array_t<double>;

Matrix heuristic_improve(Matrix X0, Matrix D, int n_improve=1) {
    auto X_buf = X0.request();
    auto D_buf = D.request();
    if (X_buf.ndim != 2 || D_buf.ndim != 2 || X_buf.shape[0] != X_buf.shape[1]) {
        throw std::runtime_error("Input matrices must be square.");
    }

    int n = D_buf.shape[0];
    auto X = Matrix(X0);  // Copy of X0
    auto X_data = static_cast<double*>(X.mutable_data());
    auto D_data = static_cast<const double*>(D.data());

    for (int iter = 0; iter < n_improve; ++iter) {
        std::vector<bool> F(n * n);
        for (int i = 0; i < n * n; ++i) {
            F[i] = (X_data[i] > D_data[i]);
        }
        
        for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                double d = D_data[i * n + j];
                X_data[i * n + j] = d;
                X_data[j * n + i] = d;

                if (F[i * n + j]) {
                    double max_val = 0.0;
                    for (int k = 0; k < n; ++k) {
                        max_val = std::max(max_val, std::abs(X_data[k * n + i] - X_data[k * n + j]));
                    }
                    X_data[i * n + j] = max_val;
                } else {
                    double min_val = std::numeric_limits<double>::infinity();
                    for (int k = 0; k < n; ++k) {
                        min_val = std::min(min_val, X_data[k * n + i] + X_data[k * n + j]);
                    }
                    X_data[i * n + j] = min_val;
                }
                X_data[j * n + i] = X_data[i * n + j];
            }
        }
    }
    return X;
}

Matrix hlwb_projection(Matrix X0, Matrix D, int n_projection=100) {
    auto X_buf = X0.request();
    auto D_buf = D.request();
    if (X_buf.ndim != 2 || D_buf.ndim != 2 || X_buf.shape[0] != X_buf.shape[1]) {
        throw std::runtime_error("Input matrices must be square.");
    }

    int n = D_buf.shape[0];
    double t = 2.0;
    double lambda = 1.0 / t;
    auto X = Matrix(X0);  // Copy of X0
    auto X_data = static_cast<double*>(X.mutable_data());
    auto D_data = static_cast<const double*>(D.data());

    for (int iter = 1; iter <= n_projection; ++iter) {
        for (int i = 0; i < n * n; ++i) {
            X_data[i] = lambda * D_data[i] + (1 - lambda) * X_data[i];
        }
        lambda = 0.382 / iter;

        for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                for (int k = 0; k < n; ++k) {
                    if (k != i && k != j) {
                        double delta = (X_data[i * n + j] - X_data[i * n + k] - X_data[j * n + k]) / 3;
                        if (delta > 0) {
                            X_data[i * n + j] -= delta;
                            X_data[j * n + i] = X_data[i * n + j];
                            X_data[i * n + k] += delta;
                            X_data[k * n + i] = X_data[i * n + k];
                            X_data[j * n + k] += delta;
                            X_data[k * n + j] = X_data[j * n + k];
                        }
                    }
                }
            }
        }
    }

    return X;
}

PYBIND11_MODULE(matrix_optimization, m) {
    m.def("heuristic_improve", &heuristic_improve, "Heuristic improvement function",
          py::arg("X0"), py::arg("D"), py::arg("n_improve") = 1);
    m.def("hlwb_projection", &hlwb_projection, "HLWB projection function",
          py::arg("X0"), py::arg("D"), py::arg("n_projection") = 100);
}
