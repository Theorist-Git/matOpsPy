#include <pybind11/pybind11.h>
#include "../include/matOps.hpp"
#include <pybind11/stl.h>       // Automatically converts STL containers (like std::vector) to Python types
#include <pybind11/operators.h> // Enables operator overloading support
#include <sstream> 

namespace py = pybind11;

PYBIND11_MODULE(pyMatOps, m) {
    m.doc() = "Python bindings for matOps library";

    py::class_<Matrix>(m, "Matrix")
        // Constructors
        .def(py::init<const std::vector<std::vector<double>> &>(),
             py::arg("container"),
             "Construct a Matrix from a 2D vector.")
        .def(py::init<size_t, size_t, double>(),
             py::arg("rows"), py::arg("cols"), py::arg("initialValue"),
             "Construct a Matrix with given dimensions and an initial value.")

        // Methods
        .def("shape", &Matrix::shape, "Returns the (rows, cols) of the matrix.")
        .def("transpose", &Matrix::transpose, "Returns the transposed matrix.")
        .def("determinant", &Matrix::determinant, "Computes the determinant of the matrix.")
        .def("inverse", &Matrix::inverse, "Computes the inverse of the matrix.")
        .def("insertRow",
             (Matrix (Matrix::*)(std::vector<double>, size_t) const) &Matrix::insertRow,
             py::arg("row"), py::arg("idx"),
             "Insert a new row (as vector) at index idx.")
        .def("insertRow",
             (Matrix (Matrix::*)(double, size_t) const) &Matrix::insertRow,
             py::arg("rowVal"), py::arg("idx"),
             "Insert a new row filled with rowVal at index idx.")
        .def("insertCol",
             (Matrix (Matrix::*)(std::vector<double>, size_t) const) &Matrix::insertCol,
             py::arg("col"), py::arg("idx"),
             "Insert a new column (as vector) at index idx.")
        .def("insertCol",
             (Matrix (Matrix::*)(double, size_t) const) &Matrix::insertCol,
             py::arg("colVal"), py::arg("idx"),
             "Insert a new column filled with colVal at index idx.")
        .def("hStack", &Matrix::hStack,
             "Horizontally concatenates two matrices (same number of rows).")
        .def("vStack", &Matrix::vStack,
             "Vertically concatenates two matrices (same number of columns).")
        .def("extractMatrix", &Matrix::extractMatrix,
             py::arg("rowSlice"), py::arg("colSlice"),
             "Extracts a submatrix given row and column index ranges (as pairs).")

        // Support for element access using the [] operator:
        .def("__getitem__", [](Matrix &self, py::tuple index) -> double {
            if (index.size() != 2)
                throw std::runtime_error("Index must be a tuple of two elements (row, col)");
            size_t i = index[0].cast<size_t>();
            size_t j = index[1].cast<size_t>();
            return self(i, j);
        })
        .def("__setitem__", [](Matrix &self, py::tuple index, double value) {
            if (index.size() != 2)
                throw std::runtime_error("Index must be a tuple of two elements (row, col)");
            size_t i = index[0].cast<size_t>();
            size_t j = index[1].cast<size_t>();
            self(i, j) = value;
        })

        // Operator overloading for arithmetic operations
        .def(py::self + py::self)
        .def(py::self + double())
        .def(double() + py::self)
        .def(py::self - py::self)
        .def(py::self - double())
        .def(double() - py::self)
        .def(py::self * py::self)
        .def(py::self * double())
        .def(double() * py::self)
        .def(py::self / double())
        .def(py::self == py::self)

        // String representations using the friend operator<<
        .def("__str__", [](const Matrix &m) {
            std::stringstream ss;
            ss << m;
            return ss.str();
        })
        .def("__repr__", [](const Matrix &m) {
            std::stringstream ss;
            ss << "Matrix(" << m << ")";
            return ss.str();
        });
}