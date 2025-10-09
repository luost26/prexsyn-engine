#pragma once

#include <functional>

#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

#include "../utils/assert.hpp"
#include "builder.hpp"

namespace py = boost::python;
namespace np = boost::python::numpy;

namespace prexsyn_engine {
class PyDictBuilder : public FeatureBuilder {
    std::vector<std::function<void(py::dict &)>> callbacks;

  public:
    using feature_type = py::dict;

    PyDictBuilder() {};

    py::dict get() {
        py::dict dict;
        for (auto &callback : callbacks) {
            callback(dict);
        }
        callbacks.clear();
        return dict;
    }

#define DEFINE_ADD_METHOD(T)                                                   \
    void add(const std::string &key, const T &value) override {                \
        callbacks.push_back(                                                   \
            [key, value](py::dict &dict) { dict[key] = value; });              \
    }                                                                          \
                                                                               \
    void add(const std::string &key, const std::vector<T> &vec) override {     \
        callbacks.push_back([key, vec](py::dict &dict) {                       \
            const long shape[1] = {(long)vec.size()};                          \
            np::ndarray np_arr =                                               \
                np::zeros(1, shape, np::dtype::get_builtin<T>());              \
            for (size_t i = 0; i < vec.size(); ++i) {                          \
                np_arr[i] = vec[i];                                            \
            }                                                                  \
            dict[key] = np_arr;                                                \
        });                                                                    \
    }                                                                          \
                                                                               \
    void add(const std::string &key, const std::vector<std::vector<T>> &mat)   \
        override {                                                             \
        Ensures(mat.size() > 0);                                               \
        callbacks.push_back([key, mat](py::dict &dict) {                       \
            auto n_rows = mat.size();                                          \
            auto n_cols = mat[0].size();                                       \
            const long shape[2] = {(long)n_rows, (long)n_cols};                \
            np::ndarray np_arr =                                               \
                np::zeros(2, shape, np::dtype::get_builtin<T>());              \
            for (size_t i = 0; i < n_rows; ++i) {                              \
                Ensures(mat[i].size() == n_cols);                              \
                for (size_t j = 0; j < n_cols; ++j) {                          \
                    np_arr[i][j] = mat[i][j];                                  \
                }                                                              \
            }                                                                  \
            dict[key] = np_arr;                                                \
        });                                                                    \
    }

    DEFINE_ADD_METHOD(long);
    DEFINE_ADD_METHOD(float);
    DEFINE_ADD_METHOD(bool);
#undef DEFINE_ADD_METHOD
};

} // namespace prexsyn_engine
