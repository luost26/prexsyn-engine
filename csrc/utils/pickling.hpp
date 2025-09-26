#pragma once

#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

namespace synthesis_backend {

namespace py = boost::python;
namespace np = boost::python::numpy;

template <class T> struct pickle_suite : py::pickle_suite {
    static boost::python::tuple getinitargs(const T &) {
        return boost::python::make_tuple();
    }
    static boost::python::tuple getstate(const T &obj) {
        auto pickle_str = obj.pickle();
        auto pickle_np = np::zeros(py::make_tuple(pickle_str.size()),
                                   np::dtype::get_builtin<char>());
        std::memcpy(pickle_np.get_data(), pickle_str.c_str(),
                    pickle_str.size());
        return boost::python::make_tuple(pickle_np);
    }
    static void setstate(T &obj, const boost::python::tuple &state) {
        np::ndarray pickle_np = boost::python::extract<np::ndarray>(state[0]);
        std::string pickle_str(pickle_np.get_data(), pickle_np.shape(0));
        std::unique_ptr<T> new_obj(T::unpickle(pickle_str));
        obj = *new_obj;
    }
};

} // namespace synthesis_backend
