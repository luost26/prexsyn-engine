#include <boost/core/noncopyable.hpp>
#include <boost/python.hpp>
#include <boost/python/numpy/ndarray.hpp>
#include <random>

#include "../feature/pydict_builder.hpp"
#include "pipeline.hpp"
#include "pipeline_v2.hpp"

namespace py = boost::python;
namespace np = boost::python::numpy;
using namespace prexsyn_engine;

py::tuple data_pipeline_get(DataPipeline<PyDictBuilder> &pipeline) {
    auto result = pipeline.get();
    return py::make_tuple(std::get<0>(result), std::get<1>(result));
}

np::dtype get_numpy_dtype(DType::Type dtype) {
    switch (dtype) {
    case DType::Float:
        return np::dtype::get_builtin<float>();
    case DType::Long:
        return np::dtype::get_builtin<long>();
    case DType::Bool:
        return np::dtype::get_builtin<bool>();
    default:
        throw std::runtime_error("Unsupported dtype");
    }
}

template <size_t capacity>
py::dict data_pipeline_v2_get(DataPipelineV2<capacity> &pipeline,
                              size_t n = 1) {
    py::dict result;
    pipeline.get(
        [&](const std::vector<typename DataBuffer<
                capacity>::ReadTransaction::ReadEntry> &entries) {
            for (const auto &entry : entries) {
                py::list shape;
                shape.append(n);
                for (const auto &dim : entry.shape) {
                    shape.append(dim);
                }
                auto arr =
                    np::zeros(py::tuple(shape), get_numpy_dtype(entry.dtype));
                auto ptr = reinterpret_cast<std::byte *>(arr.get_data());
                std::memcpy(ptr, entry.span1.data(), entry.span1.size());
                if (entry.span2.has_value()) {
                    std::memcpy(ptr + entry.span1.size(), entry.span2->data(),
                                entry.span2->size());
                }
                result[entry.name] = arr;
            }
        },
        n);

    return result;
}

BOOST_PYTHON_MODULE(pipeline) {
    np::initialize();

    py::class_<DataPipeline<PyDictBuilder>, boost::noncopyable>(
        "DataPipeline",
        py::init<size_t, std::shared_ptr<ChemicalSpaceDefinition>,
                 SynthesisGeneratorOption, std::shared_ptr<FeaturizerSet>>(
            (py::arg("num_threads"), py::arg("csd"), py::arg("gen_option"),
             py::arg("featurizer"))))
        .def("start", &DataPipeline<PyDictBuilder>::start)
        .def("get", &data_pipeline_get)
        .def("stop", &DataPipeline<PyDictBuilder>::stop);

    py::class_<DataPipelineV2<8192>, boost::noncopyable>(
        "DataPipelineV2",
        py::init<size_t, std::shared_ptr<ChemicalSpaceDefinition>,
                 SynthesisGeneratorOption, std::shared_ptr<FeaturizerSet>,
                 std::mt19937::result_type>(
            (py::arg("num_threads"), py::arg("csd"), py::arg("gen_option"),
             py::arg("featurizer"),
             py::arg("base_seed") = std::random_device{}())))
        .def("start", &DataPipelineV2<8192>::start)
        .def("get", &data_pipeline_v2_get<8192>, (py::arg("n") = 1))
        .def("stop", &DataPipelineV2<8192>::stop);
}
