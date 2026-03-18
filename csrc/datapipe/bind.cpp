#include "bind.hpp"

#include <cstddef>
#include <map>
#include <memory>
#include <optional>
#include <span>
#include <string>
#include <vector>

#include <pybind11/cast.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>

#include "../chemspace/chemspace.hpp"
#include "../descriptor/descriptor.hpp"
#include "../utility/data_type_bind.hpp"
#include "datapipe.hpp"

namespace py = pybind11;
using namespace prexsyn;
using namespace prexsyn::datapipe;

static auto get_from_data_pipeline(DataPipeline &pipeline, size_t batch_size) {
    const auto &schema = pipeline.buffer().schema();
    NamedReadBatch read{.batch_size = batch_size, .destinations = {}};
    py::dict data;
    for (const auto &col_def : schema) {
        const auto &name = col_def.name();
        const auto &dtype = data_type_to_numpy_dtype(col_def.dtype());
        std::vector<size_t> shape = col_def.shape();
        shape.insert(shape.begin(), batch_size);
        auto arr = py::array(dtype, shape);
        data[name.c_str()] = arr;
        auto span =
            std::span<std::byte>(reinterpret_cast<std::byte *>(arr.mutable_data()), arr.nbytes());
        read.destinations[name] = span;
    }
    pipeline.get(read);
    return data;
}

void def_module_datapipe(pybind11::module &m) {

    py::class_<GeneratorConfig, py::smart_holder>(m, "GeneratorConfig")
        .def(py::init<>())
        .def_readwrite("max_building_blocks", &GeneratorConfig::max_building_blocks)
        .def_readwrite("heavy_atom_limit", &GeneratorConfig::heavy_atom_limit);

    py::class_<Generator, py::smart_holder>(m, "Generator")
        .def(py::init<std::shared_ptr<chemspace::ChemicalSpace>, const Generator::Config &,
                      std::optional<size_t>>(),
             py::arg("chemical_space"), py::arg("config") = Generator::default_config,
             py::arg("random_seed") = std::nullopt)
        .def("next", &Generator::next)
        .def("next_with_product", &Generator::next_with_product);

    py::class_<DataPipeline, py::smart_holder>(m, "DataPipeline")
        .def(py::init<
                 const std::shared_ptr<chemspace::ChemicalSpace> &,
                 const std::map<std::string, std::shared_ptr<descriptor::MoleculeDescriptor>> &,
                 const std::map<std::string, std::shared_ptr<descriptor::SynthesisDescriptor>> &,
                 const GeneratorConfig &>(),
             py::arg("chemical_space"), py::arg("molecule_descriptors"),
             py::arg("synthesis_descriptors"),
             py::arg("generator_config") = Generator::default_config)
        .def("start_workers", &DataPipeline::start_workers)
        .def("stop_workers", &DataPipeline::stop_workers)
        .def("get", &get_from_data_pipeline);
}
