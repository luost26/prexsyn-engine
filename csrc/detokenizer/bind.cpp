#include "bind.hpp"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <span>
#include <stdexcept>

#include <pybind11/cast.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <pybind11/detail/common.h>
#include <pybind11/detail/using_smart_holder.h>

#include "../chemspace/chemspace.hpp"
#include "../descriptor/descriptor.hpp"
#include "detokenizer.hpp"

namespace py = pybind11;
using namespace prexsyn;
using namespace detokenizer;

using TokenNumPyArray = py::array_t<std::int64_t, py::array::c_style | py::array::forcecast>;

static std::span<const std::int64_t> single_as_span(const TokenNumPyArray &tokens) {
    if (tokens.ndim() != 2 || tokens.shape(1) != 3) {
        throw std::invalid_argument("tokens must be a 2D array with shape (length, 3)");
    }
    return {tokens.data(), static_cast<size_t>(tokens.size())};
}

static std::span<const std::int64_t> batch_as_span(const TokenNumPyArray &tokens) {
    if (tokens.ndim() != 3 || tokens.shape(2) != 3) {
        throw std::invalid_argument("tokens must be a 3D array with shape (batch_size, length, 3)");
    }
    return {tokens.data(), static_cast<size_t>(tokens.size())};
}

void def_module_detokenizer(pybind11::module &m) {
    m.def(
        "detokenize",
        [](const TokenNumPyArray &tokens, const std::shared_ptr<chemspace::ChemicalSpace> &cs,
           const descriptor::TokenDef &token_def, std::optional<size_t> max_outcomes_per_reaction) {
            return detokenize(single_as_span(tokens), cs, token_def, max_outcomes_per_reaction);
        },
        py::arg("tokens"), py::arg("chemical_space"),
        py::arg("token_def") = descriptor::kDefaultTokenDef,
        py::arg("max_outcomes_per_reaction") = std::nullopt);

    py::class_<MultiThreadedDetokenizer, py::smart_holder>(m, "MultiThreadedDetokenizer")
        .def(py::init<const std::shared_ptr<chemspace::ChemicalSpace> &,
                      const descriptor::TokenDef &, std::optional<size_t>>(),
             py::arg("chemical_space"), py::arg("token_def") = descriptor::kDefaultTokenDef,
             py::arg("max_outcomes_per_reaction") = std::nullopt)
        .def(
            "__call__",
            [](const MultiThreadedDetokenizer &detok, size_t batch_size,
               const TokenNumPyArray &tokens) { return detok(batch_size, batch_as_span(tokens)); },
            py::arg("batch_size"), py::arg("tokens"));
}
