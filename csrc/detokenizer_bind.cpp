#include "detokenizer.hpp"

#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <boost/python/numpy/dtype.hpp>
#include <boost/python/numpy/ndarray.hpp>

namespace py = boost::python;
namespace np = boost::python::numpy;
using namespace prexsyn_engine;

SynthesisVector detokenize_batch(const Detokenizer &detokenizer,
                                 const np::ndarray &token_types,
                                 const np::ndarray &bb_indices,
                                 const np::ndarray &rxn_indices) {
    Ensures(np::equivalent(token_types.get_dtype(),
                           np::dtype::get_builtin<TypeToken>()));
    Ensures(np::equivalent(bb_indices.get_dtype(),
                           np::dtype::get_builtin<BuildingBlockToken>()));
    Ensures(np::equivalent(rxn_indices.get_dtype(),
                           np::dtype::get_builtin<ReactionToken>()));
    Ensures(token_types.get_nd() == 2 && bb_indices.get_nd() == 2 &&
            rxn_indices.get_nd() == 2);
    for (size_t i = 0; i < 2; ++i)
        Ensures(token_types.shape(i) == bb_indices.shape(i) &&
                token_types.shape(i) == rxn_indices.shape(i));

    size_t batch_size = token_types.shape(0);
    size_t stride = token_types.shape(1);
    size_t numel = batch_size * stride;

    std::span<TypeToken> token_types_span(
        reinterpret_cast<TypeToken *>(token_types.get_data()), numel);
    std::span<BuildingBlockToken> bb_indices_span(
        reinterpret_cast<BuildingBlockToken *>(bb_indices.get_data()), numel);
    std::span<ReactionToken> rxn_indices_span(
        reinterpret_cast<ReactionToken *>(rxn_indices.get_data()), numel);
    return detokenizer.detokenize_many(token_types_span, bb_indices_span,
                                       rxn_indices_span, stride);
}

BOOST_PYTHON_MODULE(detokenizer) {
    np::initialize();

    py::class_<Detokenizer>(
        "Detokenizer",
        py::init<std::shared_ptr<BuildingBlockList>,
                 std::shared_ptr<ReactionList>, PostfixNotationTokenDef>(
            (py::arg("building_blocks"), py::arg("reactions"),
             py::arg("token_def") = PostfixNotationTokenDef())))
        .def("__call__", &detokenize_batch,
             (py::arg("token_types"), py::arg("bb_indices"),
              py::arg("rxn_indices")));
}
