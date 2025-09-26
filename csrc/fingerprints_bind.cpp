#include "fingerprints.hpp"

#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <boost/python/numpy/dtype.hpp>
#include <boost/python/numpy/ndarray.hpp>
#include <boost/python/tuple.hpp>
#include <omp.h>

#include "synthesis.hpp"
#include "types.hpp"
#include "utils/assert.hpp"

namespace py = boost::python;
namespace np = boost::python::numpy;
using namespace synthesis_backend;

template <std::vector<float> F(const std::optional<Mol_sptr> &)>
np::ndarray fp_func_numpy_wrapper(const Mol_sptr &mol) {
    if (!mol) {
        throw std::invalid_argument("Input molecule is None.");
    }
    std::vector<float> fp = F(mol);
    np::ndarray result = np::empty(py::make_tuple((int)fp.size()),
                                   np::dtype::get_builtin<float>());
    std::copy(fp.begin(), fp.end(),
              reinterpret_cast<float *>(result.get_data()));
    return result;
}

np::ndarray fp_func_bind(const Mol_sptr &mol, const std::string &name) {
    auto fp = fp_func<float>(name)(mol);
    np::ndarray result = np::empty(py::make_tuple((int)fp.size()),
                                   np::dtype::get_builtin<float>());
    std::copy(fp.begin(), fp.end(),
              reinterpret_cast<float *>(result.get_data()));
    return result;
}

np::ndarray get_fingerprints(const py::list &mlist,
                             const std::string &fp_type) {
    auto n = py::len(mlist);
    if (n == 0) {
        throw std::invalid_argument("Input molecule list is empty.");
    }
    auto fp0 = fp_func<float>(fp_type)(py::extract<Mol_sptr>(mlist[0]));
    auto fp_dim = fp0.size();
    np::ndarray result = np::zeros(py::make_tuple(n, (int)fp_dim),
                                   np::dtype::get_builtin<float>());

    std::vector<Mol_sptr> mols;
    mols.reserve(n);
    for (int i = 0; i < n; ++i) {
        mols.push_back(py::extract<Mol_sptr>(mlist[i]));
    }

    auto base_ptr = reinterpret_cast<float *>(result.get_data());
#pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        auto ptr = base_ptr + i * fp_dim;
        if (!mols[i]) {
            continue;
        }
        auto fp = fp_func<float>(fp_type)(mols[i]);
        std::copy(fp.begin(), fp.end(), ptr);
    }

    return result;
}

float diversity(const py::list &mlist, const std::string &fp_type) {
    auto n = py::len(mlist);
    std::vector<Mol_sptr> mols;
    mols.reserve(n);
    for (int i = 0; i < n; ++i) {
        auto mol = py::extract<Mol_sptr>(mlist[i])();
        if (mol) {
            mols.push_back(mol);
        }
    }
    if (mols.size() == 0) {
        throw std::invalid_argument("No valid molecules found.");
    }

    float total_sim = 0;
#pragma omp parallel for reduction(+ : total_sim)
    for (int i = 0; i < (int)mols.size(); ++i) {
        for (int j = i + 1; j < (int)mols.size(); ++j) {
            total_sim += tanimoto_similarity(mols[i], mols[j], fp_type);
        }
    }

    int count = (int)mols.size() * ((int)mols.size() - 1) / 2;
    return 1 - total_sim / count;
}

np::ndarray pairwise_tanimoto_similarity(py::list mlist1, py::list mlist2,
                                         const std::string &fp_type) {
    auto n1 = py::len(mlist1), n2 = py::len(mlist2);
    np::ndarray result =
        np::zeros(py::make_tuple(n1, n2), np::dtype::get_builtin<float>());
    std::vector<Mol_sptr> mols1, mols2;
    mols1.reserve(n1);
    mols2.reserve(n2);
    for (int i = 0; i < n1; ++i) {
        mols1.push_back(py::extract<Mol_sptr>(mlist1[i]));
    }
    for (int j = 0; j < n2; ++j) {
        mols2.push_back(py::extract<Mol_sptr>(mlist2[j]));
    }

#pragma omp parallel for
    for (int ij = 0; ij < n1 * n2; ++ij) {
        auto ptr = reinterpret_cast<float *>(result.get_data()) + ij;
        int i = ij / n2, j = ij % n2;
        const Mol_sptr &mol1 = mols1[i];
        const Mol_sptr &mol2 = mols2[j];
        if (!mol1 || !mol2) {
            *ptr = 0;
            continue;
        }
        *ptr = tanimoto_similarity(mol1, mol2, fp_type);
    }

    return result;
}

np::ndarray mol_to_syntheses_tanimoto_similarity(const Mol_sptr &mol,
                                                 const SynthesisVector &syns,
                                                 const std::string &fp_type) {
    size_t max_num_products = 0;
    for (const auto &syn : syns) {
        if (syn->stack_size() == 0)
            continue;
        max_num_products = std::max(max_num_products, syn->top().size());
    }

    np::ndarray result =
        np::zeros(py::make_tuple((int)syns.size(), (int)max_num_products),
                  np::dtype::get_builtin<float>());

#pragma omp parallel for
    for (int i = 0; i < (int)syns.size(); ++i) {
        const auto &syn = syns[i];
        if (syn->stack_size() == 0) {
            continue;
        }
        auto ptr =
            reinterpret_cast<float *>(result.get_data()) + i * max_num_products;
        size_t index = 0;
        for (const auto &product : syn->top()) {
            Ensures(product != nullptr);
            auto sim_this = tanimoto_similarity(mol, product, fp_type);
            if (index < max_num_products) {
                ptr[index] = sim_this;
                ++index;
            }
        }
    }
    return result;
}

BOOST_PYTHON_MODULE(fingerprints) {
    np::initialize();

    py::def("ecfp4", &fp_func_numpy_wrapper<ecfp4_fingerprint<float>>,
            (py::arg("mol")));
    py::def("fp_func", &fp_func_bind, (py::arg("mol"), py::arg("fp_type")));
    py::def("get_fingerprints", &get_fingerprints,
            (py::arg("mlist"), py::arg("fp_type")));
    py::def("diversity", &diversity, (py::arg("mlist"), py::arg("fp_type")));
    py::def("tanimoto_similarity", &tanimoto_similarity,
            (py::arg("mol1"), py::arg("mol2"), py::arg("fp_type") = "ecfp4"));
    py::def(
        "pairwise_tanimoto_similarity", &pairwise_tanimoto_similarity,
        (py::arg("mlist1"), py::arg("mlist2"), py::arg("fp_type") = "ecfp4"));
    py::def("mol_to_syntheses_tanimoto_similarity",
            &mol_to_syntheses_tanimoto_similarity,
            (py::arg("mol"), py::arg("syns"), py::arg("fp_type") = "ecfp4"));
}
