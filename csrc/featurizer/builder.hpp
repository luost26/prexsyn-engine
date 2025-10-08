#pragma once

#include <string>
#include <vector>

#include "dtype.hpp"

namespace prexsyn_engine {
class FeatureBuilder {
  public:
    virtual ~FeatureBuilder() = default;
#define DEFINE_ADD_METHODS(T)                                                  \
    virtual void add(const std::string &, const T &) = 0;                      \
    virtual void add(const std::string &, const std::vector<T> &) = 0;         \
    virtual void add(const std::string &,                                      \
                     const std::vector<std::vector<T>> &) = 0;
    DEFINE_ADD_METHODS(Long)
    DEFINE_ADD_METHODS(Float)
    DEFINE_ADD_METHODS(Bool)
#undef DEFINE_ADD_METHODS
    FeatureBuilder &erase_type() { return *this; }
};

} // namespace prexsyn_engine
