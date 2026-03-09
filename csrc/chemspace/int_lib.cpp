#include "int_lib.hpp"

#include <cstddef>
#include <istream>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>
#include <utility>

#include "../utility/serialization.hpp"

namespace prexsyn::chemspace {

std::unique_ptr<IntermediateLibrary> IntermediateLibrary::deserialize(std::istream &data) {
    boost::archive::binary_iarchive ia(data);
    size_t num_items = 0;
    ia >> num_items;
    auto int_lib = std::make_unique<IntermediateLibrary>();
    int_lib->intermediates_.reserve(num_items);
    for (size_t i = 0; i < num_items; ++i) {
        IntermediateItem item;
        ia >> item;
        int_lib->intermediates_.push_back(std::move(item));
    }
    ia >> int_lib->identifier_to_index_;
    return int_lib;
}

void IntermediateLibrary::serialize(std::ostream &stream) const {
    boost::archive::binary_oarchive oa(stream);
    oa << intermediates_.size();
    for (const auto &item : intermediates_) {
        oa << item;
    }
    oa << identifier_to_index_;
}

const IntermediateItem &IntermediateLibrary::get(Index index) const {
    if (index >= intermediates_.size()) {
        throw std::out_of_range("Intermediate index out of range");
    }
    return intermediates_[index];
}

const IntermediateItem &IntermediateLibrary::get(const std::string &identifier) const {
    auto it = identifier_to_index_.find(identifier);
    if (it == identifier_to_index_.end()) {
        throw std::out_of_range("Intermediate identifier not found: " + identifier);
    }
    return intermediates_[it->second];
}

IntermediateLibrary::Index IntermediateLibrary::add(const IntermediateEntry &entry) {
    if (identifier_to_index_.contains(entry.identifier)) {
        throw std::invalid_argument("Intermediate with the same identifier already exists: " +
                                    entry.identifier);
    }
    auto new_index = intermediates_.size();
    intermediates_.push_back(IntermediateItem{entry, new_index});
    identifier_to_index_[entry.identifier] = new_index;
    return new_index;
}

void IntermediateLibrary::clear() {
    intermediates_.clear();
    identifier_to_index_.clear();
}

} // namespace prexsyn::chemspace
