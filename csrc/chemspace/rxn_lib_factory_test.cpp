#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>

#include <gtest/gtest.h>

#include "rxn_lib_factory.hpp"

namespace {

std::filesystem::path find_project_root() {
    auto current = std::filesystem::absolute(__FILE__);
    while (current.has_parent_path()) {
        current = current.parent_path();
        if (std::filesystem::exists(current / "resources/test/reaction.json")) {
            return current;
        }
    }
    throw std::runtime_error("Could not locate project root from __FILE__");
}

std::filesystem::path write_temp_json(const std::string &contents, const std::string &filename) {
    auto path = std::filesystem::temp_directory_path() / filename;
    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("Failed to open temp file: " + path.string());
    }
    out << contents;
    out.close();
    return path;
}

} // namespace

TEST(RxnLibFactoryTest, JsonLoaderLoadsExampleReactionFile) {
    const auto path = find_project_root() / "resources/test/reaction.json";
    auto rxn_lib = prexsyn::chemspace::rxn_lib_from_json(path);

    ASSERT_NE(rxn_lib, nullptr);
    EXPECT_EQ(rxn_lib->size(), 3U);

    const auto &suzuki = rxn_lib->get("Suzuki Coupling");
    EXPECT_EQ(suzuki.reaction->num_reactants(), 2U);
    EXPECT_TRUE(suzuki.reaction->reactant_name_to_index().contains("Halides"));
    EXPECT_TRUE(suzuki.reaction->reactant_name_to_index().contains("Boronates"));

    const auto &amide = rxn_lib->get("Amide Coupling");
    EXPECT_EQ(amide.reaction->num_reactants(), 2U);
    EXPECT_TRUE(amide.reaction->reactant_name_to_index().contains("Acids"));
    EXPECT_TRUE(amide.reaction->reactant_name_to_index().contains("Amines"));

    const auto &alkynylation = rxn_lib->get("Alkynylation");
    EXPECT_EQ(alkynylation.reaction->num_reactants(), 1U);
    EXPECT_TRUE(alkynylation.reaction->reactant_name_to_index().contains("Halides"));
}

TEST(RxnLibFactoryTest, JsonLoaderThrowsOnInvalidEntryWhenIgnoreErrorsIsFalse) {
    const auto path = write_temp_json(
        R"json([
  {
	"name": "Bad",
	"reactants": {"A": "[C:1]-[Br]"}
  }
])json",
        "prexsyn_rxn_lib_invalid_throw.json");

    EXPECT_THROW(prexsyn::chemspace::rxn_lib_from_json(path, false), std::exception);
    std::filesystem::remove(path);
}

TEST(RxnLibFactoryTest, JsonLoaderSkipsInvalidEntriesWhenIgnoreErrorsIsTrue) {
    const auto path = write_temp_json(
        R"json([
  {
	"name": "Bad",
	"reactants": {"A": "[C:1]-[Br]"}
  },
  {
	"name": "Good",
	"reactants": {"A": "[C:1]-[Br]"},
	"product": "[C:1]-[C]"
  }
])json",
        "prexsyn_rxn_lib_invalid_skip.json");

    auto rxn_lib = prexsyn::chemspace::rxn_lib_from_json(path, true);

    ASSERT_NE(rxn_lib, nullptr);
    EXPECT_EQ(rxn_lib->size(), 1U);
    EXPECT_NO_THROW((void)rxn_lib->get("Good"));
    EXPECT_THROW((void)rxn_lib->get("Bad"), std::exception);

    std::filesystem::remove(path);
}
