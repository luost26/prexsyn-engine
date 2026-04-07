#include "rxn_lib_factory.hpp"

#include <cstddef>
#include <filesystem>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <csv.hpp>
#include <nlohmann/json.hpp>

#include "../chemistry/chemistry.hpp"
#include "../utility/logging.hpp"
#include "rxn_lib.hpp"

namespace prexsyn::chemspace {

static std::vector<std::string> split_line_by_tab(const std::string &line) {
    std::vector<std::string> parts;
    size_t start = 0;
    while (start < line.size()) {
        size_t end = line.find('\t', start);
        if (end == std::string::npos) {
            end = line.size();
        }
        parts.push_back(line.substr(start, end - start));
        start = end + 1;
    }
    return parts;
}

std::unique_ptr<ReactionLibrary> rxn_lib_from_plain_text(const std::filesystem::path &path,
                                                         bool ignore_errors) {
    auto rxn_lib = std::make_unique<ReactionLibrary>();

    std::ifstream infile(path);
    if (!infile) {
        throw std::runtime_error("Failed to open reaction library file: " + path.string());
    }

    logger()->info("Starting to load reactions from plain text file: {}", path.string());

    std::string line;
    unsigned line_no = 0;
    while (std::getline(infile, line)) {
        line_no++;
        line.erase(0, line.find_first_not_of(" \t\r\n"));
        line.erase(line.find_last_not_of(" \t\r\n") + 1);
        if (line.empty() || line.starts_with("#")) {
            continue; // Skip empty lines and comments
        }

        auto parts = split_line_by_tab(line);
        std::string smarts, name;
        if (parts.size() == 1) {
            smarts = std::move(parts[0]);
            name = "RXN_" + std::to_string(line_no);
        } else if (parts.size() >= 2) {
            smarts = std::move(parts[0]);
            name = std::move(parts[1]);
        }

        try {
            std::unique_ptr<Reaction> rxn;
            if (parts.size() <= 2) {
                rxn = Reaction::from_smarts(smarts);
            } else {
                std::vector<std::string> reactant_names(parts.begin() + 2, parts.end());
                rxn = Reaction::from_smarts(smarts, reactant_names);
            }
            rxn_lib->add({.reaction = std::move(rxn), .name = std::move(name)});
        } catch (const ReactionError &e) {
            logger()->warn("Line " + std::to_string(line_no) + ": " + e.what());
            if (!ignore_errors) {
                throw;
            }
            continue;
        }
    }

    logger()->info("Done. Loaded: {}", rxn_lib->size());

    return rxn_lib;
}

std::unique_ptr<ReactionLibrary> rxn_lib_from_csv(const std::filesystem::path &path,
                                                  const ReactionCSVConfig &config,
                                                  bool ignore_errors) {
    auto rxn_lib = std::make_unique<ReactionLibrary>();

    csv::CSVReader reader(path.string());
    logger()->info("Starting to load reactions from CSV: {}", path.string());
    size_t rowno = 0;
    for (auto &row : reader) {
        rowno++;
        try {
            std::string smarts, name;
            if (!row[config.smarts_column].try_get(smarts)) {
                logger()->warn("Missing SMARTS column at row {}", rowno);
                continue;
            }
            if (!row[config.name_column].try_get(name)) {
                name = "RXN_" + std::to_string(rowno);
            }

            std::string reactant_names_str;
            std::vector<std::string> reactant_names;
            if (row[config.reactant_name_column].try_get(reactant_names_str)) {
                size_t start = 0;
                while (start < reactant_names_str.size()) {
                    size_t end = reactant_names_str.find(config.reactant_name_delimiter, start);
                    if (end == std::string::npos) {
                        end = reactant_names_str.size();
                    }
                    reactant_names.push_back(reactant_names_str.substr(start, end - start));
                    start = end + config.reactant_name_delimiter.size();
                }
            }

            std::unique_ptr<Reaction> rxn;
            if (reactant_names.empty()) {
                rxn = Reaction::from_smarts(smarts);
            } else {
                rxn = Reaction::from_smarts(smarts, reactant_names);
            }
        } catch (const ReactionError &e) {
            logger()->warn("Row {}: {}", rowno, e.what());
            if (!ignore_errors) {
                throw;
            }
        }
    }

    return rxn_lib;
}

std::unique_ptr<ReactionLibrary> rxn_lib_from_json(const std::filesystem::path &path,
                                                   bool ignore_errors) {
    auto rxn_lib = std::make_unique<ReactionLibrary>();

    std::ifstream infile(path);
    if (!infile) {
        throw std::runtime_error("Failed to open reaction library file: " + path.string());
    }

    nlohmann::json data;
    try {
        infile >> data;
    } catch (const nlohmann::json::exception &e) {
        throw std::runtime_error("Failed to parse reaction JSON file: " + path.string() +
                                 ", error: " + e.what());
    }

    if (!data.is_array()) {
        throw std::runtime_error("Reaction JSON root must be an array: " + path.string());
    }

    logger()->info("Starting to load reactions from JSON: {}", path.string());

    size_t entry_no = 0;
    for (const auto &entry : data) {
        entry_no++;

        try {
            if (!entry.is_object()) {
                throw ReactionError("JSON entry must be an object");
            }

            std::string name = "RXN_" + std::to_string(entry_no);
            if (entry.contains("name")) {
                if (!entry.at("name").is_string()) {
                    throw ReactionError("Field 'name' must be a string");
                }
                name = entry.at("name").get<std::string>();
            }

            if (!entry.contains("reactants") || !entry.at("reactants").is_object()) {
                throw ReactionError("Field 'reactants' must be an object");
            }
            if (!entry.contains("product") || !entry.at("product").is_string()) {
                throw ReactionError("Field 'product' must be a string");
            }

            const auto &reactants = entry.at("reactants");
            if (reactants.empty()) {
                throw ReactionError("Field 'reactants' must not be empty");
            }

            std::vector<std::string> reactant_names;
            reactant_names.reserve(reactants.size());

            std::string reaction_smarts;
            bool first = true;
            for (auto it = reactants.begin(); it != reactants.end(); ++it) {
                if (!it.value().is_string()) {
                    throw ReactionError("Reactant SMARTS for '" + it.key() + "' must be a string");
                }
                if (!first) {
                    reaction_smarts += ".";
                }
                reaction_smarts += it.value().get<std::string>();
                reactant_names.push_back(it.key());
                first = false;
            }

            reaction_smarts += ">>";
            reaction_smarts += entry.at("product").get<std::string>();

            auto rxn = Reaction::from_smarts(reaction_smarts, reactant_names);
            rxn_lib->add({.reaction = std::move(rxn), .name = std::move(name)});
        } catch (const ReactionError &e) {
            logger()->warn("Entry {}: {}", entry_no, e.what());
            if (!ignore_errors) {
                throw;
            }
        }
    }

    logger()->info("Done. Loaded: {}", rxn_lib->size());
    return rxn_lib;
}

} // namespace prexsyn::chemspace
