#include "bed2graph/bed_reader.hpp"

#include <algorithm>
#include <charconv>
#include <fstream>
#include <map>
#include <ostream>
#include <stdexcept>
#include <string>
#include <string_view>

namespace bed2graph {
namespace {

std::vector<std::string_view> split_tabs(const std::string& line) {
    std::vector<std::string_view> fields;
    std::size_t field_start = 0;

    while (true) {
        const auto tab = line.find('\t', field_start);
        if (tab == std::string::npos) {
            fields.emplace_back(line.data() + field_start, line.size() - field_start);
            return fields;
        }
        fields.emplace_back(line.data() + field_start, tab - field_start);
        field_start = tab + 1;
    }
}

std::int64_t parse_coordinate(std::string_view value, const std::filesystem::path& path,
                              std::size_t line_number, std::string_view field_name) {
    std::int64_t coordinate = 0;
    const auto result = std::from_chars(value.data(), value.data() + value.size(), coordinate);
    if (result.ec != std::errc{} || result.ptr != value.data() + value.size()) {
        throw std::runtime_error(path.string() + ":" + std::to_string(line_number) +
                                 ": invalid " + std::string(field_name) + " coordinate");
    }
    return coordinate;
}

bool is_ignored_line(std::string_view line) {
    const auto first = line.find_first_not_of(" \t");
    if (first == std::string_view::npos || line[first] == '#') {
        return true;
    }

    line.remove_prefix(first);
    return line.starts_with("track ") || line.starts_with("track\t") ||
           line == "track" || line.starts_with("browser ") ||
           line.starts_with("browser\t") || line == "browser";
}

std::uint64_t parse_frv_integer(std::string_view value, const std::string& chromosome,
                                std::size_t index) {
    std::uint64_t parsed = 0;
    const auto result = std::from_chars(value.data(), value.data() + value.size(), parsed);
    if (value.empty() || result.ec != std::errc{} ||
        result.ptr != value.data() + value.size()) {
        throw std::runtime_error("invalid FRV for chromosome " + chromosome +
                                 " at index " + std::to_string(index) +
                                 ": expected non-negative integers in FR:variant format");
    }
    return parsed;
}

void write_quoted_dot_string(std::ostream& output, std::string_view value) {
    output << '"';
    for (const char character : value) {
        switch (character) {
            case '"':
                output << "\\\"";
                break;
            case '\\':
                output << "\\\\";
                break;
            case '\n':
                output << "\\n";
                break;
            case '\r':
                output << "\\r";
                break;
            case '\t':
                output << "\\t";
                break;
            default:
                output << character;
                break;
        }
    }
    output << '"';
}

}  // namespace

std::vector<BedRecord> read_bed_file(const std::filesystem::path& path) {
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("could not open BED file: " + path.string());
    }

    std::vector<BedRecord> records;
    std::string line;
    std::size_t line_number = 0;
    while (std::getline(input, line)) {
        ++line_number;
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        if (is_ignored_line(line)) {
            continue;
        }

        const auto fields = split_tabs(line);
        if (fields.size() < 3) {
            throw std::runtime_error(path.string() + ":" + std::to_string(line_number) +
                                     ": expected at least 3 tab-separated BED fields");
        }
        if (fields[0].empty()) {
            throw std::runtime_error(path.string() + ":" + std::to_string(line_number) +
                                     ": chromosome field must not be empty");
        }

        const auto start = parse_coordinate(fields[1], path, line_number, "start");
        const auto end = parse_coordinate(fields[2], path, line_number, "end");
        if (start < 0 || end <= start) {
            throw std::runtime_error(path.string() + ":" + std::to_string(line_number) +
                                     ": BED coordinates must satisfy 0 <= start < end");
        }

        BedRecord record{std::string(fields[0]), start, end, {}};
        for (std::size_t i = 3; i < fields.size(); ++i) {
            if (i > 3) {
                record.FRV.push_back('\t');
            }
            record.FRV.append(fields[i]);
        }
        records.emplace_back(std::move(record));
    }

    if (input.bad()) {
        throw std::runtime_error("error reading BED file: " + path.string());
    }
    return records;
}

std::map<std::string, std::vector<std::string>> group_frvs_by_chromosome(
    const std::vector<BedRecord>& records) {
    std::map<std::string, std::vector<const BedRecord*>> records_by_chromosome;
    for (const auto& record : records) {
        records_by_chromosome[record.chromosome].push_back(&record);
    }

    std::map<std::string, std::vector<std::string>> frvs_by_chromosome;
    for (auto& [chromosome, chromosome_records] : records_by_chromosome) {
        std::stable_sort(chromosome_records.begin(), chromosome_records.end(),
                         [](const BedRecord* left, const BedRecord* right) {
                             return left->start < right->start;
                         });

        auto& frvs = frvs_by_chromosome[chromosome];
        frvs.reserve(chromosome_records.size());
        for (const auto* record : chromosome_records) {
            frvs.push_back(record->FRV);
        }
    }
    return frvs_by_chromosome;
}

std::map<std::string, std::vector<std::size_t>> find_frv_indices_by_chromosome(
    const std::map<std::string, std::vector<std::string>>& frvs_by_chromosome,
    std::uint64_t fr_number) {
    std::map<std::string, std::vector<std::size_t>> indices_by_chromosome;
    for (const auto& [chromosome, frvs] : frvs_by_chromosome) {
        auto& matching_indices = indices_by_chromosome[chromosome];
        for (std::size_t index = 0; index < frvs.size(); ++index) {
            const std::string_view frv = frvs[index];
            const auto colon = frv.find(':');
            if (colon == std::string_view::npos ||
                frv.find(':', colon + 1) != std::string_view::npos) {
                throw std::runtime_error("invalid FRV for chromosome " + chromosome +
                                         " at index " + std::to_string(index) +
                                         ": expected FR:variant format");
            }

            const auto record_fr_number =
                parse_frv_integer(frv.substr(0, colon), chromosome, index);
            static_cast<void>(parse_frv_integer(frv.substr(colon + 1), chromosome, index));
            if (record_fr_number == fr_number) {
                matching_indices.push_back(index);
            }
        }
    }
    return indices_by_chromosome;
}

DirectedGraph build_graph(
    const std::map<std::string, std::vector<std::string>>& frvs_by_chromosome,
    const std::map<std::string, std::vector<std::size_t>>& matching_indices,
    std::uint64_t k) {
    DirectedGraph graph;
    for (const auto& [chromosome, indices] : matching_indices) {
        const auto chromosome_frvs = frvs_by_chromosome.find(chromosome);
        if (chromosome_frvs == frvs_by_chromosome.end() && !indices.empty()) {
            throw std::invalid_argument("matching indexes refer to unknown chromosome: " +
                                        chromosome);
        }
        if (chromosome_frvs == frvs_by_chromosome.end()) {
            continue;
        }
        const auto& frvs = chromosome_frvs->second;
        if (frvs.empty()) {
            if (!indices.empty()) {
                throw std::invalid_argument("matching indexes refer to an empty FRV array for " +
                                            chromosome);
            }
            continue;
        }

        for (const auto index : indices) {
            if (index >= frvs.size()) {
                throw std::invalid_argument("matching FRV index is out of range for chromosome " +
                                            chromosome);
            }

            const auto first = k >= index ? 0 : index - static_cast<std::size_t>(k);
            const auto last_distance = frvs.size() - 1 - index;
            const auto last = k >= last_distance
                                  ? frvs.size() - 1
                                  : index + static_cast<std::size_t>(k);

            for (auto current = first;; ++current) {
                graph.nodes.insert(frvs[current]);
                if (current < last) {
                    graph.edges.emplace(frvs[current], frvs[current + 1]);
                } else {
                    break;
                }
            }
        }
    }
    return graph;
}

void write_dot_file(const DirectedGraph& graph, const std::filesystem::path& path) {
    std::ofstream output(path);
    if (!output) {
        throw std::runtime_error("could not create DOT file: " + path.string());
    }

    output << "digraph bed2graph {\n";
    for (const auto& node : graph.nodes) {
        output << "  ";
        write_quoted_dot_string(output, node);
        output << ";\n";
    }
    for (const auto& [source, target] : graph.edges) {
        output << "  ";
        write_quoted_dot_string(output, source);
        output << " -> ";
        write_quoted_dot_string(output, target);
        output << ";\n";
    }
    output << "}\n";
    if (!output) {
        throw std::runtime_error("error writing DOT file: " + path.string());
    }
}

}  // namespace bed2graph
