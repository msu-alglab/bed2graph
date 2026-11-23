#include "bed2graph/bed_reader.hpp"

#include <charconv>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>

namespace {

std::uint64_t parse_non_negative_integer(std::string_view value, std::string_view name) {
    std::uint64_t parsed = 0;
    const auto result = std::from_chars(value.data(), value.data() + value.size(), parsed);
    if (value.empty() || result.ec != std::errc{} ||
        result.ptr != value.data() + value.size()) {
        throw std::invalid_argument(std::string(name) + " must be a non-negative integer");
    }
    return parsed;
}

}  // namespace

int main(int argc, char* argv[]) {
    if (argc != 4) {
        std::cerr << "Usage: bed2graph <input.bed> <FR> <k>\n";
        return 2;
    }

    try {
        const auto FR = parse_non_negative_integer(argv[2], "FR");
        const auto k = parse_non_negative_integer(argv[3], "k");
        const auto records = bed2graph::read_bed_file(argv[1]);
        const auto frvs_by_chromosome = bed2graph::group_frvs_by_chromosome(records);
        const auto matching_indices =
            bed2graph::find_frv_indices_by_chromosome(frvs_by_chromosome, FR);

        std::cout << "Read " << records.size() << " BED records (FR=" << FR
                  << ", k=" << k << ") across " << frvs_by_chromosome.size()
                  << " chromosomes.\n";
        const std::filesystem::path output_directory =
            std::filesystem::path(argv[1]).parent_path();
        for (const auto& [chromosome, frvs] : frvs_by_chromosome) {
            const std::filesystem::path chromosome_path(chromosome);
            if (chromosome_path.filename() != chromosome || chromosome == "." ||
                chromosome == "..") {
                throw std::invalid_argument("chromosome cannot be used as an output filename: " +
                                            chromosome);
            }

            const auto graph = bed2graph::build_graph(
                {{chromosome, frvs}}, {{chromosome, matching_indices.at(chromosome)}}, k);
            if (graph.nodes.empty()) {
                continue;
            }

            const auto dot_path = output_directory / (chromosome + ".dot");
            bed2graph::write_dot_file(graph, dot_path);
            std::cout << "Generated " << graph.nodes.size() << " nodes and "
                      << graph.edges.size()
                      << (graph.edges.size() == 1 ? " edge" : " edges") << " in "
                      << dot_path.string() << ".\n";
        }
        // for (const auto& [chromosome, indices] : matching_indices) {
        //     std::cout << chromosome << " matching FRV indices (zero-based):";
        //     for (const auto index : indices) {
        //         std::cout << ' ' << index;
        //     }
        //     std::cout << '\n';
        //}
    } catch (const std::invalid_argument& error) {
        std::cerr << "bed2graph: error: " << error.what() << '\n';
        return 2;
    } catch (const std::exception& error) {
        std::cerr << "bed2graph: error: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
