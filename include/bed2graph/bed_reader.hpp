#pragma once

#include <cstdint>
#include <filesystem>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace bed2graph {

struct BedRecord {
    std::string chromosome;
    std::int64_t start;
    std::int64_t end;
    std::string FRV;
};

struct DirectedGraph {
    std::set<std::string> nodes;
    std::set<std::pair<std::string, std::string>> edges;
};

std::vector<BedRecord> read_bed_file(const std::filesystem::path& path);
std::map<std::string, std::vector<std::string>> group_frvs_by_chromosome(
    const std::vector<BedRecord>& records);
std::map<std::string, std::vector<std::size_t>> find_frv_indices_by_chromosome(
    const std::map<std::string, std::vector<std::string>>& frvs_by_chromosome,
    std::uint64_t fr_number);
DirectedGraph build_graph(
    const std::map<std::string, std::vector<std::string>>& frvs_by_chromosome,
    const std::map<std::string, std::vector<std::size_t>>& matching_indices,
    std::uint64_t k);
void write_dot_file(const DirectedGraph& graph, const std::filesystem::path& path);

}  // namespace bed2graph
