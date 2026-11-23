#include "bed2graph/bed_reader.hpp"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <string>

namespace {

void require(bool condition, const std::string& message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

class TemporaryBedFile {
public:
    explicit TemporaryBedFile(const std::string& contents)
        : path_(std::filesystem::temp_directory_path() / "bed2graph_reader_test.bed") {
        std::ofstream output(path_);
        output << contents;
        if (!output) {
            throw std::runtime_error("could not create temporary BED test file");
        }
    }

    ~TemporaryBedFile() {
        std::error_code error;
        std::filesystem::remove(path_, error);
    }

    const std::filesystem::path& path() const {
        return path_;
    }

private:
    std::filesystem::path path_;
};

void test_reads_records_and_extra_fields() {
    TemporaryBedFile file(
        "# comment\n"
        "track name=example\n"
        "chr1\t10\t25\t1:7\n"
        "chrX\t30\t40\t2:4\n");

    const auto records = bed2graph::read_bed_file(file.path());
    require(records.size() == 2, "expected two BED records");
    require(records[0].chromosome == "chr1", "chromosome was not parsed");
    require(records[0].start == 10 && records[0].end == 25,
            "coordinates were not parsed");
    require(records[0].FRV == "1:7", "FRV was not preserved");
    require(records[1].FRV == "2:4", "FRV was not preserved");
}

void test_rejects_invalid_coordinates() {
    TemporaryBedFile file("chr1\t12\t8\n");
    try {
        static_cast<void>(bed2graph::read_bed_file(file.path()));
    } catch (const std::runtime_error& error) {
        require(std::string(error.what()).find(":1:") != std::string::npos,
                "error did not include the input line number");
        return;
    }
    throw std::runtime_error("expected invalid coordinates to be rejected");
}

void test_groups_frvs_by_chromosome_and_sorts_by_start() {
    const std::vector<bed2graph::BedRecord> records{
        {"chr1", 30, 35, "5:99"},
        {"chr2", 5, 10, "2:6"},
        {"chr1", 10, 15, "5:20"},
        {"chr1", 10, 12, "3:14"},
        {"chr1", 20, 25, "2:4"},
    };

    const auto grouped = bed2graph::group_frvs_by_chromosome(records);
    require(grouped.size() == 2, "expected records grouped into two chromosomes");
    require(grouped.at("chr1") ==
                std::vector<std::string>{"5:20", "3:14", "2:4", "5:99"},
            "FRVs were not sorted by start while preserving equal-start order");
    require(grouped.at("chr2") == std::vector<std::string>{"2:6"},
            "chromosome group contains unexpected FRVs");

    const auto matching_indices =
        bed2graph::find_frv_indices_by_chromosome(grouped, 5);
    require(matching_indices.at("chr1") == std::vector<std::size_t>{0, 3},
            "matching FR indices were not found in sorted-array positions");
    require(matching_indices.at("chr2").empty(),
            "unexpected FR match found for chromosome");
}

void test_rejects_malformed_frv() {
    const std::map<std::string, std::vector<std::string>> grouped{
        {"chr1", {"not-an-frv"}},
    };
    try {
        static_cast<void>(bed2graph::find_frv_indices_by_chromosome(grouped, 1));
    } catch (const std::runtime_error& error) {
        require(std::string(error.what()).find("chr1") != std::string::npos,
                "FRV error did not include chromosome");
        return;
    }
    throw std::runtime_error("expected malformed FRV to be rejected");
}

void test_builds_graph_with_trimmed_windows_and_writes_dot() {
    const std::map<std::string, std::vector<std::string>> frvs_by_chromosome{
        {"chr1", {"1:0", "2:0", "3:0", "4:0", "1:1"}},
        {"chr2", {"2:1", "1:2", "1:3", "3:1"}},
    };
    const std::map<std::string, std::vector<std::size_t>> matching_indices{
        {"chr1", {0, 4}},
        {"chr2", {1, 2}},
    };

    const auto graph = bed2graph::build_graph(frvs_by_chromosome, matching_indices, 1);
    require(graph.nodes ==
                std::set<std::string>{
                    "1:0", "2:0", "4:0", "1:1", "2:1", "1:2", "1:3", "3:1"},
            "graph nodes do not reflect clamped matching windows");
    require(graph.edges ==
                std::set<std::pair<std::string, std::string>>{
                    {"1:0", "2:0"}, {"4:0", "1:1"}, {"2:1", "1:2"},
                    {"1:2", "1:3"}, {"1:3", "3:1"}},
            "graph edges do not connect consecutive FRVs within each window");

    const auto dot_path =
        std::filesystem::temp_directory_path() / "bed2graph_graph_test.dot";
    bed2graph::write_dot_file(graph, dot_path);
    std::ifstream input(dot_path);
    require(static_cast<bool>(input), "DOT file was not created");
    const std::string dot((std::istreambuf_iterator<char>(input)),
                          std::istreambuf_iterator<char>());
    std::error_code error;
    std::filesystem::remove(dot_path, error);
    require(dot.find("digraph bed2graph {") != std::string::npos,
            "DOT graph declaration is missing");
    require(dot.find("\"1:0\" -> \"2:0\";") != std::string::npos,
            "DOT directed edge is missing");
}

}  // namespace

int main() {
    try {
        test_reads_records_and_extra_fields();
        test_rejects_invalid_coordinates();
        test_groups_frvs_by_chromosome_and_sorts_by_start();
        test_rejects_malformed_frv();
        test_builds_graph_with_trimmed_windows_and_writes_dot();
    } catch (const std::exception& error) {
        std::cerr << "Test failure: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
