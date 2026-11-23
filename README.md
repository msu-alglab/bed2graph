# bed2graph

`bed2graph` is a C++ command-line program that reads BED records containing
FRV values and writes a directed graph for each chromosome in Graphviz DOT
format.

## Requirements

- CMake 3.16 or newer
- A C++20-compatible compiler

## Build

From the project directory:

```sh
cmake -S . -B build
cmake --build build
```

## Input format

Each data line is tab-separated and has at least four fields:

```text
chromosome  start  end  FRV
```

The first three fields follow BED conventions: `start` and `end` are
non-negative integers with `start < end`. The fourth field is an FRV in
`FR:variant` form, where both values are non-negative integers (for example,
`2:17`). The FRV must be the fourth and final field; additional fields are
joined to it and will fail FRV validation.

Blank lines, lines starting with `#`, and BED `track` and `browser` lines are
ignored. Coordinates are used to sort records within each chromosome; equal
start coordinates retain their input order.

## Run

```sh
./build/bed2graph <input.bed> <FR> <k>
```

`FR` and `k` must be non-negative integers. For each chromosome, the program
finds every FRV whose number before the colon equals `FR`. For each match at
index `i` in the chromosome's start-sorted FRV list, it selects the range
`i-k` through `i+k`, clipped to the start and end of the list. Each selected
FRV becomes a node, and each adjacent pair in a selected range becomes a
directed edge. Duplicate nodes and edges are written once per chromosome.

The program writes `<chromosome>.dot` beside the input BED file for each
chromosome with at least one selected node. The output is a Graphviz directed
graph and can be rendered with Graphviz, for example:

```sh
dot -Tpng chr1.dot -o chr1.png
```

## Example

The [sample BED file](examples/sample.bed) contains two chromosomes. Build and
run the example with `FR=1` and `k=1`:

```sh
cmake -S . -B build
cmake --build build
./build/bed2graph examples/sample.bed 1 1
```

The command reports 10 records across two chromosomes and creates
`examples/chr1.dot` and `examples/chr2.dot`. The generated files are included
as [chr1 output](examples/chr1.dot) and [chr2 output](examples/chr2.dot).
Each matching FRV is included with its immediate neighbors in start order.
For example, `chr1` contains matches at indexes 0, 3, and 5, so the resulting
graph includes the clipped windows around all three matches.

## Tests

```sh
ctest --test-dir build --output-on-failure
```

## License

This project is distributed under the BSD 2-Clause License. See [LICENSE](LICENSE)
for the full text.
