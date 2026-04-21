
#!/usr/bin/env bash

set -euo pipefail

mkdir -p data

rm -f data/static.csv data/dynamic_fixed.csv data/dynamic_variable.csv data/dynamic_decreasing.csv

julia --project=. compare.jl

cd forest_of_trees
cargo build --release
cargo run --release
cd ..

cd proposal_array
rm -rf build
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build .
./bench_sampling
cd ..
cd ..

cd bus
cd WSS-WIRS
rm -rf build
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build .
./bench_sampling
cd ..
cd ..
cd ..

mkdir -p figures

julia --project=. plot.jl
