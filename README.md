# ☄️ LIKE push-down for FSST


## Setup

Clone the repo, along with the submodules:

```
git clone --recurse-submodules git@github.com:utndatasystems/fsst-like.git
cd fsst-like
```

Build FSST:

```
./scripts/fsst-setup.sh
```

## Build

```
mkdir -p build
cd build
cmake ..
make
```

### Build `token-vldb2026`

From repo root:

```bash
cmake -S third-party/token-vldb2026 -B third-party/token-vldb2026/build \
  -DSGTT_TOKENIZER_RUNTIME_BIN=ON
cmake --build third-party/token-vldb2026/build -j20
```

## Benchmark

First, prepare the datasets.

```
./scripts/data-setup.sh 10
```

and run

```
./build/main data/l_comment_sf10.csv %special%
```
