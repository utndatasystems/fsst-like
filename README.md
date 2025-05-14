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

## Benchmark

First, prepare the datasets.

```
./scripts/data-setup.sh
```

and run

```
./build/main data/l_comment.csv %special%
```