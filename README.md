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

## Benchmark

First, prepare the datasets.

```
./scripts/data-setup.sh
```

and run

```
make
./main
```