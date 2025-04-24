#!/bin/bash

echo "Building FSST.."

(
  # Exit on failure.
  set -e
  cd fsst
  mkdir -p build
  cd build
  cmake ..
  make
)

echo "FSST is built."