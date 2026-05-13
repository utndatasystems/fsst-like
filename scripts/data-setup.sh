#!/usr/bin/env bash
set -e

scale_factor=${1:-1}

if [[ ! "${scale_factor}" =~ ^[0-9]+([.][0-9]+)?$ || "${scale_factor}" =~ ^0+([.]0+)?$ ]]; then
  echo "usage: $0 [scale_factor]" >&2
  echo "scale_factor must be a positive number, e.g. 1 or 10" >&2
  exit 1
fi

scale_suffix=${scale_factor//./_}

duckdb -c "
INSTALL tpch;
LOAD tpch;
CALL dbgen(sf=${scale_factor});

COPY (SELECT l_comment FROM lineitem) TO 'data/l_comment_sf${scale_suffix}.csv' (HEADER);
COPY (SELECT p_type FROM part) TO 'data/p_type_sf${scale_suffix}.csv' (HEADER);
"
