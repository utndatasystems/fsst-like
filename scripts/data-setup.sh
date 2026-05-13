#!/usr/bin/env bash
set -e

mkdir -p data

duckdb -c "
INSTALL tpch;
LOAD tpch;
CALL dbgen(sf=1);

COPY (SELECT l_comment FROM lineitem) TO 'data/l_comment.csv' (HEADER FALSE);
COPY (SELECT p_type FROM part) TO 'data/p_type.csv' (HEADER FALSE);
"