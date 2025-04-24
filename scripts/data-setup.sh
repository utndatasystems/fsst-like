#!/bin/bash
set -e

# Create a temporary directory for DuckDB data.
TEMP_DIR=$(mktemp -d)
DB_PATH="$TEMP_DIR/tpch.duckdb"  # Temporary DB path

DUCKDB_CLI="$HOME/.duckdb/cli/latest/duckdb"

OUTPUT_DIR="data"
mkdir -p "$OUTPUT_DIR"

# Columns to export.
EXPORTS=(
  "lineitem:l_comment"
  "part:p_type"
  # Add more..
)

# Prepare TPC-H data.
echo "Generating TPC-H (sf=1) data.."
"$DUCKDB_CLI" "$DB_PATH" <<SQL
INSTALL tpch;
LOAD tpch;
CALL dbgen(sf = 1);
SQL

# And export to csv.
echo "⚙️  Exporting.."

for entry in "${EXPORTS[@]}"; do
  IFS=":" read -r table column <<< "$entry"
  outfile="${OUTPUT_DIR}/${column}.csv"
  
  echo "Exporting $table.$column to $outfile.."

  "$DUCKDB_CLI" "$DB_PATH" <<SQL
COPY (
  SELECT $column FROM $table
) TO '$outfile' (HEADER, DELIMITER ',');
SQL
done

# Clean up by removing the temporary database
rm -rf "$TEMP_DIR"

echo "✅ Data setup: done."