"""Resumable chunked execution helper for the PoseBusters CLI."""

from __future__ import annotations

import csv
import subprocess
from pathlib import Path
from typing import Optional


def _csv_row_count(path: Path) -> int:
    if not path.is_file() or path.stat().st_size == 0:
        return -1
    with path.open(newline="", encoding="utf-8") as handle:
        return sum(1 for _ in csv.DictReader(handle))


def _write_rows(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def run_posebusters_chunked(
    *,
    table_csv: Path,
    output_csv: Path,
    bust_bin: str,
    max_workers: int,
    posebusters_config: Optional[str],
    chunk_size: int,
) -> None:
    """Run PoseBusters in checkpointed chunks and merge results in input order."""
    with table_csv.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        input_fieldnames = list(reader.fieldnames or [])
        input_rows = list(reader)
    if not input_fieldnames or not input_rows:
        raise ValueError(f"PoseBusters input table is empty: {table_csv}")

    effective_chunk_size = len(input_rows) if chunk_size <= 0 else chunk_size
    parts_dir = output_csv.parent / f".{output_csv.stem}_parts"
    parts_dir.mkdir(parents=True, exist_ok=True)
    result_parts: list[Path] = []

    for part_index, start in enumerate(range(0, len(input_rows), effective_chunk_size)):
        chunk_rows = input_rows[start : start + effective_chunk_size]
        chunk_table = parts_dir / f"input_{part_index:04d}.csv"
        chunk_result = parts_dir / f"result_{part_index:04d}.csv"
        _write_rows(chunk_table, input_fieldnames, chunk_rows)

        if _csv_row_count(chunk_result) != len(chunk_rows):
            chunk_result.unlink(missing_ok=True)
            command = [
                bust_bin,
                "-t",
                str(chunk_table),
                "--outfmt",
                "csv",
                "--output",
                str(chunk_result),
                "--max-workers",
                str(max_workers),
            ]
            if posebusters_config:
                command.extend(["--config", posebusters_config])
            subprocess.run(command, check=True)
            actual_rows = _csv_row_count(chunk_result)
            if actual_rows != len(chunk_rows):
                raise RuntimeError(
                    f"PoseBusters chunk {part_index} returned {actual_rows} rows; "
                    f"expected {len(chunk_rows)}"
                )
        else:
            print(f"Resume: keeping completed PoseBusters chunk {chunk_result}")
        result_parts.append(chunk_result)

    merged_rows: list[dict[str, str]] = []
    result_fieldnames: list[str] = []
    for result_part in result_parts:
        with result_part.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            if not result_fieldnames:
                result_fieldnames = list(reader.fieldnames or [])
            elif list(reader.fieldnames or []) != result_fieldnames:
                raise RuntimeError(f"Inconsistent PoseBusters columns in {result_part}")
            merged_rows.extend(reader)

    if len(merged_rows) != len(input_rows):
        raise RuntimeError(
            f"Merged PoseBusters output has {len(merged_rows)} rows; "
            f"expected {len(input_rows)}"
        )
    _write_rows(output_csv, result_fieldnames, merged_rows)
