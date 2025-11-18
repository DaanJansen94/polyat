from __future__ import annotations

import argparse
import gzip
import json
import sys
from collections import defaultdict
from html import escape
from pathlib import Path
from statistics import median
from typing import Iterable

HISTOGRAM_MIN_LENGTH = 10

HISTOGRAM_MIN_LENGTH = 10


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="polyat",
        description=(
            "Quantify poly-A/T stretches (>=10/15/20 nt) across FASTQ reads "
            "and summarize counts per sample."
        ),
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Directory containing .fastq/.fastq.gz/.fq/.fq.gz files.",
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Directory where the summary table will be written.",
    )
    return parser.parse_args(argv)


def resolve_directory(path_str: str, role: str) -> Path:
    path = Path(path_str).expanduser().resolve()
    if role == "input":
        if not path.exists():
            sys.exit(f"[error] Input path does not exist: {path}")
        if not path.is_dir():
            sys.exit(f"[error] Input path is not a directory: {path}")
    else:
        path.mkdir(parents=True, exist_ok=True)
    return path


def find_fastq_files(input_dir: Path) -> list[Path]:
    fastq_files: list[Path] = []
    for entry in sorted(input_dir.iterdir()):
        if entry.is_file() and has_fastq_suffix(entry.name):
            fastq_files.append(entry)
    return fastq_files


def has_fastq_suffix(filename: str) -> bool:
    valid_suffixes = (".fastq", ".fastq.gz", ".fq", ".fq.gz")
    return any(filename.endswith(suffix) for suffix in valid_suffixes)


def open_fastq(path: Path) -> Iterable[str]:
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return open(path, "r", encoding="utf-8")


def sanitize_sample_name(file_path: Path) -> str:
    name = file_path.name
    for suffix in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return name


def count_poly_runs(
    file_path: Path,
) -> tuple[int, int, int, int, dict[int, int], list[tuple[int, int]]]:
    total = 0
    poly10 = 0
    poly15 = 0
    poly20 = 0
    histogram: dict[int, int] = defaultdict(int)
    positions: list[int] = []
    with open_fastq(file_path) as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            seq = handle.readline().strip()
            handle.readline()  # +
            handle.readline()  # quality
            if not seq:
                continue
            total += 1
            longest, start_idx = longest_poly_run(seq)
            if longest >= 10:
                poly10 += 1
            if longest >= 15:
                poly15 += 1
            if longest >= 20:
                poly20 += 1
            if longest >= HISTOGRAM_MIN_LENGTH:
                histogram[longest] += 1
                if start_idx is not None:
                    end_offset = len(seq) - (start_idx + longest)
                    positions.append(min(start_idx, end_offset))
    return total, poly10, poly15, poly20, dict(histogram), positions


def longest_poly_run(sequence: str) -> tuple[int, int | None]:
    longest = 0
    longest_start: int | None = None
    current = 0
    prev_char = ""
    current_start: int | None = None
    for index, base in enumerate(sequence.upper()):
        if base not in ("A", "T"):
            current = 0
            prev_char = ""
            current_start = None
            continue
        if base == prev_char:
            current += 1
        else:
            current = 1
            prev_char = base
            current_start = index
        if current > longest:
            longest = current
            longest_start = current_start
    return longest, longest_start


def format_percent(count: int, total: int) -> str:
    if total == 0:
        return "0.00"
    return f"{(count * 100) / total:.2f}"


def build_filterable_table(
    table_id: str, title: str, headers: list[tuple[str, str]], rows: list[list[str]]
) -> str:
    parts: list[str] = [
        "<section>",
        f"<h2>{escape(title)}</h2>",
        f"<table data-filterable='true' id='{escape(table_id)}'>",
        "<thead>",
        "<tr>",
    ]
    for label, _type in headers:
        parts.append(f"<th>{escape(label)}</th>")
    parts.append("</tr>")
    parts.append("<tr class='filters'>")
    for idx, (_label, col_type) in enumerate(headers):
        if col_type == "number":
            placeholder = "min value"
            input_type = "number"
        else:
            placeholder = "text"
            input_type = "text"
        parts.append(
            "<th>"
            f"<input data-col='{idx}' data-type='{col_type}' "
            f"type='{input_type}' placeholder='{placeholder}' />"
            "</th>"
        )
    parts.append("</tr></thead><tbody>")
    for row in rows:
        parts.append("<tr>")
        parts.extend(f"<td>{escape(value)}</td>" for value in row)
        parts.append("</tr>")
    parts.extend(["</tbody></table>", "</section>"])
    return "\n".join(parts)


def write_html_summary(
    output_dir: Path,
    summary_headers: list[tuple[str, str]],
    summary_rows: list[list[str]],
    histogram_series: dict[str, list[list[int]]],
    sample_order: list[str],
    position_headers: list[tuple[str, str]],
    position_rows: list[list[str]],
) -> None:
    output_file = output_dir / "polyA_report.html"
    css = (
        "body{font-family:Arial,sans-serif;margin:20px;background:#fefefe;}"
        "h1{margin-bottom:0.5em;}"
        "section{margin-bottom:30px;}"
        "#histogram-section label{display:block;margin-bottom:6px;font-weight:bold;}"
        "#histogram-section select{padding:6px 10px;margin-bottom:12px;}"
        "#histogram-canvas{width:100%;max-width:900px;height:auto;border:1px solid #ddd;background:#fff;}"
        ".hist-tooltip{position:absolute;padding:6px 10px;background:rgba(0,0,0,0.75);"
        "color:#fff;border-radius:4px;font-size:12px;pointer-events:none;"
        "transform:translate(-50%,-120%);white-space:nowrap;display:none;z-index:10;}"
        "table{border-collapse:collapse;width:100%;font-family:Arial,sans-serif;}"
        "th,td{border:1px solid #ccc;padding:8px;text-align:center;}"
        "th{background-color:#f4f4f4;}"
        ".filters th{background-color:#fafafa;}"
        ".filters input{width:100%;box-sizing:border-box;padding:4px;}"
        "tr:nth-child(even){background:#fafafa;}"
    )
    histogram_json = json.dumps(histogram_series)
    sample_json = json.dumps(sample_order)
    script = (
        "document.querySelectorAll('table[data-filterable]').forEach(table=>{"
        "const columnFilters=table.querySelectorAll('thead input[data-col]');"
        "const rows=table.querySelectorAll('tbody tr');"
        "function applyFilters(){"
        "rows.forEach(row=>{"
        "let visible=true;"
        "columnFilters.forEach(input=>{"
        "if(!visible)return;"
        "const value=input.value.trim();"
        "if(!value)return;"
        "const col=parseInt(input.dataset.col,10);"
        "const type=input.dataset.type;"
        "const cell=row.children[col];"
        "if(!cell)return;"
        "const cellText=cell.innerText.trim();"
        "if(type==='number'){"
        "const cellValue=parseFloat(cellText);"
        "const filterValue=parseFloat(value);"
        "if(isNaN(filterValue)||isNaN(cellValue))return;"
        "if(cellValue<filterValue){visible=false;}"
        "}else{"
        "if(!cellText.toLowerCase().includes(value.toLowerCase())){"
        "visible=false;"
        "}"
        "}"
        "});"
        "row.style.display=visible?'':'none';"
        "});"
        "}"
        "columnFilters.forEach(input=>input.addEventListener('input',applyFilters));"
        "applyFilters();"
        "});"
        f"const histogramData = {histogram_json};"
        f"const histogramSamples = {sample_json};"
        "const sampleSelect = document.getElementById('histogram-sample');"
        "const canvas = document.getElementById('histogram-canvas');"
        "const tooltip = document.getElementById('histogram-tooltip');"
        "const ctx = canvas ? canvas.getContext('2d') : null;"
        "let histogramBars = [];"
        "function showTooltip(bar, event){"
        "if(!tooltip || !bar)return;"
        "tooltip.textContent = `Length: ${bar.length} nt | Reads: ${bar.count}`;"
        "tooltip.style.display='block';"
        "tooltip.style.left = `${event.clientX}px`;"
        "tooltip.style.top = `${event.clientY}px`;"
        "}"
        "function hideTooltip(){"
        "if(tooltip){tooltip.style.display='none';}"
        "}"
        "function renderHistogram(sample){"
        "if(!ctx)return;"
        "ctx.clearRect(0,0,canvas.width,canvas.height);"
        "histogramBars = [];"
        "const entries = histogramData[sample] || [];"
        "const margin = 60;"
        "const width = canvas.width - margin * 2;"
        "const height = canvas.height - margin * 2;"
        "ctx.strokeStyle = '#333';"
        "ctx.fillStyle = '#333';"
        "ctx.lineWidth = 1;"
        "ctx.font = '12px Arial';"
        "ctx.beginPath();"
        "ctx.moveTo(margin, margin);"
        "ctx.lineTo(margin, margin + height);"
        "ctx.lineTo(margin + width, margin + height);"
        "ctx.stroke();"
        "if(entries.length === 0){"
        f"ctx.fillText('No data for selected sample (no runs >={HISTOGRAM_MIN_LENGTH} nt).', margin, margin);"
        "return;"
        "}"
        "const maxCount = Math.max(...entries.map(item => item[1]));"
        "if(maxCount === 0){"
        "ctx.fillText('All counts are zero for this sample.', margin, margin);"
        "return;"
        "}"
        "const tickCount = 5;"
        "const tickStep = maxCount / tickCount;"
        "ctx.fillStyle = '#000';"
        "ctx.textAlign = 'right';"
        "ctx.textBaseline = 'middle';"
        "for(let i=0;i<=tickCount;i++){"
        "const value = Math.round(i * tickStep);"
        "const y = margin + height - (value / maxCount) * height;"
        "ctx.beginPath();"
        "ctx.moveTo(margin - 5, y);"
        "ctx.lineTo(margin, y);"
        "ctx.stroke();"
        "ctx.fillText(String(value), margin - 8, y);"
        "ctx.strokeStyle='#eee';"
        "ctx.beginPath();"
        "ctx.moveTo(margin, y);"
        "ctx.lineTo(margin + width, y);"
        "ctx.stroke();"
        "ctx.strokeStyle='#333';"
        "}"
        "ctx.textAlign = 'center';"
        "ctx.textBaseline = 'middle';"
        "ctx.fillText('Count', margin - 35, margin - 20);"
        "ctx.fillText('Run length (nt)', margin + width / 2, margin + height + 45);"
        "const barWidth = width / entries.length;"
        "entries.forEach((item, index) => {"
        "const length = item[0];"
        "const count = item[1];"
        "const barHeight = (count / maxCount) * height;"
        "const x = margin + index * barWidth + barWidth * 0.1;"
        "const y = margin + height - barHeight;"
        "const w = barWidth * 0.8;"
        "const h = barHeight || (count > 0 ? 1 : 0);"
        "ctx.fillStyle = '#4a90e2';"
        "ctx.fillRect(x, y, w, h);"
        "histogramBars.push({x,y,width:w,height:h,count,length});"
        "if(entries.length <= 40){"
        "ctx.save();"
        "ctx.translate(x + w / 2, margin + height + 15);"
        "ctx.rotate(-Math.PI / 4);"
        "ctx.fillStyle = '#000';"
        "ctx.fillText(String(length), 0, 0);"
        "ctx.restore();"
        "}"
        "});"
        "}"
        "function initHistogram(){"
        "if(!sampleSelect || !canvas || !ctx)return;"
        "sampleSelect.innerHTML='';"
        "histogramSamples.forEach(sample => {"
        "const option=document.createElement('option');"
        "option.value=sample;"
        "option.textContent=sample;"
        "sampleSelect.appendChild(option);"
        "});"
        "if(histogramSamples.length>0){"
        "sampleSelect.value = histogramSamples[0];"
        "renderHistogram(sampleSelect.value);"
        "}else if(ctx){"
        "ctx.font='16px Arial';"
        "ctx.fillText('No histogram data available', 10, 30);"
        "}"
        "sampleSelect.addEventListener('change', () => {"
        "hideTooltip();"
        "renderHistogram(sampleSelect.value);"
        "});"
        "canvas.addEventListener('mousemove', (event) => {"
        "const rect = canvas.getBoundingClientRect();"
        "const x = event.clientX - rect.left;"
        "const y = event.clientY - rect.top;"
        "const bar = histogramBars.find(b => x>=b.x && x<=b.x+b.width && y>=b.y && y<=b.y+b.height);"
        "if(bar && bar.count>0){"
        "showTooltip(bar, event);"
        "}else{"
        "hideTooltip();"
        "}"
        "});"
        "canvas.addEventListener('mouseleave', hideTooltip);"
        "}"
        "initHistogram();"
    )
    parts = [
        "<!DOCTYPE html>",
        "<html lang='en'>",
        "<head>",
        "<meta charset='utf-8' />",
        "<title>polyA Counts</title>",
        f"<style>{css}</style>",
        "</head>",
        "<body>",
        "<h1>polyA/T Report</h1>",
        build_filterable_table(
            "polyat-summary-table", "polyA/T Summary", summary_headers, summary_rows
        ),
        "<section id='histogram-section'>",
        f"<h2>polyA/T Histogram (>={HISTOGRAM_MIN_LENGTH} nt)</h2>",
        "<label for='histogram-sample'>Sample</label>",
        "<select id='histogram-sample'></select>",
        "<canvas id='histogram-canvas' width='900' height='400'></canvas>",
        "<div id='histogram-tooltip' class='hist-tooltip'></div>",
        "</section>",
        build_filterable_table(
            "polyat-position-table",
            "polyA/T Position Offsets",
            position_headers,
            position_rows,
        ),
        f"<script>{script}</script>",
        "</body>",
        "</html>",
    ]
    output_file.write_text("\n".join(parts), encoding="utf-8")


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    input_dir = resolve_directory(args.input, "input")
    output_dir = resolve_directory(args.output, "output")
    output_file = output_dir / "polyA_counts.txt"

    fastq_files = find_fastq_files(input_dir)
    if not fastq_files:
        sys.exit(f"[error] No FASTQ/FASTQ.GZ files found in {input_dir}")

    summary_headers = [
        ("Sample", "text"),
        ("Total_Reads", "number"),
        ("PolyA/T_10+", "number"),
        ("PolyA/T_15+", "number"),
        ("PolyA/T_20+", "number"),
        ("Percent_10+", "number"),
        ("Percent_15+", "number"),
        ("Percent_20+", "number"),
    ]
    header_line = "\t".join(label for label, _type in summary_headers)
    summary_rows: list[list[str]] = []
    histogram_data: dict[str, dict[int, int]] = {}
    position_data: dict[str, list[tuple[int, int]]] = {}
    sample_order: list[str] = []

    with open(output_file, "w", encoding="utf-8") as out_handle:
        out_handle.write(header_line + "\n")
        for file_path in fastq_files:
            sample = sanitize_sample_name(file_path)
            (
                total,
                poly10,
                poly15,
                poly20,
                histogram,
                positions,
            ) = count_poly_runs(file_path)
            pct10 = format_percent(poly10, total)
            pct15 = format_percent(poly15, total)
            pct20 = format_percent(poly20, total)
            line = (
                f"{sample}\t{total}\t{poly10}\t{poly15}\t{poly20}\t"
                f"{pct10}\t{pct15}\t{pct20}"
            )
            out_handle.write(line + "\n")
            summary_rows.append(
                [
                    sample,
                    str(total),
                    str(poly10),
                    str(poly15),
                    str(poly20),
                    pct10,
                    pct15,
                    pct20,
                ]
            )
            histogram_data[sample] = histogram
            position_data[sample] = positions
            sample_order.append(sample)

    histogram_file = output_dir / "polyA_histogram.txt"
    histogram_series: dict[str, list[list[int]]] = {}
    with open(histogram_file, "w", encoding="utf-8") as hist_handle:
        hist_handle.write("Sample\tRun_Length\tRead_Count\n")
        for sample in sample_order:
            hist = histogram_data.get(sample, {})
            if not hist:
                histogram_series[sample] = []
                continue
            max_length = max(hist)
            series: list[list[int]] = []
            for length in range(HISTOGRAM_MIN_LENGTH, max_length + 1):
                count = hist.get(length, 0)
                hist_handle.write(f"{sample}\t{length}\t{count}\n")
                series.append([length, count])
            histogram_series[sample] = series

    position_headers = [
        ("Sample", "text"),
        ("Detected_Runs", "number"),
        ("Avg_nearest_end_offset", "number"),
        ("Median_nearest_end_offset", "number"),
    ]
    position_rows: list[list[str]] = []
    for sample in sample_order:
        offsets = position_data.get(sample, [])
        count = len(offsets)
        if count == 0:
            row = [sample, "0", "0.00", "0.00"]
        else:
            avg_offset = sum(offsets) / count
            median_offset = float(median(offsets))
            row = [
                sample,
                str(count),
                f"{avg_offset:.2f}",
                f"{median_offset:.2f}",
            ]
        position_rows.append(row)

    write_html_summary(
        output_dir,
        summary_headers,
        summary_rows,
        histogram_series,
        sample_order,
        position_headers,
        position_rows,
    )

    print(f"[polyat] Summary written to {output_file}")
    print(f"[polyat] Histogram written to {histogram_file}")
    print(f"[polyat] HTML report written to {output_dir / 'polyA_report.html'}")


if __name__ == "__main__":
    main()

