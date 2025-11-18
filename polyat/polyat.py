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
) -> tuple[int, int, int, int, dict[int, int], list[int]]:
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
    table_id: str,
    title: str,
    headers: list[tuple[str, str]],
    rows: list[list[str]],
    download_filename: str | None = None,
) -> str:
    parts: list[str] = [
        "<section>",
        f"<h2>{escape(title)}</h2>",
    ]
    if download_filename:
        parts.append(
            "<div class='download-controls'>"
            f"<button class='download-btn' data-download-table='{escape(table_id)}' "
            f"data-filename='{escape(download_filename)}'>Download table (TSV)</button>"
            f"<button class='download-btn secondary' data-download-table-pdf='{escape(table_id)}' "
            f"data-title='{escape(title)}'>Download table (PDF)</button>"
            "</div>"
        )
    parts.extend(
        [
        f"<table data-filterable='true' id='{escape(table_id)}'>",
        "<thead>",
        "<tr>",
        ]
    )
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
    combined_histogram: list[list[int]],
) -> None:
    output_file = output_dir / "polyA_report.html"
    css = (
        "body{font-family:Arial,sans-serif;margin:20px;background:#fefefe;}"
        "h1{margin-bottom:0.5em;}"
        "section{margin-bottom:30px;}"
        "#histogram-section label{display:block;margin-bottom:6px;font-weight:bold;}"
        "#histogram-section select{padding:6px 10px;margin-bottom:12px;}"
        "#histogram-canvas{width:100%;max-width:900px;height:auto;border:1px solid #ddd;background:#fff;}"
        "#histogram-section-all canvas{width:100%;max-width:900px;height:auto;border:1px solid #ddd;background:#fff;}"
        ".hist-tooltip{position:absolute;padding:6px 10px;background:rgba(0,0,0,0.75);"
        "color:#fff;border-radius:4px;font-size:12px;pointer-events:none;"
        "transform:translate(-50%,-120%);white-space:nowrap;display:none;z-index:10;}"
        ".download-controls{display:flex;gap:8px;flex-wrap:wrap;margin:8px 0;}"
        ".download-btn{display:inline-block;padding:6px 12px;border:1px solid #4a90e2;"
        "background:#4a90e2;color:#fff;border-radius:4px;cursor:pointer;font-size:0.9rem;}"
        ".download-btn.secondary{border-color:#555;background:#555;}"
        ".download-btn:hover{opacity:0.9;}"
        "table{border-collapse:collapse;width:100%;font-family:Arial,sans-serif;}"
        "th,td{border:1px solid #ccc;padding:8px;text-align:center;}"
        "th{background-color:#f4f4f4;}"
        ".filters th{background-color:#fafafa;}"
        ".filters input{width:100%;box-sizing:border-box;padding:4px;}"
        "tr:nth-child(even){background:#fafafa;}"
    )
    histogram_json = json.dumps(histogram_series)
    sample_json = json.dumps(sample_order)
    combined_json = json.dumps(combined_histogram)
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
        "function tableToTSV(table){"
        "return Array.from(table.rows).map(row=>Array.from(row.cells)"
        ".map(cell=>cell.innerText.replace(/\\s+/g,' ').trim()).join('\\t')).join('\\n');"
        "}"
        "function downloadBlob(content, filename, mime){"
        "let blob;"
        "if(content instanceof Blob){"
        "blob=content;"
        "}else if(content instanceof Uint8Array){"
        "blob=new Blob([content],{type:mime});"
        "}else{"
        "blob=new Blob([content],{type:mime});"
        "}"
        "const link=document.createElement('a');"
        "link.href=URL.createObjectURL(blob);"
        "link.download=filename;"
        "document.body.appendChild(link);"
        "link.click();"
        "setTimeout(()=>{URL.revokeObjectURL(link.href);link.remove();},0);"
        "}"
        "document.querySelectorAll('[data-download-table]').forEach(button=>{"
        "button.addEventListener('click',()=>{"
        "const table=document.getElementById(button.dataset.downloadTable);"
        "if(!table)return;"
        "const tsv=tableToTSV(table);"
        "const filename=button.dataset.filename||'table.tsv';"
        "downloadBlob(tsv, filename, 'text/tab-separated-values');"
        "});"
        "});"
        "document.querySelectorAll('[data-download-table-pdf]').forEach(button=>{"
        "button.addEventListener('click',()=>{"
        "const table=document.getElementById(button.dataset.downloadTablePdf);"
        "if(!table)return;"
        "const title=button.dataset.title || 'polyat_table';"
        "downloadTablePDF(table,title);"
        "});"
        "});"
        "function slugify(text){"
        "return (text||'polyat').toLowerCase().replace(/[^a-z0-9]+/g,'_').replace(/^_|_$/g,'') || 'polyat';"
        "}"
        "function getCanvasWithBackground(canvas){"
        "const exportCanvas=document.createElement('canvas');"
        "exportCanvas.width=canvas.width;"
        "exportCanvas.height=canvas.height;"
        "const ctx=exportCanvas.getContext('2d');"
        "ctx.fillStyle='#fff';"
        "ctx.fillRect(0,0,exportCanvas.width,exportCanvas.height);"
        "ctx.drawImage(canvas,0,0);"
        "return exportCanvas;"
        "}"
        "function renderTableCanvas(table){"
        "const rect=table.getBoundingClientRect();"
        "const scale=window.devicePixelRatio || 1;"
        "const width=Math.max(1,Math.round(rect.width));"
        "const height=Math.max(1,Math.round(rect.height));"
        "const canvas=document.createElement('canvas');"
        "canvas.width=width*scale;"
        "canvas.height=height*scale;"
        "const ctx=canvas.getContext('2d');"
        "ctx.scale(scale,scale);"
        "ctx.fillStyle='#fff';"
        "ctx.fillRect(0,0,width,height);"
        "const rows=Array.from(table.rows);"
        "let y=0;"
        "rows.forEach(row=>{"
        "const rowRect=row.getBoundingClientRect();"
        "const rowHeight=rowRect.height || 0;"
        "let x=0;"
        "Array.from(row.cells).forEach(cell=>{"
        "const cellRect=cell.getBoundingClientRect();"
        "const cellWidth=cellRect.width || 0;"
        "const styles=window.getComputedStyle(cell);"
        "const bg=styles.backgroundColor && styles.backgroundColor!=='rgba(0, 0, 0, 0)' ? styles.backgroundColor : '#fff';"
        "ctx.fillStyle=bg;"
        "ctx.fillRect(x,y,cellWidth,rowHeight);"
        "ctx.strokeStyle='#ccc';"
        "ctx.lineWidth=1;"
        "ctx.strokeRect(x+0.5,y+0.5,cellWidth,rowHeight);"
        "ctx.fillStyle=styles.color || '#000';"
        "const fontSize=parseFloat(styles.fontSize)||12;"
        "const fontFamily=styles.fontFamily || 'Arial';"
        "const fontWeight=styles.fontWeight && styles.fontWeight!=='normal' ? styles.fontWeight : '400';"
        "ctx.font=`${fontWeight} ${fontSize}px ${fontFamily}`;"
        "ctx.textBaseline='middle';"
        "const text=cell.innerText.trim();"
        "const padding=6;"
        "const textAlign=styles.textAlign || 'center';"
        "if(textAlign==='left'){"
        "ctx.textAlign='left';"
        "ctx.fillText(text,x+padding,y+rowHeight/2,Math.max(0,cellWidth-padding*2));"
        "}else if(textAlign==='right'){"
        "ctx.textAlign='right';"
        "ctx.fillText(text,x+cellWidth-padding,y+rowHeight/2,Math.max(0,cellWidth-padding*2));"
        "}else{"
        "ctx.textAlign='center';"
        "ctx.fillText(text,x+cellWidth/2,y+rowHeight/2,Math.max(0,cellWidth-padding*2));"
        "}"
        "x+=cellWidth;"
        "});"
        "y+=rowHeight;"
        "});"
        "return canvas;"
        "}"
        "function downloadTablePDF(table,title){"
        "const canvas=renderTableCanvas(table);"
        "const pdfBytes=canvasToPdfBytes(canvas);"
        "downloadBlob(pdfBytes, `${slugify(title)}.pdf`, 'application/pdf');"
        "}"
        f"const histogramData = {histogram_json};"
        f"const histogramSamples = {sample_json};"
        f"const combinedHistogramData = {combined_json};"
        "const sampleSelect = document.getElementById('histogram-sample');"
        "const sampleCanvas = document.getElementById('histogram-canvas');"
        "const combinedCanvas = document.getElementById('histogram-canvas-all');"
        "const tooltip = document.getElementById('histogram-tooltip');"
        "const sampleCtx = sampleCanvas ? sampleCanvas.getContext('2d') : null;"
        "const combinedCtx = combinedCanvas ? combinedCanvas.getContext('2d') : null;"
        "let sampleBars = [];"
        "let combinedBars = [];"
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
        "function drawHistogram(ctx, canvas, entries, messages){"
        "if(!ctx || !canvas)return [];"
        "ctx.clearRect(0,0,canvas.width,canvas.height);"
        "ctx.save();"
        "ctx.fillStyle='#fff';"
        "ctx.fillRect(0,0,canvas.width,canvas.height);"
        "ctx.restore();"
        "const bars=[];"
        "const margin=60;"
        "const width=canvas.width - margin*2;"
        "const height=canvas.height - margin*2;"
        "ctx.strokeStyle='#333';"
        "ctx.fillStyle='#333';"
        "ctx.lineWidth=1;"
        "ctx.font='12px Arial';"
        "ctx.beginPath();"
        "ctx.moveTo(margin, margin);"
        "ctx.lineTo(margin, margin + height);"
        "ctx.lineTo(margin + width, margin + height);"
        "ctx.stroke();"
        "if(!entries || entries.length===0){"
        "ctx.fillText(messages.empty, margin, margin);"
        "return bars;"
        "}"
        "const counts=entries.map(item=>item[1]);"
        "const maxCount=Math.max(...counts);"
        "if(!isFinite(maxCount) || maxCount<=0){"
        "ctx.fillText(messages.zero, margin, margin);"
        "return bars;"
        "}"
        "const tickCount=5;"
        "const tickStep=maxCount/tickCount;"
        "ctx.fillStyle='#000';"
        "ctx.textAlign='right';"
        "ctx.textBaseline='middle';"
        "for(let i=0;i<=tickCount;i++){"
        "const value=Math.round(i*tickStep);"
        "const y=margin + height - (value/maxCount)*height;"
        "ctx.beginPath();"
        "ctx.moveTo(margin-5,y);"
        "ctx.lineTo(margin,y);"
        "ctx.stroke();"
        "ctx.fillText(String(value), margin-8, y);"
        "ctx.strokeStyle='#eee';"
        "ctx.beginPath();"
        "ctx.moveTo(margin,y);"
        "ctx.lineTo(margin+width,y);"
        "ctx.stroke();"
        "ctx.strokeStyle='#333';"
        "}"
        "ctx.textAlign='center';"
        "ctx.textBaseline='middle';"
        "ctx.fillText('Count', margin-35, margin-20);"
        "ctx.fillText('Run length (nt)', margin + width/2, margin + height + 45);"
        "const minLength=entries[0][0];"
        "const maxLength=entries[entries.length-1][0];"
        "const lengthSpan=Math.max(1, maxLength - minLength);"
        "const desiredXTicks=Math.min(10, lengthSpan);"
        "const xStep=Math.max(1, Math.round(lengthSpan / desiredXTicks));"
        "ctx.textAlign='center';"
        "ctx.textBaseline='top';"
        "for(let value=minLength; value<=maxLength; value+=xStep){"
        "const ratio=(value - minLength)/lengthSpan;"
        "const x=margin + ratio*width;"
        "ctx.beginPath();"
        "ctx.moveTo(x, margin + height);"
        "ctx.lineTo(x, margin + height + 5);"
        "ctx.stroke();"
        "ctx.fillText(String(value), x, margin + height + 8);"
        "}"
        "const barWidth=width/entries.length;"
        "entries.forEach((item,index)=>{"
        "const length=item[0];"
        "const count=item[1];"
        "const barHeight=(count/maxCount)*height;"
        "const x=margin + index*barWidth + barWidth*0.1;"
        "const y=margin + height - barHeight;"
        "const w=barWidth*0.8;"
        "const h=barHeight || (count>0 ? 1 : 0);"
        "ctx.fillStyle='#4a90e2';"
        "ctx.fillRect(x,y,w,h);"
        "bars.push({x,y,width:w,height:h,count,length});"
        "});"
        "return bars;"
        "}"
        "function renderSampleHistogram(sample){"
        "const entries=histogramData[sample] || [];"
        "sampleBars = drawHistogram(sampleCtx, sampleCanvas, entries, {"
        f"empty:'No data for selected sample (no runs >={HISTOGRAM_MIN_LENGTH} nt).',"
        "zero:'All counts are zero for this sample.'"
        "});"
        "}"
        "function renderCombinedHistogram(){"
        "combinedBars = drawHistogram(combinedCtx, combinedCanvas, combinedHistogramData, {"
        "empty:'No combined histogram data available.',"
        "zero:'All combined counts are zero.'"
        "});"
        "}"
        "function attachCanvasHover(canvas, getBars){"
        "if(!canvas)return;"
        "canvas.addEventListener('mousemove', event=>{"
        "const rect=canvas.getBoundingClientRect();"
        "const x=event.clientX - rect.left;"
        "const y=event.clientY - rect.top;"
        "const bars=getBars?getBars():[];"
        "const bar=bars.find(b=>x>=b.x && x<=b.x+b.width && y>=b.y && y<=b.y+b.height);"
        "if(bar && bar.count>0){"
        "showTooltip(bar,event);"
        "}else{"
        "hideTooltip();"
        "}"
        "});"
        "canvas.addEventListener('mouseleave', hideTooltip);"
        "}"
        "function downloadCanvasImage(targetCanvas, filename){"
        "if(!targetCanvas)return;"
        "const exportCanvas=getCanvasWithBackground(targetCanvas);"
        "const link=document.createElement('a');"
        "link.href=exportCanvas.toDataURL('image/png');"
        "link.download=filename;"
        "document.body.appendChild(link);"
        "link.click();"
        "link.remove();"
        "}"
        "function canvasToPdfBytes(targetCanvas){"
        "if(!targetCanvas)return new Uint8Array();"
        "const exportCanvas=getCanvasWithBackground(targetCanvas);"
        "const jpegData=exportCanvas.toDataURL('image/jpeg',0.95);"
        "const base64=jpegData.split(',')[1];"
        "const binary=atob(base64);"
        "const imageBytes=new Uint8Array(binary.length);"
        "for(let i=0;i<binary.length;i++){"
        "imageBytes[i]=binary.charCodeAt(i);"
        "}"
        "const width=targetCanvas.width;"
        "const height=targetCanvas.height;"
        "const encoder=new TextEncoder();"
        "const chunks=[];"
        "const offsets=[0];"
        "let position=0;"
        "function append(data){"
        "let bytes;"
        "if(typeof data==='string'){"
        "bytes=encoder.encode(data);"
        "}else{"
        "bytes=data;"
        "}"
        "chunks.push(bytes);"
        "position+=bytes.length;"
        "}"
        "append('%PDF-1.4\\n');"
        "offsets.push(position);append('1 0 obj << /Type /Catalog /Pages 2 0 R >> endobj\\n');"
        "offsets.push(position);append('2 0 obj << /Type /Pages /Kids [3 0 R] /Count 1 >> endobj\\n');"
        "offsets.push(position);append(`3 0 obj << /Type /Page /Parent 2 0 R /MediaBox [0 0 ${width} ${height}] /Resources << /ProcSet [/PDF /ImageC] /XObject << /Im0 4 0 R >> >> /Contents 5 0 R >> endobj\\n`);"
        "offsets.push(position);append(`4 0 obj << /Type /XObject /Subtype /Image /Width ${width} /Height ${height} /ColorSpace /DeviceRGB /BitsPerComponent 8 /Filter /DCTDecode /Length ${imageBytes.length} >> stream\\n`);"
        "append(imageBytes);"
        "append('\\nendstream\\nendobj\\n');"
        "const content=`q\\n${width} 0 0 ${height} 0 0 cm\\n/Im0 Do\\nQ\\n`;"
        "offsets.push(position);append(`5 0 obj << /Length ${content.length} >> stream\\n${content}endstream\\nendobj\\n`);"
        "const xrefStart=position;"
        "append('xref\\n0 6\\n0000000000 65535 f \\n');"
        "for(let i=1;i<=5;i++){"
        "append(`${offsets[i].toString().padStart(10,'0')} 00000 n \\n`);"
        "}"
        "append(`trailer << /Size 6 /Root 1 0 R >>\\nstartxref\\n${xrefStart}\\n%%EOF`);"
        "let total=0;"
        "chunks.forEach(chunk=>total+=chunk.length);"
        "const pdfBytes=new Uint8Array(total);"
        "let offset=0;"
        "chunks.forEach(chunk=>{"
        "pdfBytes.set(chunk,offset);"
        "offset+=chunk.length;"
        "});"
        "return pdfBytes;"
        "}"
        "function downloadCanvasPDF(targetCanvas, filename){"
        "if(!targetCanvas)return;"
        "const pdfBytes=canvasToPdfBytes(targetCanvas);"
        "downloadBlob(pdfBytes, filename, 'application/pdf');"
        "}"
        "attachCanvasHover(sampleCanvas, ()=>sampleBars);"
        "attachCanvasHover(combinedCanvas, ()=>combinedBars);"
        "const sampleDownloadBtn=document.getElementById('download-sample-histogram');"
        "if(sampleDownloadBtn && sampleCanvas){"
        "sampleDownloadBtn.addEventListener('click',()=>{"
        "const suffix=sampleSelect && sampleSelect.value ? sampleSelect.value : 'sample';"
        "downloadCanvasImage(sampleCanvas, `polyA_histogram_${suffix}.png`);"
        "});"
        "}"
        "const sampleDownloadPdfBtn=document.getElementById('download-sample-histogram-pdf');"
        "if(sampleDownloadPdfBtn && sampleCanvas){"
        "sampleDownloadPdfBtn.addEventListener('click',()=>{"
        "const suffix=sampleSelect && sampleSelect.value ? sampleSelect.value : 'sample';"
        "downloadCanvasPDF(sampleCanvas, `polyA_histogram_${suffix}.pdf`);"
        "});"
        "}"
        "const combinedDownloadBtn=document.getElementById('download-combined-histogram');"
        "if(combinedDownloadBtn && combinedCanvas){"
        "combinedDownloadBtn.addEventListener('click',()=>{"
        "downloadCanvasImage(combinedCanvas, 'polyA_histogram_combined.png');"
        "});"
        "}"
        "const combinedDownloadPdfBtn=document.getElementById('download-combined-histogram-pdf');"
        "if(combinedDownloadPdfBtn && combinedCanvas){"
        "combinedDownloadPdfBtn.addEventListener('click',()=>{"
        "downloadCanvasPDF(combinedCanvas, 'polyA_histogram_combined.pdf');"
        "});"
        "}"
        "function initHistogram(){"
        "if(sampleSelect){"
        "sampleSelect.innerHTML='';"
        "histogramSamples.forEach(sample=>{"
        "const option=document.createElement('option');"
        "option.value=sample;"
        "option.textContent=sample;"
        "sampleSelect.appendChild(option);"
        "});"
        "if(histogramSamples.length>0){"
        "sampleSelect.value=histogramSamples[0];"
        "renderSampleHistogram(sampleSelect.value);"
        "}else{"
        "renderSampleHistogram('');"
        "}"
        "sampleSelect.addEventListener('change',()=>{"
        "hideTooltip();"
        "renderSampleHistogram(sampleSelect.value);"
        "});"
        "}"
        "renderCombinedHistogram();"
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
            "polyat-summary-table",
            "polyA/T Summary",
            summary_headers,
            summary_rows,
            "polyA_counts_table.tsv",
        ),
        "<section id='histogram-section'>",
        f"<h2>polyA/T Histogram (>={HISTOGRAM_MIN_LENGTH} nt)</h2>",
        "<label for='histogram-sample'>Sample</label>",
        "<select id='histogram-sample'></select>",
        "<div class='download-controls'>"
        "<button id='download-sample-histogram' class='download-btn'>Download histogram (PNG)</button>"
        "<button id='download-sample-histogram-pdf' class='download-btn secondary'>Download histogram (PDF)</button>"
        "</div>",
        "<canvas id='histogram-canvas' width='900' height='400'></canvas>",
        "<div id='histogram-tooltip' class='hist-tooltip'></div>",
        "</section>",
        "<section id='histogram-section-all'>",
        "<h2>Combined polyA/T Histogram</h2>",
        "<div class='download-controls'>"
        "<button id='download-combined-histogram' class='download-btn'>Download combined histogram (PNG)</button>"
        "<button id='download-combined-histogram-pdf' class='download-btn secondary'>Download combined histogram (PDF)</button>"
        "</div>",
        "<canvas id='histogram-canvas-all' width='900' height='400'></canvas>",
        "</section>",
        build_filterable_table(
            "polyat-position-table",
            "polyA/T Position Offsets",
            position_headers,
            position_rows,
            "polyA_offsets_table.tsv",
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
    position_data: dict[str, list[int]] = {}
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
    combined_counts: dict[int, int] = defaultdict(int)
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
                combined_counts[length] += count
            histogram_series[sample] = series

    combined_histogram: list[list[int]] = []
    if combined_counts:
        max_length = max(combined_counts)
        for length in range(HISTOGRAM_MIN_LENGTH, max_length + 1):
            combined_histogram.append([length, combined_counts.get(length, 0)])

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
        combined_histogram,
    )

    print(f"[polyat] Summary written to {output_file}")
    print(f"[polyat] Histogram written to {histogram_file}")
    print(f"[polyat] HTML report written to {output_dir / 'polyA_report.html'}")


if __name__ == "__main__":
    main()

