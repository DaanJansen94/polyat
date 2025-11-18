# polyat

polyat is a command-line tool that scans FASTQ files for poly-A/T homopolymers typically derived from carrier RNA (such as those in the QIAamp Viral RNA Mini Kit). Each read is evaluated once; the tool reports how many reads contain ≥10, ≥15, or ≥20 identical A/T bases together with their relative percentages. The run produces both a concise `polyA_counts.txt` QC table and an interactive `polyA_report.html` summary with per-column filters and a poly-A/T histogram.

## Features

- Accepts `.fastq`, `.fastq.gz`, `.fq`, and `.fq.gz` files
- Streams data without decompressing entire files to disk
- Counts each read at most once per length threshold (10/15/20 nt)
- Writes a tab-separated summary table suitable for downstream QC

## Installation

Clone and install from source code:

```bash
git clone https://github.com/DaanJansen94/polyat.git
cd polyat
pip install .
```

## Usage

General command:

```bash
polyat -i /path/to/fastq_dir -o /path/to/output_dir
```

Arguments:

```
-i / --input   Required input directory containing FASTQ/FQ files
-o / --output  Required output directory (created if missing)
```

`polyat` always writes `polyA_counts.txt` inside the output directory, counting
each read at most once per threshold (≥10/15/20 nt).

## Citation

If you use polyat in your research, please cite:

```
Jansen, D., Laumen, J., Siebenmann, E., & Vercauteren, K. (2025). polyat: Poly-A/T summarization (v0.1.1). [https://doi.org/10.5281/zenodo.17640044](https://doi.org/10.5281/zenodo.17640044).
```

## License

This project is licensed under the GNU General Public License v3.0 (GPL-3.0) - see the [LICENSE](LICENSE) file for details.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## Support

If you encounter any problems or have questions, please open an issue on GitHub.
