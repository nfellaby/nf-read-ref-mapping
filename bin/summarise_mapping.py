#!/usr/bin/env python3
"""Summarise read mapping across samples: mapping rates, MAPQ, and how mapped reads are
distributed across each reference (per contig and per window).

Inputs
------
--manifest   Header-less TSV: sample, platform, reference, layout, prefix
--stats-dir  Directory with <prefix>.flagstat.json, .coverage.tsv, .mapq.tsv and .windows.tsv

Outputs (in the working directory)
-------
mapping_summary.tsv       one row per sample/platform
contig_distribution.tsv   primary mapped reads, depth and breadth per contig
window_distribution.tsv   primary read starts per window (zeros included)
mapq_distribution.tsv     MAPQ histogram of primary mapped reads
mapping_summary.html      self-contained report
"""

import argparse
import csv
import datetime
import html
import json
import math
import sys
from collections import defaultdict
from pathlib import Path

MAX_PLOT_BARS = 1000
MAPQ_BINS = [(0, 0), (1, 9), (10, 19), (20, 29), (30, 39), (40, 49), (50, 59), (60, 255)]

SUMMARY_COLUMNS = [
    "sample", "platform", "reference", "layout",
    "reads", "mapped", "unmapped", "mapped_pct",
    "mapq_ge_min", "mapq_ge_min_pct", "properly_paired_pct",
    "secondary", "supplementary",
    "contigs_with_reads", "contigs_total",
    "top_contig", "top_contig_pct_of_mapped",
    "reference_length", "covered_bases", "breadth_pct", "mean_depth",
]
CONTIG_COLUMNS = [
    "sample", "platform", "reference", "contig", "length",
    "primary_reads", "pct_of_mapped", "pct_of_reference_length", "enrichment",
    "covered_bases", "breadth_pct", "mean_depth", "mean_mapq",
]
WINDOW_COLUMNS = ["sample", "platform", "reference", "contig", "start", "end", "reads"]
MAPQ_COLUMNS = ["sample", "platform", "reference", "mapq", "reads"]


def pct(part, whole):
    return 100.0 * part / whole if whole else 0.0


def read_manifest(path):
    rows = []
    with open(path) as handle:
        for line in handle:
            line = line.rstrip("\r\n")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) != 5:
                sys.exit(f"Malformed manifest line (expected 5 fields): {line!r}")
            rows.append(dict(zip(["sample", "platform", "reference", "layout", "prefix"], fields)))
    if not rows:
        sys.exit("Manifest is empty: no mapping statistics to report")
    return sorted(rows, key=lambda r: (r["platform"], r["sample"], r["reference"]))


def read_coverage(path):
    contigs = []
    with open(path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            contigs.append(
                {
                    "contig": row["#rname"],
                    "length": int(row["endpos"]),
                    "covered_bases": int(row["covbases"]),
                    "mean_depth": float(row["meandepth"]),
                    "mean_mapq": float(row["meanmapq"]),
                }
            )
    return contigs


def read_counts(path, key_fields):
    counts = defaultdict(int)
    with open(path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) != key_fields + 1:
                continue
            key = fields[0] if key_fields == 1 else tuple(fields[:key_fields])
            counts[key] += int(fields[-1])
    return counts


def load_sample(entry, stats_dir, window_size, min_mapq):
    base = Path(stats_dir) / entry["prefix"]
    flagstat = json.loads(Path(f"{base}.flagstat.json").read_text())["QC-passed reads"]
    contigs = read_coverage(f"{base}.coverage.tsv")
    mapq = {int(q): n for q, n in read_counts(f"{base}.mapq.tsv", 1).items()}
    window_counts = read_counts(f"{base}.windows.tsv", 2)

    reads = flagstat["primary"]
    mapped = flagstat["primary mapped"]
    reference_length = sum(c["length"] for c in contigs)

    # Primary reads per contig from the window counts
    per_contig = defaultdict(int)
    for (contig, _bin), n in window_counts.items():
        per_contig[contig] += n

    contig_rows = []
    for c in contigs:
        n = per_contig.get(c["contig"], 0)
        share_mapped = pct(n, mapped)
        share_length = pct(c["length"], reference_length)
        contig_rows.append(
            {
                **c,
                "primary_reads": n,
                "pct_of_mapped": share_mapped,
                "pct_of_reference_length": share_length,
                "enrichment": share_mapped / share_length if share_length else 0.0,
                "breadth_pct": pct(c["covered_bases"], c["length"]),
            }
        )

    windows = []
    for c in contigs:
        n_windows = max(1, math.ceil(c["length"] / window_size))
        for b in range(n_windows):
            windows.append(
                {
                    "contig": c["contig"],
                    "start": b * window_size + 1,
                    "end": min((b + 1) * window_size, c["length"]),
                    "reads": window_counts.get((c["contig"], str(b)), 0),
                }
            )

    mapq_ge = sum(n for q, n in mapq.items() if q >= min_mapq)
    covered = sum(c["covered_bases"] for c in contigs)
    top = max(contig_rows, key=lambda r: r["primary_reads"], default=None)
    paired = entry["layout"] == "paired"

    summary = {
        "sample": entry["sample"],
        "platform": entry["platform"],
        "reference": entry["reference"],
        "layout": entry["layout"],
        "reads": reads,
        "mapped": mapped,
        "unmapped": reads - mapped,
        "mapped_pct": pct(mapped, reads),
        "mapq_ge_min": mapq_ge,
        "mapq_ge_min_pct": pct(mapq_ge, reads),
        "properly_paired_pct": pct(flagstat["properly paired"], reads) if paired else "NA",
        "secondary": flagstat["secondary"],
        "supplementary": flagstat["supplementary"],
        "contigs_with_reads": sum(1 for r in contig_rows if r["primary_reads"] > 0),
        "contigs_total": len(contig_rows),
        "top_contig": top["contig"] if top and top["primary_reads"] else "NA",
        "top_contig_pct_of_mapped": top["pct_of_mapped"] if top and top["primary_reads"] else 0.0,
        "reference_length": reference_length,
        "covered_bases": covered,
        "breadth_pct": pct(covered, reference_length),
        "mean_depth": (
            sum(c["mean_depth"] * c["length"] for c in contigs) / reference_length if reference_length else 0.0
        ),
    }
    return {"summary": summary, "contigs": contig_rows, "windows": windows, "mapq": mapq}


def write_tsv(path, columns, rows):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: (f"{v:.4f}" if isinstance(v, float) else v) for k, v in row.items()})


# ---------------------------------------------------------------- HTML helpers

esc = html.escape


def fmt_int(value):
    return f"{value:,}"


def fmt_pct(value):
    return value if isinstance(value, str) else f"{value:.2f}%"


def bar(value):
    width = max(0.0, min(100.0, value))
    return f'<div class="bar"><span style="width:{width:.2f}%"></span></div>'


def window_plot(windows, window_size):
    """Read starts along the reference, contigs in reference order, alternating colours."""
    if not windows:
        return '<p class="muted">No windows.</p>'
    group = max(1, math.ceil(len(windows) / MAX_PLOT_BARS))
    bars = []
    for i in range(0, len(windows), group):
        chunk = windows[i:i + group]
        bars.append((chunk[0]["contig"], sum(w["reads"] for w in chunk), chunk[0], chunk[-1]))

    contig_order = {}
    for contig, *_ in bars:
        contig_order.setdefault(contig, len(contig_order))

    peak = max(b[1] for b in bars) or 1
    width, height = 1000, 140
    bar_w = width / len(bars)
    rects = []
    for i, (contig, n, first, last) in enumerate(bars):
        h = (n / peak) * (height - 4)
        cls = "a" if contig_order[contig] % 2 == 0 else "b"
        title = f"{contig}:{first['start']:,}-{last['end']:,} · {n:,} reads"
        rects.append(
            f'<rect class="{cls}" x="{i * bar_w:.2f}" y="{height - h:.2f}" '
            f'width="{max(bar_w, 0.6):.2f}" height="{h:.2f}"><title>{esc(title)}</title></rect>'
        )
    span = window_size * group
    span_label = f"{span / 1e6:g} Mb" if span >= 1e6 else f"{span / 1e3:g} kb"
    return (
        f'<svg class="plot" viewBox="0 0 {width} {height}" preserveAspectRatio="none" role="img" '
        f'aria-label="Read starts along the reference">{"".join(rects)}</svg>'
        f'<p class="muted small">Primary read starts per {span_label} bin, contigs in reference order '
        f"(alternating shades). Peak bin: {peak:,} reads. Hover a bar for its location.</p>"
    )


def mapq_plot(mapq):
    counts = [sum(n for q, n in mapq.items() if lo <= q <= hi) for lo, hi in MAPQ_BINS]
    labels = ["0" if lo == hi else (f"{lo}+" if hi == 255 else f"{lo}-{hi}") for lo, hi in MAPQ_BINS]
    peak = max(counts) or 1
    width, height, pad = 400, 130, 22
    bar_w = width / len(counts)
    parts = []
    for i, (n, label) in enumerate(zip(counts, labels)):
        h = (n / peak) * (height - pad - 4)
        x = i * bar_w
        parts.append(
            f'<rect class="a" x="{x + 4:.1f}" y="{height - pad - h:.1f}" width="{bar_w - 8:.1f}" '
            f'height="{h:.1f}"><title>MAPQ {label}: {n:,} reads</title></rect>'
            f'<text x="{x + bar_w / 2:.1f}" y="{height - 6}" text-anchor="middle">{label}</text>'
        )
    return (
        f'<svg class="mapq" viewBox="0 0 {width} {height}" role="img" aria-label="MAPQ distribution">'
        f'{"".join(parts)}</svg>'
    )


def contig_table(contigs, top_n):
    ranked = sorted(contigs, key=lambda r: r["primary_reads"], reverse=True)
    shown = [r for r in ranked if r["primary_reads"] > 0][:top_n]
    if not shown:
        return '<p class="muted">No mapped reads.</p>'
    rows = "".join(
        f"<tr><td class='mono'>{esc(r['contig'])}</td><td class='num'>{fmt_int(r['length'])}</td>"
        f"<td class='num'>{fmt_int(r['primary_reads'])}</td><td class='num'>{fmt_pct(r['pct_of_mapped'])}</td>"
        f"<td class='barcell'>{bar(r['pct_of_mapped'])}</td>"
        f"<td class='num'>{r['enrichment']:.2f}</td><td class='num'>{fmt_pct(r['breadth_pct'])}</td>"
        f"<td class='num'>{r['mean_depth']:.2f}</td><td class='num'>{r['mean_mapq']:.1f}</td></tr>"
        for r in shown
    )
    remaining = sum(1 for r in ranked if r["primary_reads"] > 0) - len(shown)
    more = f'<p class="muted small">{remaining} more contigs with reads in contig_distribution.tsv.</p>' if remaining > 0 else ""
    return (
        '<div class="scroll"><table><thead><tr><th>Contig</th><th class="num">Length</th>'
        '<th class="num">Reads</th><th class="num">% of mapped</th><th></th>'
        '<th class="num">Enrichment</th><th class="num">Breadth</th><th class="num">Mean depth</th>'
        f'<th class="num">Mean MAPQ</th></tr></thead><tbody>{rows}</tbody></table></div>{more}'
    )


def render_html(samples, window_size, min_mapq, top_n):
    generated = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")

    summary_rows = "".join(
        f"<tr><td><a href='#s{i}'>{esc(s['summary']['sample'])}</a></td>"
        f"<td>{esc(s['summary']['platform'])}</td><td class='mono'>{esc(s['summary']['reference'])}</td>"
        f"<td class='num'>{fmt_int(s['summary']['reads'])}</td>"
        f"<td class='num'>{fmt_pct(s['summary']['mapped_pct'])}</td><td class='barcell'>{bar(s['summary']['mapped_pct'])}</td>"
        f"<td class='num'>{fmt_pct(s['summary']['mapq_ge_min_pct'])}</td>"
        f"<td class='num'>{fmt_pct(s['summary']['properly_paired_pct'])}</td>"
        f"<td class='num'>{s['summary']['contigs_with_reads']:,} / {s['summary']['contigs_total']:,}</td>"
        f"<td class='num'>{fmt_pct(s['summary']['breadth_pct'])}</td></tr>"
        for i, s in enumerate(samples)
    )

    sections = "".join(
        f"""
  <details id="s{i}"{' open' if len(samples) <= 3 else ''}>
    <summary><b>{esc(s['summary']['sample'])}</b> · {esc(s['summary']['platform'])} ·
      <span class="mono">{esc(s['summary']['reference'])}</span> ·
      {fmt_pct(s['summary']['mapped_pct'])} mapped</summary>
    <h3>Distribution along the reference</h3>
    {window_plot(s['windows'], window_size)}
    <div class="split">
      <div><h3>Top contigs</h3>{contig_table(s['contigs'], top_n)}</div>
      <div class="narrow"><h3>MAPQ (primary mapped reads)</h3>{mapq_plot(s['mapq'])}</div>
    </div>
  </details>"""
        for i, s in enumerate(samples)
    )

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Read mapping summary</title>
<style>
  :root {{
    --bg: #ffffff; --fg: #1d1d1f; --muted: #6b6b70; --line: #e3e3e6; --card: #f6f6f8;
    --a: #2f6fbf; --b: #8fb3e3; --track: #e6ebf2;
    color-scheme: light;
  }}
  @media (prefers-color-scheme: dark) {{
    :root {{
      --bg: #17181a; --fg: #ececee; --muted: #9a9aa1; --line: #2c2d31; --card: #202125;
      --a: #6aa3ea; --b: #3b5f8f; --track: #2b3240;
      color-scheme: dark;
    }}
  }}
  body {{ margin: 0; background: var(--bg); color: var(--fg);
         font: 15px/1.5 system-ui, -apple-system, "Segoe UI", Roboto, sans-serif; }}
  main {{ max-width: 1200px; margin: 0 auto; padding: 2rem 1.25rem 3rem; }}
  h1 {{ font-size: 1.6rem; margin: 0 0 .25rem; }}
  h2 {{ font-size: 1.15rem; margin: 2rem 0 .75rem; }}
  h3 {{ font-size: .95rem; margin: 1.2rem 0 .5rem; color: var(--muted); }}
  a {{ color: var(--a); }}
  .muted {{ color: var(--muted); }}
  .small {{ font-size: .85rem; margin-top: .3rem; }}
  .mono {{ font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace; font-size: .85rem; }}
  .scroll {{ overflow-x: auto; }}
  table {{ border-collapse: collapse; width: 100%; font-size: .9rem; }}
  th, td {{ padding: .4rem .55rem; border-bottom: 1px solid var(--line); text-align: left; white-space: nowrap; }}
  th {{ font-weight: 600; color: var(--muted); }}
  .num {{ text-align: right; font-variant-numeric: tabular-nums; }}
  .barcell {{ width: 16%; min-width: 90px; }}
  .bar {{ background: var(--track); border-radius: 4px; height: 9px; overflow: hidden; }}
  .bar span {{ display: block; height: 100%; background: var(--a); }}
  details {{ background: var(--card); border-radius: 10px; padding: .8rem 1.1rem; margin: .8rem 0; }}
  summary {{ cursor: pointer; }}
  svg.plot {{ width: 100%; height: 140px; display: block; background: var(--bg); border-radius: 6px; }}
  svg.mapq {{ width: 100%; max-width: 400px; display: block; }}
  svg rect.a {{ fill: var(--a); }}
  svg rect.b {{ fill: var(--b); }}
  svg text {{ fill: var(--muted); font-size: 11px; }}
  .split {{ display: grid; grid-template-columns: minmax(0, 2fr) minmax(0, 1fr); gap: 1.5rem; }}
  @media (max-width: 800px) {{ .split {{ grid-template-columns: 1fr; }} }}
  .note {{ font-size: .88rem; color: var(--muted); max-width: 80ch; }}
</style>
</head>
<body>
<main>
  <h1>Read mapping summary</h1>
  <p class="muted">Generated {generated} · {len(samples)} sample/platform combinations ·
    window {window_size:,} bp · MAPQ threshold {min_mapq}</p>

  <h2>Overview</h2>
  <div class="scroll">
    <table>
      <thead><tr>
        <th>Sample</th><th>Platform</th><th>Reference</th><th class="num">Reads</th>
        <th class="num">Mapped</th><th></th><th class="num">MAPQ &ge; {min_mapq}</th>
        <th class="num">Properly paired</th><th class="num">Contigs hit</th><th class="num">Breadth</th>
      </tr></thead>
      <tbody>{summary_rows}</tbody>
    </table>
  </div>

  <h2>Per sample</h2>
  {sections}

  <p class="note">Reads are primary alignments; paired-end counts include both mates. Enrichment is a contig's
  share of mapped reads divided by its share of the reference length, so values well above 1 mark contigs that
  attract more reads than their size predicts. Breadth and depth come from samtools coverage.</p>
</main>
</body>
</html>
"""


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--stats-dir", default=".")
    parser.add_argument("--window-size", type=int, default=1_000_000)
    parser.add_argument("--min-mapq", type=int, default=20)
    parser.add_argument("--top-contigs", type=int, default=25)
    args = parser.parse_args()

    samples = [
        load_sample(entry, args.stats_dir, args.window_size, args.min_mapq)
        for entry in read_manifest(args.manifest)
    ]

    def tag(sample, rows):
        s = sample["summary"]
        return [{"sample": s["sample"], "platform": s["platform"], "reference": s["reference"], **r} for r in rows]

    write_tsv("mapping_summary.tsv", SUMMARY_COLUMNS, [s["summary"] for s in samples])
    write_tsv("contig_distribution.tsv", CONTIG_COLUMNS, [r for s in samples for r in tag(s, s["contigs"])])
    write_tsv("window_distribution.tsv", WINDOW_COLUMNS, [r for s in samples for r in tag(s, s["windows"])])
    write_tsv(
        "mapq_distribution.tsv",
        MAPQ_COLUMNS,
        [r for s in samples for r in tag(s, [{"mapq": q, "reads": n} for q, n in sorted(s["mapq"].items())])],
    )
    Path("mapping_summary.html").write_text(render_html(samples, args.window_size, args.min_mapq, args.top_contigs))


if __name__ == "__main__":
    main()
