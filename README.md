# nf-read-ref-mapping
Maps Illumina and/or ONT reads to a per-sample reference with minimap2, then summarises how the
mapped reads are distributed across each reference.

## Samplesheet

CSV or TSV with columns `sample,platform,fastq_1,fastq_2,reference`.

- `platform`: `illumina` (minimap2 `-x sr`) or `ont` (`-x map-ont`)
- `fastq_2`: leave empty for single-end reads
- `reference`: FASTA (gzipped or not); each reference is indexed once per preset and shared

Each `sample` + `platform` combination must be unique.

## Run

```bash
nextflow run . -profile conda --input samplesheet.csv --outdir mapping_results
```

| Parameter | Default | |
|---|---|---|
| `--window_size` | `1000000` | Bin size for the distribution along each reference |
| `--min_mapq` | `20` | Threshold reported as "MAPQ >= x" |
| `--top_contigs` | `25` | Contigs listed per sample in the HTML report |
| `--preset_illumina` / `--preset_ont` | `sr` / `map-ont` | minimap2 presets (e.g. `lr:hq` for R10 SUP) |
| `--save_bam` | `true` | Publish sorted, indexed BAMs |

Extra minimap2 options go in `ext.args` for `MINIMAP2_ALIGN` in `conf/modules.config`.

## Outputs

- `mapping_summary.tsv`: reads, mapping rate, MAPQ >= threshold, properly paired %, contigs hit, top contig, breadth and depth
- `contig_distribution.tsv`: primary reads, share of mapped reads, enrichment, breadth, depth and MAPQ per contig
- `window_distribution.tsv`: primary read starts per window, zeros included
- `mapq_distribution.tsv`: MAPQ histogram of primary mapped reads
- `mapping_summary.html`: overview table, plus per-sample distribution plot, top contigs and MAPQ histogram
- `alignments/<platform>/`, `stats/<platform>/`, `pipeline_info/software_versions.yml`

Enrichment = a contig's share of mapped reads / its share of reference length. Values well above 1 mark
contigs attracting more reads than their size predicts.
