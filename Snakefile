from __future__ import annotations

import csv
import contextlib
import re
from pathlib import Path
from typing import Iterable
from urllib.parse import unquote

from openpyxl import Workbook
from openpyxl.utils import get_column_letter

try:
    import gffutils
except ImportError:
    gffutils = None


configfile: "config.yaml"


# ---------------------------
# Config and validation
# ---------------------------

def normalize_path(raw_path: str) -> Path:
    return Path(raw_path).expanduser().resolve()


def parse_bool(value, key_name: str, default: bool) -> bool:
    if value is None:
        return default
    if not isinstance(value, bool):
        raise ValueError(f"'{key_name}' must be true/false")
    return value


def parse_positive_int(value, key_name: str, default: int) -> int:
    parsed = default if value is None else value
    if not isinstance(parsed, int) or parsed <= 0:
        raise ValueError(f"'{key_name}' must be a positive integer")
    return parsed


MAX_THREADS = parse_positive_int(config.get("threads"), "threads", 8)
PREPARE_READS_THREADS = parse_positive_int(config.get("prepare_reads_threads"), "prepare_reads_threads", 4)
MAP_THREADS = parse_positive_int(config.get("map_threads"), "map_threads", MAX_THREADS)

SRA_IDS_FILE_RAW = config.get("sra_ids_file")
if not isinstance(SRA_IDS_FILE_RAW, str) or not SRA_IDS_FILE_RAW.strip():
    raise ValueError("Missing required 'sra_ids_file' in config")
SRA_IDS_FILE = normalize_path(SRA_IDS_FILE_RAW)
if not SRA_IDS_FILE.is_file():
    raise FileNotFoundError(f"SRA IDs file not found: {SRA_IDS_FILE}")

GENOMES_RAW = config.get("genomes")
if not isinstance(GENOMES_RAW, list) or not GENOMES_RAW:
    raise ValueError("'genomes' must be a non-empty list")

INCLUDE_ZERO_ROWS = parse_bool(config.get("include_zero_rows"), "include_zero_rows", False)
KEEP_AUX = parse_bool(config.get("keep_aux"), "keep_aux", True)
KEEP_MAPPING = parse_bool(config.get("keep_mapping"), "keep_mapping", True)

GENOME_BY_ID: dict[str, dict[str, Path]] = {}
for idx, genome in enumerate(GENOMES_RAW, start=1):
    if not isinstance(genome, dict):
        raise ValueError(f"genomes[{idx}] must be an object")

    genome_id = genome.get("genome_id")
    fasta = genome.get("fasta")
    gff3 = genome.get("gff3")

    if not isinstance(genome_id, str) or not genome_id.strip():
        raise ValueError(f"genomes[{idx}].genome_id is required")
    if not isinstance(fasta, str) or not fasta.strip():
        raise ValueError(f"genomes[{idx}].fasta is required")
    if not isinstance(gff3, str) or not gff3.strip():
        raise ValueError(f"genomes[{idx}].gff3 is required")

    if genome_id in GENOME_BY_ID:
        raise ValueError(f"Duplicated genome_id in config: {genome_id}")

    fasta_path = normalize_path(fasta)
    gff3_path = normalize_path(gff3)

    if not fasta_path.is_file():
        raise FileNotFoundError(f"FASTA file not found for {genome_id}: {fasta_path}")
    if not gff3_path.is_file():
        raise FileNotFoundError(f"GFF3 file not found for {genome_id}: {gff3_path}")

    GENOME_BY_ID[genome_id] = {"fasta": fasta_path, "gff3": gff3_path}


def read_sra_ids(path: Path) -> list[str]:
    sra_ids: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        entry = line.strip()
        if entry:
            sra_ids.append(entry)
    if not sra_ids:
        raise ValueError(f"No SRA IDs found in file: {path}")
    return sra_ids


SRA_IDS = read_sra_ids(SRA_IDS_FILE)
GENOME_IDS = sorted(GENOME_BY_ID.keys())


def maybe_temp(path: str, keep: bool) -> str:
    return path if keep else temp(path)


# ---------------------------
# Python helpers (only where shell is not enough)
# ---------------------------

GFF_FEATURES_CACHE: dict[Path, list[tuple[str, int, int, str, str, str]]] = {}


def _first_attr(attrs: dict[str, list[str]], keys: tuple[str, ...]) -> str:
    for key in keys:
        values = attrs.get(key)
        if values:
            return unquote(values[0])
    return ""


def _feature_type_rank(feature_type: str) -> int:
    ft = feature_type.lower()
    if ft == "gene":
        return 0
    if ft == "cds":
        return 1
    if ft in {"rrna", "trna", "ncrna", "tmrna", "snrna", "snorna"}:
        return 2
    return 3


def _select_representative_features(
    rows: list[tuple[str, int, int, str, str, str, str]],
) -> list[tuple[str, int, int, str, str, str]]:
    # Keep a single representative per exact locus to avoid gene/CDS/RNA/exon duplicates.
    best_by_locus: dict[tuple[str, int, int, str], tuple[int, int, tuple[str, int, int, str, str, str, str]]] = {}
    for idx, row in enumerate(rows):
        seqid, start0, end1, feature_id, feature_name, strand, feature_type = row
        locus = (seqid, start0, end1, strand)
        candidate = (_feature_type_rank(feature_type), idx, row)
        current = best_by_locus.get(locus)
        if current is None or candidate < current:
            best_by_locus[locus] = candidate

    selected = [entry[2] for entry in sorted(best_by_locus.values(), key=lambda x: x[1])]
    return [(seqid, start0, end1, feature_id, feature_name, strand) for seqid, start0, end1, feature_id, feature_name, strand, _ in selected]


def _iter_gff_features_fallback(gff_path: Path) -> Iterable[tuple[str, int, int, str, str, str, str]]:
    with gff_path.open("r", encoding="utf-8") as handle:
        for lineno, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue

            seqid, _, feature_type, start, end, _, strand, _, attrs = fields[:9]
            if not start.isdigit() or not end.isdigit():
                continue

            start_i = int(start)
            end_i = int(end)
            if start_i > end_i:
                continue

            feature_id = ""
            feature_name = ""

            for attr in attrs.split(";"):
                if "=" not in attr:
                    continue
                key, value = attr.split("=", 1)
                if key == "ID" and not feature_id:
                    feature_id = unquote(value)
                if key in {"Name", "NAME"} and not feature_name:
                    feature_name = unquote(value)

            if not feature_id:
                feature_id = f"feature_{lineno}"
            if not feature_name:
                feature_name = feature_id

            clean_strand = strand if strand in {"+", "-"} else "."
            yield seqid, start_i - 1, end_i, feature_id, feature_name, clean_strand, feature_type


def _iter_gff_features_gffutils(gff_path: Path) -> Iterable[tuple[str, int, int, str, str, str, str]]:
    db = gffutils.create_db(
        str(gff_path),
        dbfn=":memory:",
        force=True,
        keep_order=True,
        merge_strategy="create_unique",
        sort_attribute_values=True,
        disable_infer_genes=True,
        disable_infer_transcripts=True,
    )

    for lineno, feature in enumerate(db.all_features(order_by=("seqid", "start", "end")), start=1):
        if feature.start is None or feature.end is None:
            continue

        start_i = int(feature.start)
        end_i = int(feature.end)
        if start_i > end_i:
            continue

        attrs = feature.attributes
        feature_id = _first_attr(attrs, ("ID",))
        if not feature_id:
            feature_id = feature.id or f"feature_{lineno}"

        feature_name = _first_attr(attrs, ("Name", "NAME"))
        if not feature_name:
            feature_name = feature_id

        strand = feature.strand if feature.strand in {"+", "-"} else "."
        yield feature.seqid, start_i - 1, end_i, feature_id, feature_name, strand, feature.featuretype


def iter_gff_features(gff_path: Path) -> Iterable[tuple[str, int, int, str, str, str]]:
    cached = GFF_FEATURES_CACHE.get(gff_path)
    if cached is not None:
        yield from cached
        return

    if gffutils is not None:
        rows_raw = list(_iter_gff_features_gffutils(gff_path))
    else:
        rows_raw = list(_iter_gff_features_fallback(gff_path))

    rows = _select_representative_features(rows_raw)

    GFF_FEATURES_CACHE[gff_path] = rows
    yield from rows


def safe_sheet_name(name: str, used: set[str]) -> str:
    clean = re.sub(r"[\\/*?:\[\]]", "_", name)[:31] or "sheet"
    base = clean
    idx = 1
    while clean in used:
        suffix = f"_{idx}"
        clean = (base[: 31 - len(suffix)] + suffix) if len(base) + len(suffix) > 31 else (base + suffix)
        idx += 1
    used.add(clean)
    return clean


def auto_fit_columns(ws) -> None:
    for col_cells in ws.columns:
        max_len = 0
        col_idx = col_cells[0].column
        for cell in col_cells:
            value = "" if cell.value is None else str(cell.value)
            if len(value) > max_len:
                max_len = len(value)
        ws.column_dimensions[get_column_letter(col_idx)].width = min(max_len + 2, 80)


def generate_excel_for_genome(
    genome_id: str,
    tables_dir: Path,
    out_xlsx: Path,
    sra_ids: list[str],
    mapping_counts: dict[str, int],
) -> None:
    tsv_files = sorted(tables_dir.glob("*.tsv"))
    if not tsv_files:
        raise RuntimeError(f"No TSV tables found to convert in {tables_dir}")

    out_xlsx.parent.mkdir(parents=True, exist_ok=True)
    wb = Workbook()
    wb.remove(wb.active)
    used: set[str] = set()

    summary_ws = wb.create_sheet(title="General_mapping")
    summary_ws.append(["SRA_ID", genome_id])
    for sra_id in sra_ids:
        summary_ws.append([sra_id, mapping_counts.get(sra_id, 0)])
    auto_fit_columns(summary_ws)

    for tsv in tsv_files:
        ws = wb.create_sheet(title=safe_sheet_name(tsv.stem, used))
        with tsv.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.reader(handle, delimiter="\t")
            for row in reader:
                ws.append(row)
        auto_fit_columns(ws)

    wb.save(out_xlsx)


# ---------------------------
# Workflow rules
# ---------------------------
# Pipeline order:
# 01 reference indexes -> 02 download reads -> 03 trim/QC -> 04 mapping
# -> 05 mapped counts -> 06-10 feature coverage intermediates
# -> 11 per-sample feature table -> 12 per-genome Excel

rule all:
    input:
        expand("results/excel/{genome_id}.xlsx", genome_id=GENOME_IDS)


rule r01_prepare_reference:
    input:
        fasta=lambda wc: str(GENOME_BY_ID[wc.genome_id]["fasta"])
    output:
        amb=maybe_temp("results/tmp/bwa_index/{genome_id}.amb", KEEP_AUX),
        ann=maybe_temp("results/tmp/bwa_index/{genome_id}.ann", KEEP_AUX),
        bwt=maybe_temp("results/tmp/bwa_index/{genome_id}.bwt", KEEP_AUX),
        pac=maybe_temp("results/tmp/bwa_index/{genome_id}.pac", KEEP_AUX),
        sa=maybe_temp("results/tmp/bwa_index/{genome_id}.sa", KEEP_AUX)
    log:
        "results/logs/r01_prepare_reference/{genome_id}.log"
    threads: 1
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{log}")"
        exec >> "{log}" 2>&1
        mkdir -p results/tmp/bwa_index
        if [[ ! -s results/tmp/bwa_index/{wildcards.genome_id}.bwt ]]; then
            bwa index -p results/tmp/bwa_index/{wildcards.genome_id} {input.fasta}
        fi
        if [[ ! -s {input.fasta}.fai ]]; then
            samtools faidx {input.fasta}
        fi
        """


rule r02_download_reads:
    output:
        raw_single=maybe_temp("results/reads/{sra_id}/{sra_id}.fastq.gz", KEEP_AUX),
        raw_r1=maybe_temp("results/reads/{sra_id}/{sra_id}_1.fastq.gz", KEEP_AUX),
        raw_r2=maybe_temp("results/reads/{sra_id}/{sra_id}_2.fastq.gz", KEEP_AUX)
    log:
        "results/logs/r02_download_reads/{sra_id}.log"
    threads: PREPARE_READS_THREADS
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{log}")"
        exec >> "{log}" 2>&1
        mkdir -p results/reads/{wildcards.sra_id}

        if [[ -s {output.raw_single} ]]; then
            :
        elif [[ -s {output.raw_r1} && -s {output.raw_r2} ]]; then
            :
        else
            rm -f {output.raw_single} {output.raw_r1} {output.raw_r2}
            rm -f results/reads/{wildcards.sra_id}/{wildcards.sra_id}.fastq \
                  results/reads/{wildcards.sra_id}/{wildcards.sra_id}_1.fastq \
                  results/reads/{wildcards.sra_id}/{wildcards.sra_id}_2.fastq
            fasterq-dump {wildcards.sra_id} \
                --outdir results/reads/{wildcards.sra_id} \
                --temp results/reads/ \
                --threads {threads} \
                --skip-technical

            for fq in \
                results/reads/{wildcards.sra_id}/{wildcards.sra_id}.fastq \
                results/reads/{wildcards.sra_id}/{wildcards.sra_id}_1.fastq \
                results/reads/{wildcards.sra_id}/{wildcards.sra_id}_2.fastq; do
                if [[ -s "$fq" ]]; then
                    pigz -f -p {threads} "$fq"
                else
                    rm -f "$fq"
                fi
            done
        fi

        [[ -f {output.raw_single} ]] || : > {output.raw_single}
        [[ -f {output.raw_r1} ]] || : > {output.raw_r1}
        [[ -f {output.raw_r2} ]] || : > {output.raw_r2}
        """


rule r03_prepare_reads:
    input:
        raw_single="results/reads/{sra_id}/{sra_id}.fastq.gz",
        raw_r1="results/reads/{sra_id}/{sra_id}_1.fastq.gz",
        raw_r2="results/reads/{sra_id}/{sra_id}_2.fastq.gz"
    output:
        trimmed_single=maybe_temp("results/qc/{sra_id}/{sra_id}.trimmed.fastq.gz", KEEP_AUX),
        trimmed_r1=maybe_temp("results/qc/{sra_id}/{sra_id}_1.trimmed.fastq.gz", KEEP_AUX),
        trimmed_r2=maybe_temp("results/qc/{sra_id}/{sra_id}_2.trimmed.fastq.gz", KEEP_AUX),
        fastp_html=maybe_temp("results/qc/{sra_id}/{sra_id}.fastp.html", KEEP_AUX),
        fastp_json=maybe_temp("results/qc/{sra_id}/{sra_id}.fastp.json", KEEP_AUX),
        layout=maybe_temp("results/qc/{sra_id}/layout.txt", KEEP_AUX)
    log:
        "results/logs/r03_prepare_reads/{sra_id}.log"
    threads: PREPARE_READS_THREADS
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{log}")"
        exec >> "{log}" 2>&1
        mkdir -p results/qc/{wildcards.sra_id}

        if [[ -s {input.raw_r1} && -s {input.raw_r2} ]]; then
            fastp \
                -i {input.raw_r1} \
                -I {input.raw_r2} \
                -o {output.trimmed_r1} \
                -O {output.trimmed_r2} \
                -h {output.fastp_html} \
                -j {output.fastp_json} \
                -w {threads}
            : > {output.trimmed_single}
            echo "PE" > {output.layout}
        elif [[ -s {input.raw_single} ]]; then
            fastp \
                -i {input.raw_single} \
                -o {output.trimmed_single} \
                -h {output.fastp_html} \
                -j {output.fastp_json} \
                -w {threads}
            : > {output.trimmed_r1}
            : > {output.trimmed_r2}
            echo "SE" > {output.layout}
        else
            echo "Error: no valid FASTQ(.gz) files found for {wildcards.sra_id}" >&2
            exit 1
        fi
        """


rule r04_map_reads:
    input:
        ref_fasta=lambda wc: str(GENOME_BY_ID[wc.genome_id]["fasta"]),
        ref_bwt="results/tmp/bwa_index/{genome_id}.bwt",
        ref_ann="results/tmp/bwa_index/{genome_id}.ann",
        ref_amb="results/tmp/bwa_index/{genome_id}.amb",
        ref_pac="results/tmp/bwa_index/{genome_id}.pac",
        ref_sa="results/tmp/bwa_index/{genome_id}.sa",
        trimmed_single="results/qc/{sra_id}/{sra_id}.trimmed.fastq.gz",
        trimmed_r1="results/qc/{sra_id}/{sra_id}_1.trimmed.fastq.gz",
        trimmed_r2="results/qc/{sra_id}/{sra_id}_2.trimmed.fastq.gz",
        layout="results/qc/{sra_id}/layout.txt"
    output:
        cram=maybe_temp("results/mapping/{genome_id}/{sra_id}.cram", KEEP_MAPPING),
        crai=maybe_temp("results/mapping/{genome_id}/{sra_id}.cram.crai", KEEP_MAPPING)
    log:
        "results/logs/r04_map_reads/{genome_id}/{sra_id}.log"
    threads: MAP_THREADS
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{log}")"
        exec >> "{log}" 2>&1
        mkdir -p results/mapping/{wildcards.genome_id}
        printf "[%s] genome=%s sra_id=%s\n" "$(date '+%Y-%m-%d %H:%M:%S')" "{wildcards.genome_id}" "{wildcards.sra_id}"

        layout=$(tr -d '[:space:]' < {input.layout} | tr '[:lower:]' '[:upper:]')
        if [[ "$layout" == "PE" ]]; then
            (bwa mem -t {threads} results/tmp/bwa_index/{wildcards.genome_id} {input.trimmed_r1} {input.trimmed_r2} \
                | samtools sort -@ {threads} -O CRAM --reference {input.ref_fasta} -o {output.cram})
        elif [[ "$layout" == "SE" ]]; then
            (bwa mem -t {threads} results/tmp/bwa_index/{wildcards.genome_id} {input.trimmed_single} \
                | samtools sort -@ {threads} -O CRAM --reference {input.ref_fasta} -o {output.cram})
        else
            echo "Error: invalid layout '$layout' for {wildcards.sra_id}" >&2
            exit 1
        fi

        samtools index {output.cram}
        """


rule r05_mapped_count:
    input:
        cram="results/mapping/{genome_id}/{sra_id}.cram",
        ref_fasta=lambda wc: str(GENOME_BY_ID[wc.genome_id]["fasta"])
    output:
        maybe_temp("results/counts/{genome_id}/{sra_id}.txt", KEEP_AUX)
    log:
        "results/logs/r05_mapped_count/{genome_id}/{sra_id}.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{log}")"
        exec >> "{log}" 2>&1
        mkdir -p results/counts/{wildcards.genome_id}
        samtools view -c -F 260 -T {input.ref_fasta} {input.cram} > {output}
        """


rule r06_mapped_primary_cram:
    input:
        cram="results/mapping/{genome_id}/{sra_id}.cram",
        ref_fasta=lambda wc: str(GENOME_BY_ID[wc.genome_id]["fasta"])
    output:
        maybe_temp("results/tmp/feature_tables/{genome_id}/{sra_id}/mapped_primary.cram", KEEP_AUX)
    log:
        "results/logs/r06_mapped_primary_cram/{genome_id}/{sra_id}.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{log}")"
        exec >> "{log}" 2>&1
        mkdir -p results/tmp/feature_tables/{wildcards.genome_id}/{wildcards.sra_id}
        samtools view -C -F 260 -T {input.ref_fasta} -o {output} {input.cram}
        """


rule r07_mapped_primary_index:
    input:
        "results/tmp/feature_tables/{genome_id}/{sra_id}/mapped_primary.cram"
    output:
        maybe_temp("results/tmp/feature_tables/{genome_id}/{sra_id}/mapped_primary.cram.crai", KEEP_AUX)
    log:
        "results/logs/r07_mapped_primary_index/{genome_id}/{sra_id}.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{log}")"
        exec >> "{log}" 2>&1
        samtools index {input}
        """


rule r08_reads_bed:
    input:
        mapped_primary="results/tmp/feature_tables/{genome_id}/{sra_id}/mapped_primary.cram",
        mapped_primary_crai="results/tmp/feature_tables/{genome_id}/{sra_id}/mapped_primary.cram.crai",
        ref_fasta=lambda wc: str(GENOME_BY_ID[wc.genome_id]["fasta"])
    output:
        maybe_temp("results/tmp/feature_tables/{genome_id}/{sra_id}/reads.bed", KEEP_AUX)
    log:
        "results/logs/r08_reads_bed/{genome_id}/{sra_id}.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{log}")"
        exec >> "{log}" 2>&1
        samtools view -b -T {input.ref_fasta} {input.mapped_primary} \
            | bedtools bamtobed -i - > {output}
        """


rule r09_features_bed:
    input:
        gff=lambda wc: str(GENOME_BY_ID[wc.genome_id]["gff3"])
    output:
        maybe_temp("results/tmp/feature_tables/{genome_id}/{sra_id}/features_info.bed", KEEP_AUX)
    log:
        "results/logs/r09_features_bed/{genome_id}/{sra_id}.log"
    run:
        log_path = Path(str(log[0]))
        log_path.parent.mkdir(parents=True, exist_ok=True)
        with log_path.open("a", encoding="utf-8") as log_handle, contextlib.redirect_stdout(log_handle), contextlib.redirect_stderr(log_handle):
            print(f"Building feature BED for genome={wildcards.genome_id}, sra_id={wildcards.sra_id}")
            out_path = Path(output[0])
            out_path.parent.mkdir(parents=True, exist_ok=True)
            with out_path.open("w", encoding="utf-8") as out:
                for row in iter_gff_features(Path(input.gff)):
                    out.write("\t".join(map(str, row)) + "\n")


rule r10_feature_coverage:
    input:
        features="results/tmp/feature_tables/{genome_id}/{sra_id}/features_info.bed",
        reads="results/tmp/feature_tables/{genome_id}/{sra_id}/reads.bed"
    output:
        maybe_temp("results/tmp/feature_tables/{genome_id}/{sra_id}/features_coverage.tsv", KEEP_AUX)
    log:
        "results/logs/r10_feature_coverage/{genome_id}/{sra_id}.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{log}")"
        exec >> "{log}" 2>&1
        bedtools coverage -a {input.features} -b {input.reads} > {output}
        """


rule r11_feature_table:
    input:
        coverage="results/tmp/feature_tables/{genome_id}/{sra_id}/features_coverage.tsv"
    output:
        maybe_temp("results/gene_tables/{genome_id}/{sra_id}.tsv", KEEP_AUX)
    log:
        "results/logs/r11_feature_table/{genome_id}/{sra_id}.log"
    run:
        log_path = Path(str(log[0]))
        log_path.parent.mkdir(parents=True, exist_ok=True)
        with log_path.open("a", encoding="utf-8") as log_handle, contextlib.redirect_stdout(log_handle), contextlib.redirect_stderr(log_handle):
            print(f"Building feature table for genome={wildcards.genome_id}, sra_id={wildcards.sra_id}")
            out_path = Path(output[0])
            out_path.parent.mkdir(parents=True, exist_ok=True)

            with out_path.open("w", encoding="utf-8", newline="") as out_handle:
                writer = csv.writer(out_handle, delimiter="\t")
                writer.writerow(
                    [
                        "SRA_ID",
                        "Number of mapped reads",
                        "FeatureID",
                        "Genome coordinates",
                        "Feature name",
                        "Coverage over feature",
                    ]
                )

                with Path(input.coverage).open("r", encoding="utf-8") as cov_handle:
                    for line in cov_handle:
                        cols = line.rstrip("\n").split("\t")
                        if len(cols) < 10:
                            continue
                        seqid, start0, end1, feature_id, feature_name, strand = cols[:6]
                        reads_int = int(cols[6])
                        coverage_fraction = cols[9]
                        coord = f"{seqid}:{int(start0)+1}-{end1}({strand})"

                        if not INCLUDE_ZERO_ROWS and reads_int == 0:
                            continue

                        writer.writerow(
                            [
                                wildcards.sra_id,
                                reads_int,
                                feature_id,
                                coord,
                                feature_name,
                                coverage_fraction,
                            ]
                        )


rule r12_excel_per_genome:
    input:
        tables=expand("results/gene_tables/{{genome_id}}/{sra_id}.tsv", sra_id=SRA_IDS),
        counts=expand("results/counts/{{genome_id}}/{sra_id}.txt", sra_id=SRA_IDS)
    output:
        "results/excel/{genome_id}.xlsx"
    log:
        "results/logs/r12_excel_per_genome/{genome_id}.log"
    run:
        log_path = Path(str(log[0]))
        log_path.parent.mkdir(parents=True, exist_ok=True)
        with log_path.open("a", encoding="utf-8") as log_handle, contextlib.redirect_stdout(log_handle), contextlib.redirect_stderr(log_handle):
            print(f"Building Excel report for genome={wildcards.genome_id}")
            mapping_counts: dict[str, int] = {}
            counts_base = Path("results/counts") / wildcards.genome_id
            for sra_id in SRA_IDS:
                count_file = counts_base / f"{sra_id}.txt"
                mapping_counts[sra_id] = int(count_file.read_text(encoding="utf-8").strip() or "0")

            generate_excel_for_genome(
                genome_id=wildcards.genome_id,
                tables_dir=Path("results/gene_tables") / wildcards.genome_id,
                out_xlsx=Path(output[0]),
                sra_ids=SRA_IDS,
                mapping_counts=mapping_counts,
            )
