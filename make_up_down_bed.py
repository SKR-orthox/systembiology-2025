from __future__ import annotations

from pathlib import Path
import pandas as pd

# ====== Settings (edit here) ======
INPUT_NAME = "input.txt"  # input filename
FDR_CUTOFF = 0.05                  # set 1.0 to disable FDR filtering
ABS_LFC_CUTOFF = 1.0               # set 0 to split only by sign (+/-)
ASSUME_1BASED_START = True         # True if Start is 1-based (BED is 0-based)
# ================================


def find_col_contains(
    df: pd.DataFrame,
    keywords: list[str],
    prefer: list[str] | None = None
) -> str | None:
    """Return the first column name that contains any keyword (case-insensitive)."""
    cols = list(df.columns)
    low = [c.lower() for c in cols]

    def _match(keys: list[str]) -> str | None:
        keys = [k.lower() for k in keys]
        for c, cl in zip(cols, low):
            if any(k in cl for k in keys):
                return c
        return None

    if prefer:
        hit = _match(prefer)
        if hit:
            return hit
    return _match(keywords)


def normalize_chr_ucsc(v) -> str:
    """
    Normalize chromosome names to UCSC style:
    - '1' -> 'chr1'
    - 'chr1' -> 'chr1'
    - 'X'/'chrX' -> 'chrX'
    - 'Y'/'chrY' -> 'chrY'
    - 'M'/'MT'/'chrM' -> 'chrM'
    """
    if pd.isna(v):
        return v
    s = str(v).strip()
    if s == "":
        return s

    sl = s.lower()

    # Already UCSC-like
    if sl.startswith("chr"):
        rest = s[3:]
        if rest.upper() in {"M", "MT"}:
            return "chrM"
        if rest.upper() == "X":
            return "chrX"
        if rest.upper() == "Y":
            return "chrY"
        if rest.isdigit():
            return "chr" + str(int(rest))
        return s

    # Numeric
    if s.isdigit():
        return "chr" + str(int(s))

    # Sex chromosomes
    if s.upper() == "X":
        return "chrX"
    if s.upper() == "Y":
        return "chrY"

    # Mitochondria
    if s.upper() in {"M", "MT"}:
        return "chrM"

    return s


def main():
    # 入力ファイル探索: カレントディレクトリ優先、なければスクリプト位置
    cwd = Path.cwd()
    script_dir = Path(__file__).resolve().parent
    in_path = cwd / INPUT_NAME
    if not in_path.exists():
        in_path = script_dir / INPUT_NAME
    if not in_path.exists():
        raise FileNotFoundError(
            f"Input file not found: {INPUT_NAME}\n"
            f"Checked: {cwd} or {script_dir}"
        )

    df = pd.read_csv(in_path, sep="\t", dtype=str)

    # 必須列を自動検出（SeqMonk/DESeq2出力の表を想定）
    chr_col = find_col_contains(df, ["chromosome", "chr"], prefer=["chromosome"])
    start_col = find_col_contains(df, ["start"], prefer=["start"])
    end_col = find_col_contains(df, ["end"], prefer=["end"])
    fdr_col = find_col_contains(df, ["fdr", "padj", "qvalue", "q-value"], prefer=["fdr"])
    lfc_col = find_col_contains(
        df,
        keywords=["log2 fold change", "log2foldchange", "log2fc"],
        prefer=["shrunk log2 fold change", "shrunken log2 fold change", "shrunk"]
    )

    missing = [name for name, col in {
        "Chromosome": chr_col,
        "Start": start_col,
        "End": end_col,
        "log2FC": lfc_col
    }.items() if col is None]
    if missing:
        raise ValueError(
            "Required columns not found: " + ", ".join(missing) + "\n"
            "Please check the header names."
        )

    # Chromosome name normalization (UCSC)
    df[chr_col] = df[chr_col].apply(normalize_chr_ucsc)

    # Convert to numeric
    df[start_col] = pd.to_numeric(df[start_col], errors="coerce")
    df[end_col] = pd.to_numeric(df[end_col], errors="coerce")
    df[lfc_col] = pd.to_numeric(df[lfc_col], errors="coerce")
    if fdr_col is not None:
        df[fdr_col] = pd.to_numeric(df[fdr_col], errors="coerce")

    # Drop rows missing required values
    df = df.dropna(subset=[chr_col, start_col, end_col, lfc_col])
    if fdr_col is not None:
        df = df.dropna(subset=[fdr_col])

    # FDR filtering (optional)
    if fdr_col is not None and FDR_CUTOFF < 1.0:
        df = df[df[fdr_col] <= FDR_CUTOFF]

    # Split UP/DOWN by log2FC
    if ABS_LFC_CUTOFF > 0:
        up = df[df[lfc_col] >= ABS_LFC_CUTOFF].copy()
        dn = df[df[lfc_col] <= -ABS_LFC_CUTOFF].copy()
    else:
        up = df[df[lfc_col] > 0].copy()
        dn = df[df[lfc_col] < 0].copy()

    # BED conversion
    def to_bed(x: pd.DataFrame) -> pd.DataFrame:
        bed = x[[chr_col, start_col, end_col]].copy()

        # BED uses 0-based start
        if ASSUME_1BASED_START:
            bed[start_col] = bed[start_col] - 1

        bed[start_col] = bed[start_col].astype(int)
        bed[end_col] = bed[end_col].astype(int)

        # Fix negative start and invalid intervals
        bed[start_col] = bed[start_col].clip(lower=0)
        bed.columns = ["chr", "start", "end"]
        bed = bed[bed["end"] > bed["start"]]
        return bed

    up_bed = to_bed(up)
    dn_bed = to_bed(dn)

    # Output
    stem = in_path.stem
    out_dir = in_path.parent
    up_path = out_dir / f"{stem}__UP.bed"
    dn_path = out_dir / f"{stem}__DOWN.bed"

    up_bed.to_csv(up_path, sep="\t", index=False, header=False)
    dn_bed.to_csv(dn_path, sep="\t", index=False, header=False)

    print("=== DONE ===")
    print(f"Input:      {in_path}")
    print(f"Chr col:    {chr_col}")
    print(f"Start/End:  {start_col}, {end_col}")
    print(f"log2FC col: {lfc_col}")
    print(f"FDR col:    {fdr_col}")
    print(f"Filter:     FDR<={FDR_CUTOFF} |LFC|>={ABS_LFC_CUTOFF}  1-based start={ASSUME_1BASED_START}")
    print(f"UP BED:     {up_path} (n={len(up_bed)})")
    print(f"DOWN BED:   {dn_path} (n={len(dn_bed)})")
    print("Chromosome normalized to UCSC style: chr1.. chrX chrY chrM")


if __name__ == "__main__":
    main()