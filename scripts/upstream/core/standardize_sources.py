#!/usr/bin/env python3
"""Prepare HM3/no-MHC LDSC inputs for the ascertainment pilot.

The script keeps one traceable first record per HM3 SNP, removes autosomal
non-SNP alleles, strand-ambiguous alleles and low-information variants, and
writes both a full standardized table and an LDSC table with Z and N.
"""

from __future__ import annotations

import argparse
import gzip
import math
import os
import re
from collections import Counter
from typing import Dict, Iterable, Optional, Sequence


CONFIG = {
    "MDD_CLIN_PGC2025": dict(snp=("ID",), chrom=("CHROM",), bp=("POS",), a1=("EA",), a2=("NEA",), beta=("BETA",), se=("SE",), p=("PVAL",), n=("NEFF",), info=("IMPINFO",), build="GRCh37"),
    "MDD_EHR_PGC2025": dict(snp=("ID",), chrom=("CHROM",), bp=("POS",), a1=("EA",), a2=("NEA",), beta=("BETA",), se=("SE",), p=("PVAL",), n=("NEFF",), info=("IMPINFO",), build="GRCh37"),
    "MDD_QUEST_PGC2025": dict(snp=("ID",), chrom=("CHROM",), bp=("POS",), a1=("EA",), a2=("NEA",), beta=("BETA",), se=("SE",), p=("PVAL",), n=("NEFF",), info=("IMPINFO",), build="GRCh37"),
    "MDD_FINNGEN_R13": dict(snp=("SNP",), chrom=("CHR",), bp=("BP",), a1=("effect_allele",), a2=("other_allele",), beta=("beta",), se=("se",), p=("pval",), n=("N",), eaf=("eaf",), build="GRCh38"),
    "BD_PGC2021": dict(snp=("ID",), chrom=("CHROM",), bp=("POS",), a1=("A1",), a2=("A2",), beta=("BETA",), se=("SE",), p=("PVAL",), n=("NEFFDIV2",), n_factor=2.0, info=("IMPINFO",), build="GRCh37"),
    "BD_CLIN_PGC4": dict(snp=("SNP",), chrom=("CHR",), bp=("BP",), a1=("A1",), a2=("A2",), beta=("BETA",), odds_ratio=("OR",), se=("SE",), p=("P",), n=("Neff_half",), n_factor=2.0, info=("INFO",), build="GRCh37"),
    "BD_FINNGEN_R13": dict(snp=("RSIDS",), chrom=("CHROM",), bp=("POS",), a1=("ALT",), a2=("REF",), beta=("BETA",), se=("SEBETA",), p=("PVAL",), n_fixed=36114.340, eaf=("AF_ALT",), build="GRCh38"),
    "SCZ_PGC2022": dict(snp=("ID",), chrom=("CHROM",), bp=("POS",), a1=("A1",), a2=("A2",), beta=("BETA",), se=("SE",), p=("PVAL",), n=("NEFF",), info=("IMPINFO",), build="GRCh37"),
    "SCZ_FINNGEN_R13": dict(snp=("RSIDS",), chrom=("CHROM",), bp=("POS",), a1=("ALT",), a2=("REF",), beta=("BETA",), se=("SEBETA",), p=("PVAL",), n_fixed=28812.977, eaf=("AF_ALT",), build="GRCh38"),
    "AD_GCST90704646_MAIN": dict(snp=("rs_id", "variant_id"), chrom=("chromosome",), bp=("base_pair_location",), a1=("effect_allele",), a2=("other_allele",), beta=("beta",), se=("standard_error",), p=("p_value",), n=("Neff_total",), eaf=("effect_allele_frequency",), build="GRCh38"),
    "AD_GCST90704647_NOPROXY": dict(snp=("rs_id", "variant_id"), chrom=("chromosome",), bp=("base_pair_location",), a1=("effect_allele",), a2=("other_allele",), beta=("beta",), se=("standard_error",), p=("p_value",), n=("Neff_total",), eaf=("effect_allele_frequency",), build="GRCh38"),
    "AD_GCST90704648_NOBIOBANK": dict(snp=("rs_id", "variant_id"), chrom=("chromosome",), bp=("base_pair_location",), a1=("effect_allele",), a2=("other_allele",), beta=("beta",), se=("standard_error",), p=("p_value",), n=("Neff_total",), eaf=("effect_allele_frequency",), build="GRCh38"),
    "AD_GCST007511_KUNKLE": dict(snp=("MarkerName",), chrom=("Chromosome",), bp=("Position",), a1=("Effect_allele",), a2=("Non_Effect_allele",), beta=("Beta",), se=("SE",), p=("Pvalue",), n_fixed=57692.5199762225, build="GRCh37"),
    "EPILEPSY_ILAE2023_EUR": dict(snp=("rs_id",), chrom=("chromosome",), bp=("base_pair_location",), a1=("effect_allele",), a2=("other_allele",), beta=("beta",), se=("standard_error",), p=("p_value",), n_fixed=66832.986585, eaf=("effect_allele_frequency",), build="GRCh37", beta_is_z=True),
    "EPILEPSY_MVP2024_EUR": dict(snp=("SNP_ID",), chrom=("chrom",), bp=("pos",), a1=("ea",), a2=("ref",), alt=("alt",), odds_ratio=("or",), ci=("ci",), p=("pval",), n_fixed=24332.120012, eaf=("af",), build="UNVERIFIED", a2_from_alt=True, se_from_ci=True),
    "EPILEPSY_FINNGEN_R13": dict(snp=("RSIDS",), chrom=("CHROM",), bp=("POS",), a1=("ALT",), a2=("REF",), beta=("BETA",), se=("SEBETA",), p=("PVAL",), n_fixed=62405.792541, eaf=("AF_ALT",), build="GRCh38"),
    "EPILEPSY_UKB_X345_GCST90435926": dict(snp=("variant_id",), chrom=("chromosome",), bp=("base_pair_location",), a1=("effect_allele",), a2=("other_allele",), beta=("beta",), se=("standard_error",), p=("p_value",), n_fixed=20089.415662409814, eaf=("effect_allele_frequency",), build="GRCh37"),
    "EPILEPSY_EHR_META_2OF3": dict(snp=("SNP",), chrom=("CHR",), bp=("BP",), a1=("A1",), a2=("A2",), beta=("BETA",), se=("SE",), p=("P",), n=("N",), build="UNVERIFIED"),
    "EPILEPSY_EHR_META_3OF3": dict(snp=("SNP",), chrom=("CHR",), bp=("BP",), a1=("A1",), a2=("A2",), beta=("BETA",), se=("SE",), p=("P",), n=("N",), build="UNVERIFIED"),
    "PD_FINNGEN_R13": dict(snp=("RSIDS",), chrom=("CHROM",), bp=("POS",), a1=("ALT",), a2=("REF",), beta=("BETA",), se=("SEBETA",), p=("PVAL",), n_fixed=24960.579624, eaf=("AF_ALT",), build="GRCh38"),
    "NEURO_PD_GCST009324": dict(snp=("rsid", "hm_rsid", "variant_id"), chrom=("chromosome", "hm_chrom"), bp=("base_pair_location", "hm_pos"), a1=("effect_allele", "hm_effect_allele"), a2=("other_allele", "hm_other_allele"), beta=("beta", "hm_beta"), se=("standard_error",), p=("p_value",), n_fixed=27481.700617, eaf=("effect_allele_frequency", "hm_effect_allele_frequency"), build="GRCh38"),
    "NEURO_ALS_GCST90027164": dict(snp=("hm_rsid", "rsid", "variant_id"), chrom=("hm_chrom", "chromosome"), bp=("hm_pos", "base_pair_location"), a1=("hm_effect_allele", "effect_allele"), a2=("hm_other_allele", "other_allele"), beta=("hm_beta", "beta"), se=("standard_error",), p=("p_value",), n_fixed=87380.842082, eaf=("hm_effect_allele_frequency", "effect_allele_frequency"), build="GRCh38"),
    "NEURO_LBD_GCST90001390": dict(snp=("hm_rsid", "rsid", "variant_id"), chrom=("hm_chrom", "chromosome"), bp=("hm_pos", "base_pair_location"), a1=("hm_effect_allele", "effect_allele"), a2=("hm_other_allele", "other_allele"), beta=("hm_beta", "beta"), se=("standard_error",), p=("p_value",), n_fixed=6306.411001, eaf=("hm_effect_allele_frequency", "effect_allele_frequency"), build="GRCh38"),
    "NEURO_MS_GCST90558092": dict(snp=("rsid",), chrom=("chromosome",), bp=("base_pair_location",), a1=("effect_allele",), a2=("other_allele",), beta=("beta",), se=("standard_error",), p=("p_value",), n_fixed=62542.506185, eaf=("effect_allele_frequency",), build="GRCh37", coordinate_snp=True),
    "NEURO_MIGRAINE_GCST90271641": dict(snp=("rsid",), chrom=("chromosome",), bp=("base_pair_location",), a1=("effect_allele",), a2=("other_allele",), beta=("beta",), se=("standard_error",), p=("p_value",), n_fixed=98918.684800, eaf=("effect_allele_frequency",), build="GRCh37", coordinate_snp=True),
    "NEURO_ISCHEMIC_STROKE_GCST90104540": dict(snp=("rsid", "hm_rsid", "variant_id"), chrom=("chromosome", "hm_chrom"), bp=("base_pair_location", "hm_pos"), a1=("effect_allele", "hm_effect_allele"), a2=("other_allele", "hm_other_allele"), beta=("beta", "hm_beta"), se=("standard_error",), p=("p_value",), n_fixed=236505.811112, eaf=("effect_allele_frequency", "hm_effect_allele_frequency"), build="GRCh38"),
    "NEURO_RLS_GCST90435387": dict(snp=("rsid", "hm_rsid", "variant_id"), chrom=("chromosome", "hm_chrom"), bp=("base_pair_location", "hm_pos"), a1=("effect_allele", "hm_effect_allele"), a2=("other_allele", "hm_other_allele"), beta=("beta", "hm_beta"), se=("standard_error",), p=("p_value",), n_fixed=31451.024318, eaf=("effect_allele_frequency", "hm_effect_allele_frequency"), build="GRCh38"),
}


def norm(value: str) -> str:
    return re.sub(r"[^a-z0-9]", "", value.lower())


def number(value: Optional[str]) -> Optional[float]:
    if value is None:
        return None
    value = value.strip().strip('"').replace(",", "")
    if value in {"", ".", "NA", "NaN", "nan", "NULL", "null", "-"}:
        return None
    try:
        result = float(value)
    except ValueError:
        return None
    return result if math.isfinite(result) else None


def open_text(path: str):
    if path.endswith((".gz", ".bgz")):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="")
    return open(path, "rt", encoding="utf-8", errors="replace", newline="")


def read_header(handle) -> tuple[list[str], str]:
    for line in handle:
        clean = line.rstrip("\r\n")
        if not clean or clean.startswith("##"):
            continue
        clean = clean.lstrip("#")
        if "\t" in clean:
            return clean.split("\t"), "\t"
        return clean.split(), " "
    raise ValueError("no header found")


def index(header: Sequence[str], aliases: Iterable[str]) -> Optional[int]:
    lookup = {norm(name): i for i, name in enumerate(header)}
    for alias in aliases:
        if norm(alias) in lookup:
            return lookup[norm(alias)]
    return None


def value(row: Sequence[str], position: Optional[int]) -> str:
    return row[position] if position is not None and position < len(row) else ""


def hm3_no_mhc_reference(path: str) -> tuple[set[str], dict[tuple[int, int], str]]:
    ids: set[str] = set()
    coordinates: dict[tuple[int, int], str] = {}
    if not os.path.isdir(path):
        raise FileNotFoundError(path)
    for name in sorted(os.listdir(path)):
        if not name.endswith(".l2.ldscore.gz"):
            continue
        with gzip.open(os.path.join(path, name), "rt", encoding="utf-8", errors="replace") as handle:
            next(handle, None)
            for line in handle:
                fields = line.rstrip("\r\n").split("\t")
                if len(fields) >= 3 and re.fullmatch(r"rs[0-9]+", fields[1], re.I):
                    snp = fields[1].lower()
                    ids.add(snp)
                    try:
                        coordinates[(int(fields[0]), int(fields[2]))] = snp
                    except ValueError:
                        pass
    if not ids:
        raise RuntimeError("no HM3/no-MHC SNP IDs found in " + path)
    if not coordinates:
        raise RuntimeError("no HM3/no-MHC coordinates found in " + path)
    return ids, coordinates


def first_hm3_rsid(raw: str, hm3: set[str]) -> str:
    for candidate in re.findall(r"rs[0-9]+", raw or "", flags=re.I):
        candidate = candidate.lower()
        if candidate in hm3:
            return candidate
    return ""


def fmt(value: float) -> str:
    return "%.12g" % value


def ci_bounds(raw: str) -> tuple[Optional[float], Optional[float]]:
    """Parse a two-sided confidence interval such as '0.8,1.2'."""
    text = (raw or "").strip().strip("[]()")
    parts = re.split(r"\s*,\s*", text)
    if len(parts) != 2:
        return None, None
    low, high = number(parts[0]), number(parts[1])
    if low is None or high is None or low <= 0 or high <= 0 or high < low:
        return None, None
    return low, high


def run(args: argparse.Namespace) -> None:
    cfg = CONFIG[args.key]
    hm3, hm3_by_coordinate = hm3_no_mhc_reference(args.hm3_dir)
    counters = Counter()
    seen: set[str] = set()
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    os.makedirs(os.path.dirname(args.ldsc_output), exist_ok=True)
    os.makedirs(os.path.dirname(args.audit), exist_ok=True)
    output_part = args.output + ".part"
    ldsc_part = args.ldsc_output + ".part"
    try:
        with open_text(args.input) as src, gzip.open(output_part, "wt", encoding="utf-8", newline="") as full, gzip.open(ldsc_part, "wt", encoding="utf-8", newline="") as ldsc:
            header, sep = read_header(src)
            positions = {key: index(header, cfg.get(key, ())) for key in ("snp", "chrom", "bp", "a1", "a2", "alt", "beta", "odds_ratio", "ci", "se", "p", "n", "info", "eaf")}
            required = ["chrom", "bp", "a1", "a2", "p"]
            if not cfg.get("coordinate_snp"):
                required.insert(0, "snp")
            if cfg.get("beta_is_z"):
                required.append("beta")
            elif positions["beta"] is None and positions["odds_ratio"] is None:
                required.append("beta or odds_ratio")
            elif positions["se"] is None and positions["ci"] is None:
                required.append("se")
            missing = [key for key in required if positions[key] is None]
            if cfg.get("n_fixed") is None and positions["n"] is None:
                missing.append("n or n_fixed")
            if missing:
                raise ValueError("{}: missing columns {}; header={}".format(args.key, missing, header))
            full.write("SNP\tCHR\tBP\tA1\tA2\tBETA\tSE\tP\tN\n")
            ldsc.write("SNP\tA1\tA2\tZ\tN\n")
            for raw in src:
                line = raw.rstrip("\r\n")
                if not line:
                    continue
                counters["rows_read"] += 1
                row = line.split(sep) if sep == "\t" else line.split()
                if len(row) != len(header):
                    counters["malformed_rows"] += 1
                    continue
                a1 = value(row, positions["a1"]).strip().strip('"').upper()
                a2 = value(row, positions["a2"]).strip().strip('"').upper()
                if cfg.get("a2_from_alt") and a1 == a2:
                    a2 = value(row, positions["alt"]).strip().upper()
                if not re.fullmatch(r"[ACGT]", a1) or not re.fullmatch(r"[ACGT]", a2) or a1 == a2:
                    counters["invalid_or_indel_alleles"] += 1
                    continue
                if {a1, a2} in ({"A", "T"}, {"C", "G"}):
                    counters["strand_ambiguous"] += 1
                    continue
                chrom_raw = value(row, positions["chrom"]).strip().strip('"').lower().removeprefix("chr")
                bp_raw = value(row, positions["bp"]).strip().strip('"')
                try:
                    chrom = int(float(chrom_raw))
                    bp = int(float(bp_raw)) + int(cfg.get("bp_offset", 0))
                except ValueError:
                    counters["invalid_coordinate"] += 1
                    continue
                if not 1 <= chrom <= 22 or bp < 1:
                    counters["non_autosomal_or_invalid_coordinate"] += 1
                    continue
                snp = first_hm3_rsid(value(row, positions["snp"]), hm3)
                if not snp and cfg["build"] == "GRCh37":
                    snp = hm3_by_coordinate.get((chrom, bp), "")
                    if snp:
                        counters["hm3_rsid_from_grch37_coordinate"] += 1
                if not snp:
                    counters["not_hm3_or_no_rsid"] += 1
                    continue
                if snp in seen:
                    counters["duplicate_hm3_snp"] += 1
                    continue
                if cfg["build"] == "GRCh37" and chrom == 6 and 28_477_797 <= bp <= 33_448_354:
                    counters["mhc_removed_by_coordinate"] += 1
                    continue
                beta = number(value(row, positions["beta"]))
                if beta is None and positions["odds_ratio"] is not None:
                    odds_ratio = number(value(row, positions["odds_ratio"]))
                    if odds_ratio is not None and odds_ratio > 0:
                        beta = math.log(odds_ratio)
                        counters["beta_from_log_or"] += 1
                se = number(value(row, positions["se"]))
                p = number(value(row, positions["p"]))
                if cfg.get("n_fixed") is not None:
                    n = float(cfg["n_fixed"])
                else:
                    n = number(value(row, positions["n"]))
                    if n is not None:
                        n *= float(cfg.get("n_factor", 1.0))
                info = number(value(row, positions["info"])) if positions["info"] is not None else None
                if info is not None and info < 0.9:
                    counters["info_below_0_9"] += 1
                    continue
                eaf = number(value(row, positions["eaf"])) if positions["eaf"] is not None else None
                if eaf is not None and not 0 <= eaf <= 1:
                    counters["invalid_eaf"] += 1
                    continue
                if eaf is not None and min(eaf, 1.0 - eaf) < 0.01:
                    counters["maf_below_0_01"] += 1
                    continue
                if cfg.get("beta_is_z"):
                    z = beta
                    counters["z_from_source_beta"] += 1 if z is not None else 0
                    invalid_stats = z is None or p is None or p <= 0 or p > 1 or n is None or n <= 0
                else:
                    if se is None and cfg.get("se_from_ci"):
                        low, high = ci_bounds(value(row, positions["ci"]))
                        if low is not None and high is not None:
                            se = (math.log(high) - math.log(low)) / (2.0 * 1.959963984540054)
                            counters["se_from_95ci"] += 1
                    if p is not None and p <= 0 and beta is not None and se is not None and se > 0:
                        p = math.erfc(abs(beta / se) / math.sqrt(2.0))
                        counters["p_recalculated_from_beta_se"] += 1
                    invalid_stats = beta is None or se is None or se <= 0 or p is None or p <= 0 or p > 1 or n is None or n <= 0
                    z = beta / se if not invalid_stats else None
                if invalid_stats:
                    counters["invalid_core_statistics"] += 1
                    continue
                if not math.isfinite(z):
                    counters["invalid_z"] += 1
                    continue
                seen.add(snp)
                counters["rows_written"] += 1
                se_text = fmt(se) if se is not None else "NA"
                full.write("\t".join((snp, str(chrom), str(bp), a1, a2, fmt(beta), se_text, fmt(p), fmt(n))) + "\n")
                ldsc.write("\t".join((snp, a1, a2, fmt(z), fmt(n))) + "\n")
    except Exception:
        for path in (output_part, ldsc_part):
            if os.path.exists(path):
                os.unlink(path)
        raise
    os.replace(output_part, args.output)
    os.replace(ldsc_part, args.ldsc_output)
    with open(args.audit, "w", encoding="utf-8", newline="") as handle:
        handle.write("metric\tvalue\n")
        for key, val in sorted(counters.items()):
            handle.write("{}\t{}\n".format(key, val))
        handle.write("key\t{}\ninput\t{}\noutput\t{}\nldsc_output\t{}\nhm3_no_mhc_dir\t{}\n".format(args.key, args.input, args.output, args.ldsc_output, args.hm3_dir))
        settings = "HM3/no-MHC;autosomes 1-22;single-nucleotide alleles;strand-ambiguous removed;INFO>=0.9 when present;MAF>=0.01 when present;MHC hg19 coordinate exclusion for GRCh37"
        if cfg.get("beta_is_z"):
            settings += ";source beta interpreted as Z;source SE is missing by design;fixed effective N"
        if cfg.get("se_from_ci"):
            settings += ";beta=log(OR);SE derived from reported 95% CI when source SE is unavailable;fixed effective N"
        handle.write("settings\t{}\n".format(settings))
    print("STANDARDIZE_COMPLETE\t{}\t{}".format(args.key, counters["rows_written"]), flush=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--key", required=True, choices=sorted(CONFIG))
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--ldsc-output", required=True)
    parser.add_argument("--audit", required=True)
    parser.add_argument("--hm3-dir", required=True)
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
