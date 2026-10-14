"""Unit tests of call_dmr and annotate (no R needed)."""
import gzip

import numpy as np
import pandas as pd

from dmrdetector.steps import annotate
from dmrdetector.steps import call_dmr
from dmrdetector.util import rstats


def test_r_round_is_r4_not_numpy():
    # R >= 4.0.0: round(0.15, 1) = 0.1, round(0.25, 1) = 0.2, round(0.35, 1) =
    # 0.3, round(0.45, 1) = 0.4 (checked with R 4.4)
    assert list(rstats.r_round(np.array([0.15, 0.25, 0.35, 0.45, -0.15]),
                               1)) == [0.1, 0.2, 0.3, 0.4, -0.1]
    assert list(rstats.r_round(np.array([0.5, 1.5, 2.5]))) == [0.0, 2.0, 2.0]


def _sites(core, pos, diff):
    return pd.DataFrame({
        "CHROM": "chr1",
        "START": pos,
        "CORE_START": core,
        "COVERAGE": 10,
        "SAMPLE_METH": 0.9,
        "REF_METH": 0.4,
        "DIFF": diff
    })


def test_cpg_cluster_uses_core_start_in_site_order():
    # two nested regions (core 1000 and 5000) interleave in site order; v1 gap
    # uses core START and only a forward jump > 150 starts a new cluster
    df = _sites([1000, 5000, 1000, 5000], [1100, 1100, 1200, 1200], 0.5)
    m = call_dmr.merge_sites(df, 150, pos_col="CORE_START")
    # gaps: Inf, +4000 (new), -4000 (same), +4000 (new) -> {1000}, {5000, 1000},
    #   {5000}
    assert list(m["NUM_CPG"]) == [1, 2, 1]
    assert list(m["START"]) == [1000, 1000, 5000] and list(
        m["END"]) == [1001, 5001, 5001]


def test_empty_result_writes_header_only(tmp_path):
    p = tmp_path / "raw.tsv.gz"
    call_dmr.write_raw(pd.DataFrame({c: [] for c in call_dmr.RAW_COLS}), str(p))
    assert gzip.open(p, "rt").read() == "\t".join(call_dmr.RAW_COLS) + "\n"


GTF = "\n".join([
    "#!genome-build GRCh38",
    '1\tx\tgene\t100\t200\t.\t+\t.\tgene_id "G1"; gene_name "LNC"; '
    'gene_biotype "lncRNA";',
    '1\tx\tgene\t150\t400\t.\t+\t.\tgene_id "G2"; gene_name "PC"; '
    'gene_biotype "protein_coding";',
    '1\tx\tgene\t150\t400\t.\t-\t.\tgene_id "G0"; gene_biotype "weird";',
    '1\tx\texon\t150\t400\t.\t-\t.\tgene_id "G9"; gene_name "EXON";',
    '1\tx\tgene\t1000\t1100\t.\t+\t.\tgene_id "G3"; gene_name "A"; '
    'gene_biotype "miRNA";',
    '1\tx\tgene\t1000\t1100\t.\t+\t.\tgene_id "G4"; gene_name "B"; '
    'gene_biotype "miRNA";',
    'MT\tx\tgene\t1\t100\t.\t+\t.\tgene_id "G5"; gene_name "MT-X"; '
    'gene_biotype "protein_coding";',
]) + "\n"


def test_assign_genes_rank_overlap_ties_and_closed_coordinates(tmp_path):
    f = tmp_path / "t.gtf"
    f.write_text(GTF)
    g = annotate.load_genes(str(f))
    q = pd.DataFrame({
        "CHROM": ["chr1", "chr1", "chr1", "chr1", "chr1", "chrM", "chr2"],
        "START": [160, 100, 201, 1050, 400, 50, 150],
        "END": [170, 100, 210, 1060, 500, 60, 160]
    })
    out = annotate.assign_genes(q, g)
    # protein_coding beats lncRNA; 1-bp overlap at the gene start counts
    # (closed); 201-210 hits PC only; equal rank/overlap/length -> (start,
    # gene_id) order, G3 before G4; END touching gene end counts; chrM -> MT
    assert list(out) == ["PC", "LNC", "PC", "A", "PC", "MT-X", None]


def test_bsd_sum(tmp_path):
    p = tmp_path / "f"
    p.write_bytes(b"abc")
    assert annotate.bsd_sum(str(p)) == (16556, 1)  # `sum` of "abc" is 16556 1


def test_multiplicity_counts_query_lines_like_tabix():
    q = pd.DataFrame({
        "CHROM": "chr1",
        "S": [100, 100, 150, 400],
        "E": [200, 200, 160, 400]
    })  # last line empty
    cnt = call_dmr._multiplicity(pd.Series(["chr1"] * 5),
                                 np.array([99, 100, 155, 199, 200]), q)
    assert list(cnt) == [0, 2, 3, 2, 0]


def _icr_dmrs(region_rows):
    """Returns the ICR clusters of six CpGs inside the given panel rows."""
    sites = _sites(0, [1000, 1010, 1020, 1030, 1040, 1050], 0.5)
    regions = pd.DataFrame(region_rows,
                           columns=["CHROM", "REGION_START", "REGION_END"])
    found = call_dmr.find_sites_in_regions(sites, regions)
    return call_dmr.merge_sites(found, 150, by_region=True)


def test_identical_icr_rows_form_one_group_with_double_counts():
    one = _icr_dmrs([("chr1", 990, 1060)])
    two = _icr_dmrs([("chr1", 990, 1060), ("chr1", 990, 1060)])
    assert list(one["NUM_CPG"]) == [6]
    # one group, each CpG once per row
    assert list(two["NUM_CPG"]) == [12] and len(two) == 1


def test_overlapping_icr_rows_with_other_coordinates_repeat_the_dmr_row():
    # a gene row and an antisense row: different coordinates, same CpGs
    two = _icr_dmrs([("chr1", 990, 1060), ("chr1", 900, 1100)])
    assert list(two["NUM_CPG"]) == [6, 6]
    assert len(two.drop_duplicates(["CHROM", "START", "END"])) == 1


def test_the_cpg_at_end_extend_is_not_read():
    # the query is [START_EXTEND, END_EXTEND) like a BED line for tabix -R
    q = pd.DataFrame({"CHROM": ["chr1"], "S": [100], "E": [200]})
    lo, hi = call_dmr._union_intervals(q)["chr1"]
    assert (lo[0], hi[0]) == (100, 199)
    start = np.array([99, 100, 199, 200])
    cnt = call_dmr._multiplicity(pd.Series(["chr1"] * 4), start, q)
    assert list(cnt) == [0, 1, 1, 0]


def test_annotated_table_keeps_the_digits_of_the_raw_table(tmp_path):
    # The default pandas parser reads 0.000535864718615665 as
    # 0.0005358647186156 and drops two digits of the annotated P_FDR.
    raw = pd.DataFrame([{
        "CHROM": "chr1", "START": 150, "END": 160, "PANEL": "CpG",
        "PANEL_SCOPE": "autosomal", "NUM_CPG": 9, "NUM_READS": 20.0,
        "MEAN_SAMPLE": 0.9, "MEAN_REF": 0.4, "MEAN_DIFF": 0.5,
        "P_VALUE": 0.000111476322925955, "P_FDR": 0.000535864718615665,
        "DIRECTION": "hyper"
    }])[call_dmr.RAW_COLS]
    p = tmp_path / "raw.tsv.gz"
    call_dmr.write_raw(raw, str(p))
    f = tmp_path / "t.gtf"
    f.write_text(GTF)
    tab = annotate.annotate([str(p)], str(f), None, None)
    out = tmp_path / "annotated.tsv"
    annotate.write_annotated(tab, str(out))
    row = out.read_text().splitlines()[1].split("\t")
    assert row[-3:] == ["0.000111476322925955", "0.000535864718615665", "hyper"]
