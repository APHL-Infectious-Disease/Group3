"""
Legionella Pipeline Dashboard
Multi-sample support with per-sample caching, CheckM integration.
"""

import argparse
import sys
import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from io import StringIO
import json
import base64
import hashlib
from pathlib import Path
from bs4 import BeautifulSoup

# ======================================================
# CLI ARGUMENT PARSING
# ======================================================

def parse_cli_args():
    """
    Nextflow passes a directory of per-sample results, or individual file paths.
    Support both --results_dir (preferred for multi-sample) and legacy single-file flags.
    """
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--results_dir", default=None,
                        help="Root results directory from Nextflow (ont/<run_id>/)")
    # Legacy single-sample flags kept for backwards compat
    parser.add_argument("--elgato",    default=None)
    parser.add_argument("--amrfinder", default=None)
    parser.add_argument("--checkm",    default=None)
    parser.add_argument("--nanoplot",  default=None)
    parser.add_argument("--consensus", default=None)

    try:
        idx = sys.argv.index("--")
        arg_list = sys.argv[idx + 1:]
    except ValueError:
        arg_list = sys.argv[1:]

    args, _ = parser.parse_known_args(arg_list)
    return args


CLI_ARGS = parse_cli_args()


def discover_samples(results_dir: Path) -> dict[str, dict[str, Path | None]]:
    """
    Walk results_dir looking for per-sample sub-directories.
    Expected layout:
        results_dir/
          <sample_id>/
            nanoplot/nanoplot/NanoPlot-report.html
            checkm2/checkm2_out/quality_report.tsv
            medaka/consensus.fasta
            el_gato/report.json
            amrfinder/<sample_id>_amr.tsv
    Returns {sample_id: {slot: Path}}
    """
    samples = {}
    if not results_dir.is_dir():
        return samples

    for sample_dir in sorted(results_dir.iterdir()):
        if not sample_dir.is_dir():
            continue
        sid = sample_dir.name
        slot_map = {
            "nanoplot":  _find(sample_dir, "nanoplot/NanoPlot-report.html"),
            "checkm":    _find(sample_dir, "checkm/checkm.tsv"),
            "assembly":  _find(sample_dir, "medaka/consensus.fasta"),
            "elgato":    _find(sample_dir, "el_gato/report.json/report.json"),
            "amrfinder": _find_glob(sample_dir, "amrfinder/*_amr.tsv"),
        }
        samples[sid] = slot_map
    return samples


def _find(base: Path, rel: str) -> Path | None:
    p = base / rel
    return p if p.exists() else None


def _find_glob(base: Path, pattern: str) -> Path | None:
    hits = list(base.glob(pattern))
    return hits[0] if hits else None


def build_legacy_sample(args) -> dict | None:
    """Build a single fake sample from legacy --elgato / --checkm2 etc. flags."""
    slot_map = {}
    for slot, attr, key in [
        ("nanoplot",  "nanoplot",  "nanoplot"),
        ("checkm",    "checkm",    "checkm"),
        ("assembly",  "consensus", "assembly"),
        ("elgato",    "elgato",    "elgato"),
        ("amrfinder", "amrfinder", "amrfinder"),
    ]:
        val = getattr(args, attr, None)
        if val:
            p = Path(val)
            slot_map[slot] = p if p.exists() else None
        else:
            slot_map[slot] = None

    if any(v for v in slot_map.values()):
        return {"legacy_sample": slot_map}
    return None


# ======================================================
# PAGE CONFIG
# ======================================================

st.set_page_config(
    page_title="Legionella Pipeline Dashboard",
    page_icon=None,
    layout="wide",
    initial_sidebar_state="expanded",
)

st.markdown("""
<style>
    @import url('https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;600&family=IBM+Plex+Sans:wght@300;400;600&display=swap');

    html, body, [class*="css"] {
        font-family: 'IBM Plex Sans', sans-serif;
    }
    .main-header {
        text-align: center;
        background: linear-gradient(135deg, #0a0f1e 0%, #0d2137 50%, #0a1628 100%);
        padding: 2.5rem 2rem;
        border-radius: 12px;
        margin-bottom: 1.5rem;
        border: 1px solid #1e3a5f;
        position: relative;
        overflow: hidden;
    }
    .main-header::before {
        content: "";
        position: absolute; inset: 0;
        background: repeating-linear-gradient(
            90deg,
            transparent, transparent 39px,
            rgba(30,90,140,0.12) 39px, rgba(30,90,140,0.12) 40px
        ),
        repeating-linear-gradient(
            0deg,
            transparent, transparent 39px,
            rgba(30,90,140,0.12) 39px, rgba(30,90,140,0.12) 40px
        );
    }
    .main-header h1 {
        color: #e2f0ff;
        font-family: 'IBM Plex Mono', monospace;
        font-size: 1.8rem;
        font-weight: 600;
        letter-spacing: 0.05em;
        margin: 0;
        position: relative;
    }
    .main-header p {
        color: #7ab3d4;
        font-size: 0.85rem;
        letter-spacing: 0.12em;
        margin-top: 0.5rem;
        position: relative;
        font-family: 'IBM Plex Mono', monospace;
    }
    .sample-bar {
        background: #0d1f33;
        border: 1px solid #1e3a5f;
        border-radius: 8px;
        padding: 1rem 1.25rem;
        margin-bottom: 1.5rem;
        display: flex;
        align-items: center;
        gap: 1rem;
    }
    .stat-card {
        background: #0d1f33;
        border: 1px solid #1e3a5f;
        padding: 1.1rem 1rem;
        border-radius: 8px;
        text-align: center;
    }
    .metric-value {
        font-size: 1.6rem;
        font-weight: 600;
        color: #60b4f0;
        font-family: 'IBM Plex Mono', monospace;
        line-height: 1.1;
    }
    .metric-label {
        color: #6a8fa8;
        font-size: 0.75rem;
        margin-top: 0.3rem;
        text-transform: uppercase;
        letter-spacing: 0.08em;
    }
    .quality-good    { color: #34d399 !important; }
    .quality-warning { color: #fbbf24 !important; }
    .quality-poor    { color: #f87171 !important; }
    .elgato-summary {
        background: linear-gradient(135deg, #052e1c, #065f46);
        border: 1px solid #059669;
        padding: 1.5rem;
        border-radius: 10px;
        margin: 1rem 0;
        text-align: center;
    }
    .elgato-st {
        font-size: 2.2rem;
        font-weight: 600;
        color: #6ee7b7;
        font-family: 'IBM Plex Mono', monospace;
    }
    .file-status {
        display: inline-block;
        padding: 0.2rem 0.65rem;
        border-radius: 12px;
        font-size: 0.75rem;
        font-weight: 600;
        margin: 0.2rem;
        font-family: 'IBM Plex Mono', monospace;
    }
    .file-uploaded   { background: #064e3b; color: #34d399; border: 1px solid #059669; }
    .file-missing    { background: #3b0a0a; color: #f87171; border: 1px solid #dc2626; }
    .file-autoloaded { background: #0c2a4a; color: #60b4f0; border: 1px solid #2563eb; }
    .autoload-banner {
        background: #0c2a4a;
        border: 1px solid #2563eb;
        color: #93c5fd;
        padding: 0.6rem 1rem;
        border-radius: 6px;
        margin-bottom: 1rem;
        font-size: 0.85rem;
        font-family: 'IBM Plex Mono', monospace;
    }
    .checkm2-grid {
        display: grid;
        grid-template-columns: repeat(auto-fill, minmax(160px, 1fr));
        gap: 0.75rem;
        margin: 1rem 0;
    }
</style>
""", unsafe_allow_html=True)


# ======================================================
# PARSING HELPERS
# ======================================================

@st.cache_data
def parse_nanoplot_html(content: str) -> dict:
    soup = BeautifulSoup(content, "html.parser")
    stats = {}
    for table in soup.find_all("table"):
        for row in table.find_all("tr"):
            cells = row.find_all(["td", "th"])
            if len(cells) >= 2:
                k = cells[0].get_text().strip()
                v = cells[1].get_text().strip()
                if k and v and k != "General summary":
                    stats[k] = v
    return {"stats": stats, "html_content": content}


@st.cache_data
def parse_assembly_fasta(content: str) -> dict:
    sequences = []
    header, seq = "", ""
    for line in content.splitlines():
        if line.startswith(">"):
            if header:
                sequences.append({"header": header, "length": len(seq)})
            header, seq = line[1:], ""
        else:
            seq += line.strip()
    if header:
        sequences.append({"header": header, "length": len(seq)})

    lengths = [s["length"] for s in sequences]
    total   = sum(lengths)
    lengths_sorted = sorted(lengths, reverse=True)
    n50, cum = 0, 0
    for l in lengths_sorted:
        cum += l
        if cum >= total / 2:
            n50 = l
            break

    return {
        "sequences": sequences,
        "stats": {
            "total_contigs":   len(sequences),
            "total_length":    total,
            "longest_contig":  max(lengths) if lengths else 0,
            "shortest_contig": min(lengths) if lengths else 0,
            "mean_length":     float(np.mean(lengths)) if lengths else 0.0,
            "n50":             n50,
        },
    }


@st.cache_data
def parse_checkm_tsv(content: str) -> dict:
    """
    Parse CheckM TSV (-o 2 format) with dashed separator lines.
    Sample format has header/separator lines with dashes, then data row.
    """
    lines = [l.strip() for l in content.splitlines() if l.strip()]
    data = {}

    # Find the data line (not starting with dashes, not header)
    data_line = None
    for line in lines:
        if not line.startswith('-') and 'Bin Id' not in line and line:
            data_line = line
            break
    
    if not data_line:
        return {}

    # Split on multiple spaces to handle the formatted table structure
    import re
    parts = re.split(r'\s{2,}', data_line)
    
    if len(parts) < 6:
        return {}
    
    try:
        # Map the CheckM -o 2 output columns to our data structure
        # Format: Bin Id, Marker lineage, # genomes, # markers, # marker sets, 
        #         Completeness, Contamination, Strain heterogeneity, Genome size (bp), etc.
        
        data["bin_id"] = parts[0]
        data["lineage"] = parts[1].replace('(', '').replace(')', '') if len(parts) > 1 else ""
        data["completeness"] = float(parts[5]) if len(parts) > 5 else None
        data["contamination"] = float(parts[6]) if len(parts) > 6 else None  
        data["strain_heterogeneity"] = float(parts[7]) if len(parts) > 7 else None
        
        # Assembly stats from CheckM extended output
        if len(parts) > 8:
            data["genome_size_bp"] = int(parts[8]) if parts[8].isdigit() else None
        if len(parts) > 10:
            data["num_contigs"] = int(parts[10]) if parts[10].isdigit() else None
        if len(parts) > 12:
            data["n50_contigs"] = int(parts[12]) if parts[12].isdigit() else None
        if len(parts) > 14:
            data["mean_contig_len"] = int(float(parts[14])) if parts[14].replace('.','').isdigit() else None
        if len(parts) > 16:
            data["longest_contig"] = int(parts[16]) if parts[16].isdigit() else None
            
    except (ValueError, IndexError) as e:
        # Fallback: try to extract basic completeness/contamination with regex
        import re
        comp_match = re.search(r'(\d+\.\d+).*?(\d+\.\d+)', data_line)
        if comp_match:
            try:
                data["completeness"] = float(comp_match.group(1))
                data["contamination"] = float(comp_match.group(2))
            except ValueError:
                pass

    return data


@st.cache_data
def parse_elgato_json(content: str) -> dict:
    return json.loads(content)


@st.cache_data
def parse_amrfinder_tsv(content: str):
    return pd.read_csv(StringIO(content), sep="\t")


def get_file_hash(f) -> str:
    if f is None:
        return ""
    if hasattr(f, "getvalue"):
        return hashlib.md5(f.getvalue()).hexdigest()
    if isinstance(f, Path):
        return hashlib.md5(f.read_bytes()).hexdigest()
    return ""


def read_file(f) -> str | None:
    if f is None:
        return None
    if isinstance(f, Path):
        return f.read_text(encoding="utf-8", errors="replace")
    try:
        return f.read().decode("utf-8")
    except Exception:
        return None


# ======================================================
# PER-SAMPLE CACHE  (stored in session_state)
# ======================================================

def ensure_sample_cache(sample_id: str):
    if "sample_cache" not in st.session_state:
        st.session_state.sample_cache = {}
    if sample_id not in st.session_state.sample_cache:
        st.session_state.sample_cache[sample_id] = {
            "parsed": {},
            "file_hashes": {},
        }


def get_sample_cache(sample_id: str) -> dict:
    ensure_sample_cache(sample_id)
    return st.session_state.sample_cache[sample_id]


def parse_and_cache(sample_id: str, slot: str, source) -> None:
    """Parse a file and store into per-sample cache."""
    cache = get_sample_cache(sample_id)
    content = read_file(source)
    if content is None:
        return
    try:
        if slot == "nanoplot":
            cache["parsed"]["nanoplot"]  = parse_nanoplot_html(content)
        elif slot == "assembly":
            cache["parsed"]["assembly"]  = parse_assembly_fasta(content)
        elif slot == "elgato":
            cache["parsed"]["elgato"]    = parse_elgato_json(content)
        elif slot == "checkm":
            cache["parsed"]["checkm"]    = parse_checkm_tsv(content)
        elif slot == "amrfinder":
            cache["parsed"]["amrfinder"] = parse_amrfinder_tsv(content)
    except Exception as e:
        st.warning(f"[{sample_id}] Could not parse {slot}: {e}")


def load_sample_files(sample_id: str, slot_map: dict[str, Path | None]) -> None:
    """Load all file slots for a sample into cache (idempotent via hash check)."""
    cache = get_sample_cache(sample_id)
    for slot, path in slot_map.items():
        if path is None:
            continue
        h = get_file_hash(path)
        if h != cache["file_hashes"].get(slot):
            cache["file_hashes"][slot] = h
            parse_and_cache(sample_id, slot, path)


# ======================================================
# MAIN DASHBOARD
# ======================================================

class LegionellaPipelineDashboard:

    def __init__(self):
        if "all_samples" not in st.session_state:
            st.session_state.all_samples = {}    # {sample_id: slot_map}
        if "active_sample" not in st.session_state:
            st.session_state.active_sample = None
        if "manual_slots" not in st.session_state:
            st.session_state.manual_slots = {}   # slot -> UploadedFile (manual upload)
        if "manual_hashes" not in st.session_state:
            st.session_state.manual_hashes = {}

    # --------------------------------------------------
    # AUTO-DISCOVER SAMPLES FROM CLI / results_dir
    # --------------------------------------------------

    def auto_discover(self):
        if CLI_ARGS.results_dir:
            rd = Path(CLI_ARGS.results_dir)
            samples = discover_samples(rd)
            if samples:
                st.session_state.all_samples.update(samples)
                if st.session_state.active_sample is None:
                    st.session_state.active_sample = next(iter(samples))
                return

        legacy = build_legacy_sample(CLI_ARGS)
        if legacy:
            st.session_state.all_samples.update(legacy)
            if st.session_state.active_sample is None:
                st.session_state.active_sample = "legacy_sample"

    # --------------------------------------------------
    # SIDEBAR: sample selector + manual upload
    # --------------------------------------------------

    def render_sidebar(self):
        with st.sidebar:
            st.markdown("## Legionella Dashboard")
            st.divider()

            # ---- Sample selector ----
            all_sids = list(st.session_state.all_samples.keys())

            if all_sids:
                st.markdown("### Sample")
                chosen = st.selectbox(
                    "Select sample",
                    options=all_sids,
                    index=all_sids.index(st.session_state.active_sample)
                          if st.session_state.active_sample in all_sids else 0,
                    key="sample_selectbox",
                    label_visibility="collapsed",
                )
                st.session_state.active_sample = chosen
                st.markdown(f"**{len(all_sids)}** sample(s) loaded")
                st.divider()

            # ---- Manual upload section ----
            with st.expander("Upload Files", expanded=not all_sids):
                manual_sid = st.text_input("Sample ID (for uploaded files)",
                                           value="manual_sample")

                nano = st.file_uploader("NanoPlot HTML",  type=["html"])
                chkm = st.file_uploader("CheckM TSV",     type=["tsv", "txt"])
                asm  = st.file_uploader("Consensus FASTA", type=["fasta","fa","fna"])
                elg  = st.file_uploader("El Gato JSON",   type=["json"])
                amr  = st.file_uploader("AMRFinder TSV",  type=["tsv","txt"])

                uploads = {
                    "nanoplot":  nano,
                    "checkm":    chkm,
                    "assembly":  asm,
                    "elgato":    elg,
                    "amrfinder": amr,
                }

                any_new = False
                for slot, uf in uploads.items():
                    if uf is not None:
                        h = get_file_hash(uf)
                        key = f"{manual_sid}_{slot}"
                        if h != st.session_state.manual_hashes.get(key):
                            st.session_state.manual_hashes[key] = h
                            if manual_sid not in st.session_state.all_samples:
                                st.session_state.all_samples[manual_sid] = {}
                            st.session_state.all_samples[manual_sid][slot] = uf
                            uf.seek(0)
                            parse_and_cache(manual_sid, slot, uf)
                            any_new = True

                if any_new and st.session_state.active_sample is None:
                    st.session_state.active_sample = manual_sid

            # ---- Status badges for active sample ----
            if st.session_state.active_sample:
                st.divider()
                st.markdown("### File Status")
                sid = st.session_state.active_sample
                slot_map = st.session_state.all_samples.get(sid, {})
                cache = get_sample_cache(sid)
                labels = {
                    "nanoplot":  "NanoPlot",
                    "checkm":    "CheckM",
                    "assembly":  "Consensus",
                    "elgato":    "El Gato",
                    "amrfinder": "AMRFinder",
                }
                for slot, label in labels.items():
                    src = slot_map.get(slot)
                    parsed = cache["parsed"].get(slot)
                    if parsed is not None:
                        cls, icon = "file-uploaded", "+"
                    elif src is not None:
                        cls, icon = "file-autoloaded", "~"
                    else:
                        cls, icon = "file-missing", "-"
                    st.markdown(
                        f'<span class="file-status {cls}">{icon} {label}</span>',
                        unsafe_allow_html=True,
                    )

    # --------------------------------------------------
    # HEADER
    # --------------------------------------------------

    def render_header(self):
        sid = st.session_state.active_sample or ""
        st.markdown(f"""
        <div class="main-header">
            <h1>Legionella Pipeline Dashboard</h1>
            <p>NANOPLOT QC &nbsp;·&nbsp; RAVEN-MEDAKA ASSEMBLY &nbsp;·&nbsp;
               CHECKM2 &nbsp;·&nbsp; AMRFINDER &nbsp;·&nbsp; EL GATO TYPING</p>
        </div>
        """, unsafe_allow_html=True)

    # --------------------------------------------------
    # LOAD ACTIVE SAMPLE
    # --------------------------------------------------

    def load_active_sample(self) -> dict:
        sid = st.session_state.active_sample
        if not sid:
            return {}
        slot_map = st.session_state.all_samples.get(sid, {})
        # Ensure file Paths are loaded into cache
        path_slots = {k: v for k, v in slot_map.items() if isinstance(v, Path)}
        if path_slots:
            load_sample_files(sid, path_slots)
        return get_sample_cache(sid)["parsed"]

    # --------------------------------------------------
    # TABS
    # --------------------------------------------------

    def render_tabs(self, parsed: dict):
        sid = st.session_state.active_sample

        # Sample pill
        if sid:
            st.markdown(
                f'<div class="autoload-banner">Viewing sample: <strong>{sid}</strong></div>',
                unsafe_allow_html=True,
            )

        tab1, tab2, tab3 = st.tabs(["Read QC", "Assembly QC", "Results & Annotations"])

        with tab1:
            self._render_read_qc(parsed)
        with tab2:
            self._render_assembly_qc(parsed)
        with tab3:
            self._render_results(parsed)

    # ---- Tab 1: Read QC ----------------------------------------

    def _render_read_qc(self, parsed: dict):
        st.markdown("## Read Quality Control")
        nano_data = parsed.get("nanoplot")

        if nano_data is None:
            st.info("No NanoPlot report loaded for this sample.")
            return

        stats = nano_data.get("stats", {})
        if stats:
            st.markdown("### NanoPlot Summary Statistics")
            cols = st.columns(min(4, len(stats)))
            for i, (k, v) in enumerate(stats.items()):
                cols[i % 4].metric(k, v)
        else:
            st.warning("No summary statistics found in the NanoPlot report.")

        if st.checkbox("Show embedded NanoPlot report"):
            html = nano_data.get("html_content", "")
            b64  = base64.b64encode(html.encode()).decode()
            st.markdown(
                f'<iframe src="data:text/html;base64,{b64}" '
                f'width="100%" height="800" frameborder="0"></iframe>',
                unsafe_allow_html=True,
            )

    # ---- Tab 2: Assembly QC ------------------------------------

    def _render_assembly_qc(self, parsed: dict):
        st.markdown("## Assembly Quality Control")

        checkm = parsed.get("checkm")
        if checkm:
            st.markdown("### CheckM Assembly Assessment")
            self._render_checkm_stats(checkm)
        else:
            st.info("No CheckM report loaded for this sample.")

        st.divider()

        asm = parsed.get("assembly")
        if asm:
            st.markdown("### Contig Statistics (from FASTA)")
            s = asm["stats"]
            cols = st.columns(3)
            pairs = [
                ("Total Contigs",     f'{s["total_contigs"]:,}'),
                ("Total Length (bp)", f'{s["total_length"]:,}'),
                ("N50 (bp)",          f'{s["n50"]:,}'),
                ("Longest (bp)",      f'{s["longest_contig"]:,}'),
                ("Shortest (bp)",     f'{s["shortest_contig"]:,}'),
                ("Mean Length (bp)",  f'{s["mean_length"]:,.0f}'),
            ]
            for idx, (label, val) in enumerate(pairs):
                cols[idx % 3].markdown(
                    f'<div class="stat-card">'
                    f'<div class="metric-value">{val}</div>'
                    f'<div class="metric-label">{label}</div></div>',
                    unsafe_allow_html=True,
                )
            lengths = [s["length"] for s in asm["sequences"]]
            fig = go.Figure(go.Histogram(x=lengths, nbinsx=30,
                                         marker_color="#2563eb", opacity=0.8))
            fig.update_layout(
                title="Contig Length Distribution",
                xaxis_title="Length (bp)", yaxis_title="Count", height=350,
                paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(13,31,51,0.8)",
                font=dict(color="#93c5fd"),
            )
            st.plotly_chart(fig, use_container_width=True)

    def _render_checkm_stats(self, d: dict):
        """Render CheckM completeness, contamination, and assembly stats."""
        comp = d.get("completeness")
        cont = d.get("contamination")
        sh = d.get("strain_heterogeneity")

        if comp is None:
            st.warning("CheckM QC values could not be parsed from the TSV.")
            st.json(d)
            return

        # Quality verdict
        if comp >= 95 and cont <= 5:
            st.success("High Quality — completeness >= 95%, contamination <= 5%")
        elif comp >= 90 and cont <= 10:
            st.warning("Medium Quality — review assembly parameters")
        else:
            st.error("Low Quality — consider reassembly or additional filtering")

        # Core QC metrics
        c1, c2, c3 = st.columns(3)
        comp_cls = "quality-good" if comp >= 95 else "quality-warning" if comp >= 90 else "quality-poor"
        cont_cls = "quality-good" if cont <= 5  else "quality-warning" if cont <= 10  else "quality-poor"

        for col, label, val, cls in [
            (c1, "Completeness",        f"{comp:.2f}%", comp_cls),
            (c2, "Contamination",       f"{cont:.2f}%", cont_cls),
            (c3, "Strain Heterogeneity",f"{sh:.2f}%" if sh is not None else "—", ""),
        ]:
            col.markdown(
                f'<div class="stat-card">'
                f'<div class="metric-value {cls}">{val}</div>'
                f'<div class="metric-label">{label}</div></div>',
                unsafe_allow_html=True,
            )

        # Assembly-level stats from CheckM extended output (-o 2)
        assembly_fields = {
            "Genome Size (bp)":      d.get("genome_size_bp"),
            "# Contigs":             d.get("num_contigs"),
            "N50 Contigs (bp)":      d.get("n50_contigs"),
            "Mean Contig Len (bp)":  d.get("mean_contig_len"),
            "Longest Contig (bp)":   d.get("longest_contig"),
        }
        available = {k: v for k, v in assembly_fields.items() if v is not None}

        if available:
            st.markdown("#### Assembly Statistics (from CheckM)")
            cols = st.columns(min(3, len(available)))
            for i, (label, val) in enumerate(available.items()):
                cols[i % 3].markdown(
                    f'<div class="stat-card">'
                    f'<div class="metric-value">{val:,}</div>'
                    f'<div class="metric-label">{label}</div></div>',
                    unsafe_allow_html=True,
                )

        # Lineage / bin info
        if d.get("bin_id") or d.get("lineage"):
            st.caption(f"Bin: {d.get('bin_id','')}   |   Lineage: {d.get('lineage','')}")

    # ---- Tab 3: Results & Annotations --------------------------

    def _render_results(self, parsed: dict):
        st.markdown("## Results & Annotations")

        elgato = parsed.get("elgato")
        if elgato:
            st.markdown("### El Gato Sequence Typing")
            self._render_elgato(elgato)
        else:
            st.info("No El Gato report loaded for this sample.")

        st.divider()

        amr = parsed.get("amrfinder")
        if amr is not None:
            st.markdown("### AMR Gene Annotations")
            self._render_amrfinder(amr)
        else:
            st.info("No AMRFinder results loaded for this sample.")

    def _render_elgato(self, data: dict):
        try:
            mlst = data["mlst"]
            hits = data["mode_specific"]["BLAST_hit_locations"]
        except KeyError as e:
            st.error(f"Unexpected El Gato JSON structure — missing key: {e}")
            st.json(data)
            return

        st.markdown(f"""
        <div class="elgato-summary">
            <div class="elgato-st">ST: {mlst.get("st", "ND")}</div>
            <p style="color:#6ee7b7;margin:0.3rem 0 0;">Sample: {data.get("id","")}</p>
        </div>
        """, unsafe_allow_html=True)

        allele_cols = ["flaA","pilE","asd","mip","mompS","proA","neuA_neuAH"]
        mlst_df = pd.DataFrame([{
            "sample_id": data.get("id",""),
            "ST": mlst.get("st","ND"),
            **{c: mlst.get(c,"ND") for c in allele_cols}
        }])

        rows = []
        for locus, lhits in hits.items():
            for hit in lhits:
                rows.append({
                    "locus":  locus,
                    "allele": hit[0],
                    "contig": hit[1],
                    "start":  int(hit[2]),
                    "end":    int(hit[3]),
                    "length": int(hit[4]),
                })
        hits_df = pd.DataFrame(rows).sort_values("start") if rows else pd.DataFrame()

        c1, c2 = st.columns([1, 2])
        with c1:
            st.subheader("MLST Profile")
            st.table(mlst_df.set_index("sample_id").T)
        with c2:
            st.subheader("BLAST Hits")
            if not hits_df.empty:
                st.dataframe(
                    hits_df.style.format({"start":"{:,}","end":"{:,}","length":"{:,}"}),
                    use_container_width=True,
                )
            else:
                st.info("No BLAST hits found.")

        if st.checkbox("Show raw El Gato JSON"):
            st.json(data)

    def _render_amrfinder(self, df: pd.DataFrame):
        if df.empty:
            st.success("No AMR genes detected.")
            return

        c1, c2, c3 = st.columns(3)
        c1.metric("Total AMR Genes", len(df))
        if "Class" in df.columns:
            c2.metric("Drug Classes", df["Class"].nunique())
            top = df["Class"].value_counts().index[0] if not df["Class"].value_counts().empty else "—"
            c3.metric("Top Class", top)

        display_cols = [
            "Gene symbol", "Sequence name", "Class", "Subclass",
            "% Coverage of reference sequence", "% Identity to reference sequence",
        ]
        show = [c for c in display_cols if c in df.columns]
        st.dataframe(df[show].round(2), use_container_width=True)

        if "Class" in df.columns:
            counts = df["Class"].value_counts()
            fig = go.Figure(go.Bar(
                x=counts.index, y=counts.values, marker_color="#f87171"
            ))
            fig.update_layout(
                title="AMR Class Distribution",
                xaxis_title="Class", yaxis_title="Count", height=350,
                paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(13,31,51,0.8)",
                font=dict(color="#93c5fd"),
            )
            st.plotly_chart(fig, use_container_width=True)

    # --------------------------------------------------
    # ENTRYPOINT
    # --------------------------------------------------

    def run(self):
        self.auto_discover()
        self.render_sidebar()
        self.render_header()

        if not st.session_state.all_samples:
            st.info("Upload pipeline output files in the sidebar to get started.")
            return

        if not st.session_state.active_sample:
            st.info("Select a sample from the sidebar.")
            return

        parsed = self.load_active_sample()
        self.render_tabs(parsed)


# ======================================================
# MAIN
# ======================================================

if __name__ == "__main__":
    LegionellaPipelineDashboard().run()
