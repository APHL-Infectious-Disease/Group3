"""
Legionella Pipeline Dashboard
A Streamlit dashboard for collating QC, assembly, SBT typing, and more
from a long-read Legionella pneumophila sequencing pipeline.

el gato: https://github.com/CDCgov/el_gato

Can be launched manually (file upload UI) or automatically by the Nextflow
pipeline with pre-loaded file paths passed as CLI arguments:

    streamlit run dashboard.py -- \
        --elgato    results/ont/el_gato/report.json \
        --amrfinder results/ont/amrfinder/sample_amr.tsv \
        --checkm    results/ont/checkm/checkm.tsv \
        --nanoplot  results/ont/nanoplot/nanoplot/NanoPlot-report.html \
        --consensus results/ont/medaka/consensus.fasta

All arguments are optional — the app functions with any subset of files.
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
import graphviz

# ======================================================
# CLI ARGUMENT PARSING
# Streamlit passes extra args after '--' to sys.argv.
# We parse them here before any Streamlit calls.
# ======================================================

def parse_cli_args():
    """Parse file paths — works whether launched by Nextflow or manually."""
    import os
    
    # Streamlit passes post-'--' args but may mangle sys.argv
    # Use st.query_params as fallback, or read directly from the raw argv
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--elgato",    default=None)
    parser.add_argument("--amrfinder", default=None)
    parser.add_argument("--checkm",    default=None)
    parser.add_argument("--nanoplot",  default=None)
    parser.add_argument("--consensus", default=None)

    # Try everything after '--' first, then fall back to all args
    try:
        idx = sys.argv.index("--")
        arg_list = sys.argv[idx + 1:]
    except ValueError:
        # '--' was stripped by Streamlit — try parsing all argv
        arg_list = sys.argv[1:]

    args, _ = parser.parse_known_args(arg_list)

    result = {}
    for key in ("elgato", "amrfinder", "checkm", "nanoplot", "consensus"):
        val = getattr(args, key)
        if val and Path(val).exists():
            result[key] = Path(val)
        else:
            result[key] = None

    return result


CLI_PATHS = parse_cli_args()

# ======================================================
# PAGE CONFIG
# ======================================================

st.set_page_config(
    page_title="Legionella Pipeline Dashboard",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

st.markdown("""
<style>
    .main-header {
        text-align: center;
        background: linear-gradient(135deg, #1e3a8a, #3b82f6);
        padding: 2rem;
        border-radius: 15px;
        margin-bottom: 2rem;
        color: white;
    }
    .upload-section {
        background: linear-gradient(135deg, #f8fafc, #e2e8f0);
        padding: 1.5rem;
        border-radius: 10px;
        margin: 1rem 0;
        border: 2px dashed #64748b;
    }
    .file-status {
        display: inline-block;
        padding: 0.25rem 0.75rem;
        border-radius: 15px;
        font-size: 0.8rem;
        font-weight: 600;
        margin: 0.2rem;
    }
    .file-uploaded  { background-color: #10b981; color: white; }
    .file-missing   { background-color: #ef4444; color: white; }
    .file-autoloaded { background-color: #3b82f6; color: white; }
    .module-card {
        background: linear-gradient(135deg, #0f172a, #1e293b);
        color: white;
        padding: 1.5rem;
        border-radius: 12px;
        margin: 1rem 0;
        border: 1px solid #475569;
    }
    .stat-card {
        background: white;
        padding: 1rem;
        border-radius: 8px;
        box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        text-align: center;
        margin: 0.5rem 0;
    }
    .metric-value { font-size: 2rem; font-weight: bold; color: #1e40af; }
    .metric-label { color: #64748b; font-size: 0.9rem; }
    .quality-good    { color: #059669; font-weight: bold; }
    .quality-warning { color: #d97706; font-weight: bold; }
    .quality-poor    { color: #dc2626; font-weight: bold; }
    .elgato-summary {
        background: linear-gradient(135deg, #065f46, #047857);
        color: white;
        padding: 1.5rem;
        border-radius: 10px;
        margin: 1rem 0;
    }
    .elgato-st { font-size: 2rem; font-weight: bold; text-align: center; }
    .autoload-banner {
        background: linear-gradient(135deg, #1e40af, #3b82f6);
        color: white;
        padding: 0.75rem 1.25rem;
        border-radius: 8px;
        margin-bottom: 1rem;
        font-size: 0.95rem;
    }
</style>
""", unsafe_allow_html=True)


# ======================================================
# PARSING HELPERS  (cached so re-uploads don't reparse)
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
    Parse CheckM TSV (-o 2 format).  Returns a dict with keys:
      completeness, contamination, strain_heterogeneity, bin_id, lineage,
      and also assembly-level fields if present (genome_size, num_contigs,
      n50_contigs, mean_contig_length, longest_contig, coding_density).
    """
    import re
    lines = [l for l in content.splitlines() if l.strip() and not l.startswith("---")]
    data = {}

    for line in lines:
        # skip the header row
        if "Completeness" in line and "Contamination" in line:
            continue
        parts = line.split()
        if not parts:
            continue

        # The -o 2 format has many columns; grab by regex for robustness
        floats = re.findall(r"\d+\.\d+", line)
        if len(floats) >= 3:
            data["completeness"]          = float(floats[0])
            data["contamination"]         = float(floats[1])
            data["strain_heterogeneity"]  = float(floats[2])

        # Assembly columns present in -o 2 extended output
        # Genome size, # contigs, N50 (contigs), mean contig length, longest contig
        ints = re.findall(r"\b(\d{4,})\b", line)   # integers ≥ 4 digits = assembly stats
        if len(ints) >= 1:
            data["genome_size_bp"]      = int(ints[0]) if len(ints) > 0 else None
            data["num_contigs"]         = int(ints[1]) if len(ints) > 1 else None
            data["n50_contigs"]         = int(ints[2]) if len(ints) > 2 else None
            data["mean_contig_len"]     = int(ints[3]) if len(ints) > 3 else None
            data["longest_contig"]      = int(ints[4]) if len(ints) > 4 else None

        # Bin ID and lineage (first two text tokens)
        text_parts = [p for p in parts if not re.fullmatch(r"[\d.]+", p)]
        if text_parts:
            data["bin_id"]  = text_parts[0]
            data["lineage"] = " ".join(text_parts[1:3]) if len(text_parts) > 1 else ""

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


def read_file(f) -> str | bytes | None:
    """Read from an UploadedFile or a Path, return decoded string."""
    if f is None:
        return None
    if isinstance(f, Path):
        return f.read_text(encoding="utf-8", errors="replace")
    try:
        return f.read().decode("utf-8")
    except Exception:
        return None


# ======================================================
# MAIN DASHBOARD CLASS
# ======================================================

class LegionellaPipelineDashboard:

    def __init__(self):
        self._init_session_state()

    # --------------------------------------------------
    # SESSION STATE
    # --------------------------------------------------

    def _init_session_state(self):
        defaults = {
            "uploaded_files": {k: None for k in
                               ("nanoplot", "assembly", "elgato",
                                "checkm", "amrfinder", "assembly_graph")},
            "parsed_data":  {},
            "file_hashes":  {},
            "auto_loaded":  False,
        }
        for k, v in defaults.items():
            if k not in st.session_state:
                st.session_state[k] = v

    # --------------------------------------------------
    # AUTO-LOAD FROM CLI PATHS (Nextflow integration)
    # --------------------------------------------------

    def auto_load_cli_files(self):
        """
        If the pipeline injected file paths via CLI args, load them into
        session state once (idempotent — checked via hash).
        """
        mapping = {
            "nanoplot":  CLI_PATHS.get("nanoplot"),
            "elgato":    CLI_PATHS.get("elgato"),
            "checkm":    CLI_PATHS.get("checkm"),
            "amrfinder": CLI_PATHS.get("amrfinder"),
            "assembly":  CLI_PATHS.get("consensus"),   # consensus.fasta → assembly slot
        }

        any_loaded = False
        for slot, path in mapping.items():
            if path is None:
                continue
            h = get_file_hash(path)
            if h != st.session_state.file_hashes.get(f"cli_{slot}"):
                st.session_state.file_hashes[f"cli_{slot}"] = h
                st.session_state.uploaded_files[slot] = path
                self._parse_file(slot, path)
                any_loaded = True

        if any_loaded:
            st.session_state.auto_loaded = True

    def _parse_file(self, slot: str, source):
        """Parse a file (Path or UploadedFile) into parsed_data."""
        content = read_file(source)
        if content is None:
            return
        try:
            if slot == "nanoplot":
                st.session_state.parsed_data["nanoplot"] = parse_nanoplot_html(content)
            elif slot == "assembly":
                st.session_state.parsed_data["assembly"] = parse_assembly_fasta(content)
            elif slot == "elgato":
                st.session_state.parsed_data["elgato"]   = parse_elgato_json(content)
            elif slot == "checkm":
                st.session_state.parsed_data["checkm"]   = parse_checkm_tsv(content)
            elif slot == "amrfinder":
                st.session_state.parsed_data["amrfinder"] = parse_amrfinder_tsv(content)
            elif slot == "assembly_graph":
                st.session_state.parsed_data["assembly_graph"] = content
        except Exception as e:
            st.warning(f"Could not parse {slot}: {e}")

    # --------------------------------------------------
    # HEADER
    # --------------------------------------------------

    def render_header(self):
        st.markdown("""
        <div class="main-header">
            <h1>🧬 Legionella Pipeline Dashboard</h1>
            <p style="font-size:1.1rem;opacity:0.9;">
                NanoPlot QC &nbsp;•&nbsp; Raven-Medaka Assembly &nbsp;•&nbsp;
                CheckM Quality &nbsp;•&nbsp; AMRFinder &nbsp;•&nbsp; El Gato Typing
            </p>
        </div>
        """, unsafe_allow_html=True)

        if st.session_state.auto_loaded:
            loaded = [k for k, v in CLI_PATHS.items() if v is not None]
            st.markdown(
                f'<div class="autoload-banner">⚡ <strong>Auto-loaded from pipeline output:</strong> '
                f'{", ".join(loaded)}</div>',
                unsafe_allow_html=True
            )

    # --------------------------------------------------
    # FILE UPLOAD SECTION  (manual upload — always shown)
    # --------------------------------------------------

    def render_file_upload_section(self):
        with st.expander("📂 Upload / Replace Files", expanded=not st.session_state.auto_loaded):
            self._render_status_badges()

            col1, col2, col3 = st.columns(3)

            with col1:
                st.markdown("**Quality Control**")
                nano = st.file_uploader("NanoPlot Report (HTML)", type=["html"], key="nanoplot_up")
                chk  = st.file_uploader("CheckM Output (TSV/TXT)", type=["tsv","txt"], key="checkm_up")

            with col2:
                st.markdown("**Assembly**")
                asm  = st.file_uploader("Consensus FASTA", type=["fasta","fa","fna"], key="assembly_up")
                gv   = st.file_uploader("Assembly Graph (GV/DOT)", type=["gv","dot"], key="graph_up")

            with col3:
                st.markdown("**Analysis Results**")
                elg  = st.file_uploader("El Gato Report (JSON)", type=["json"], key="elgato_up")
                amr  = st.file_uploader("AMRFinder Results (TSV)", type=["tsv","txt"], key="amr_up")

            # Process any newly uploaded files
            for slot, uploaded in [
                ("nanoplot", nano), ("checkm", chk), ("assembly", asm),
                ("assembly_graph", gv), ("elgato", elg), ("amrfinder", amr)
            ]:
                if uploaded is not None:
                    h = get_file_hash(uploaded)
                    if h != st.session_state.file_hashes.get(slot):
                        st.session_state.file_hashes[slot] = h
                        st.session_state.uploaded_files[slot] = uploaded
                        uploaded.seek(0)
                        self._parse_file(slot, uploaded)

    def _render_status_badges(self):
        labels = {
            "nanoplot": "NanoPlot", "assembly": "Assembly FASTA",
            "elgato": "El Gato JSON", "checkm": "CheckM TSV",
            "amrfinder": "AMRFinder TSV", "assembly_graph": "Assembly Graph",
        }
        html = ""
        for slot, label in labels.items():
            src  = st.session_state.uploaded_files.get(slot)
            auto = isinstance(src, Path)
            if src is not None and auto:
                cls, icon = "file-autoloaded", "⚡"
            elif src is not None:
                cls, icon = "file-uploaded", "✅"
            else:
                cls, icon = "file-missing", "❌"
            html += f'<span class="file-status {cls}">{icon} {label}</span>'
        st.markdown(html, unsafe_allow_html=True)

        loaded = sum(1 for v in st.session_state.uploaded_files.values() if v is not None)
        total  = len(st.session_state.uploaded_files)
        st.progress(loaded / total)
        st.caption(f"{loaded}/{total} files loaded")

    # --------------------------------------------------
    # MAIN MODULES
    # --------------------------------------------------

    def render_analysis_modules(self):
        if not any(st.session_state.uploaded_files.values()):
            st.info("Upload pipeline output files above to begin analysis.")
            return

        tab1, tab2, tab3 = st.tabs(["📊 Read QC", "🧱 Assembly QC", "🔬 Results & Annotations"])
        with tab1:
            self._render_read_qc()
        with tab2:
            self._render_assembly_qc()
        with tab3:
            self._render_results()

    # ---- Tab 1: Read QC ----------------------------------------

    def _render_read_qc(self):
        st.markdown("## Read Quality Control")
        nano_data = st.session_state.parsed_data.get("nanoplot")

        if nano_data is None:
            st.info("Upload a NanoPlot HTML report to view read quality metrics.")
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

    def _render_assembly_qc(self):
        st.markdown("## Assembly Quality Control")

        # --- CheckM section (primary assembly stats source) ------
        checkm = st.session_state.parsed_data.get("checkm")
        if checkm:
            st.markdown("### CheckM Assembly Assessment")
            self._render_checkm_stats(checkm)
        else:
            st.info("Upload a CheckM TSV to view assembly quality and completeness.")

        st.divider()

        # --- FASTA contig stats (supplementary) ------------------
        asm = st.session_state.parsed_data.get("assembly")
        if asm:
            st.markdown("### Contig Statistics (from FASTA)")
            s = asm["stats"]
            c1, c2, c3 = st.columns(3)
            pairs = [
                ("Total Contigs",     f'{s["total_contigs"]:,}'),
                ("Total Length (bp)", f'{s["total_length"]:,}'),
                ("N50 (bp)",          f'{s["n50"]:,}'),
                ("Longest (bp)",      f'{s["longest_contig"]:,}'),
                ("Shortest (bp)",     f'{s["shortest_contig"]:,}'),
                ("Mean Length (bp)",  f'{s["mean_length"]:,.0f}'),
            ]
            for idx, (label, val) in enumerate(pairs):
                col = [c1, c2, c3][idx % 3]
                col.markdown(
                    f'<div class="stat-card">'
                    f'<div class="metric-value">{val}</div>'
                    f'<div class="metric-label">{label}</div></div>',
                    unsafe_allow_html=True,
                )

            lengths = [s["length"] for s in asm["sequences"]]
            fig = go.Figure(go.Histogram(x=lengths, nbinsx=30,
                                         marker_color="#10b981", opacity=0.7))
            fig.update_layout(title="Contig Length Distribution",
                              xaxis_title="Length (bp)", yaxis_title="Count", height=350)
            st.plotly_chart(fig, use_container_width=True)

        # --- Assembly graph --------------------------------------
        gv_content = st.session_state.parsed_data.get("assembly_graph")
        if gv_content:
            st.markdown("### Assembly Graph")
            try:
                st.graphviz_chart(gv_content)
            except Exception:
                st.warning("Could not render graph — showing raw DOT source.")
                st.code(gv_content, language="dot")

    def _render_checkm_stats(self, d: dict):
        """Render CheckM completeness, contamination, and assembly stats."""
        comp  = d.get("completeness")
        cont  = d.get("contamination")
        sh    = d.get("strain_heterogeneity")

        if comp is None:
            st.warning("CheckM QC values could not be parsed from the TSV.")
            st.json(d)
            return

        # Quality verdict
        if comp >= 95 and cont <= 5:
            st.success("🟢 **High Quality** — completeness ≥ 95%, contamination ≤ 5%")
        elif comp >= 90 and cont <= 10:
            st.warning("🟡 **Medium Quality** — review assembly parameters")
        else:
            st.error("🔴 **Low Quality** — consider reassembly or additional filtering")

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

    def _render_results(self):
        st.markdown("## Results & Annotations")

        elgato = st.session_state.parsed_data.get("elgato")
        if elgato:
            st.markdown("### El Gato Sequence Typing")
            self._render_elgato(elgato)
        else:
            st.info("Upload an El Gato JSON report to view sequence typing results.")

        st.divider()

        amr = st.session_state.parsed_data.get("amrfinder")
        if amr is not None:
            st.markdown("### AMR Gene Annotations")
            self._render_amrfinder(amr)
        else:
            st.info("Upload an AMRFinder TSV to view resistance gene annotations.")

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
            <div class="elgato-st">Sequence Type: {mlst.get("st", "ND")}</div>
            <p style="text-align:center;margin-top:0.5rem;">Sample: {data.get("id","")}</p>
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
            st.success("🟢 No AMR genes detected.")
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
                x=counts.index, y=counts.values, marker_color="#dc2626"
            ))
            fig.update_layout(title="AMR Class Distribution",
                              xaxis_title="Class", yaxis_title="Count", height=350)
            st.plotly_chart(fig, use_container_width=True)

    # --------------------------------------------------
    # ENTRYPOINT
    # --------------------------------------------------

    def run(self):
        self.auto_load_cli_files()

        # Refuse to start if launched with no pipeline files at all
        if not any(CLI_PATHS.values()) and not any(st.session_state.uploaded_files.values()):
            st.error("Dashboard was launched without any pipeline output files. "
                 "This app should be started automatically by the ONT pipeline.")
            st.stop()

        self.render_header()
        self.render_file_upload_section()
        self.render_analysis_modules()


# ======================================================
# MAIN
# ======================================================

if __name__ == "__main__":
    LegionellaPipelineDashboard().run()