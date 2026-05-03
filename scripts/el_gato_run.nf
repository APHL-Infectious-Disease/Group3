#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ======================================================
// PARAMETERS
// ======================================================
def timestamp = new Date().format('yyMMddHHmmss')

params.read1       = null
params.read2       = null
params.fasta       = null
params.samplesheet = null
params.sample_id   = "sample"
params.outdir      = "results"
params.run_id      = timestamp
params.checkm2_db  = "/mnt/checkm2_db"
params.dashboard   = "${projectDir}/dashboard.py"

// ======================================================
// ILLUMINA MODE — FASTQ
// ======================================================

process EL_GATO_FASTQ {

    publishDir "${params.outdir}/illumina_fastq/${params.run_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    path "${sample_id}"

    script:
    """
    el_gato.py --read1 ${read1} --read2 ${read2} --out ${sample_id}
    """
}

// ======================================================
// ILLUMINA MODE — FASTA
// ======================================================

process EL_GATO_FASTA {

    publishDir "${params.outdir}/illumina_fasta/${params.run_id}", mode: 'copy'

    input:
    path fasta

    output:
    path "${fasta.baseName}"

    script:
    """
    el_gato.py --assembly ${fasta} --out ${fasta.baseName}
    """
}

// ======================================================
// ONT PIPELINE — LINEAR ASSEMBLY DAG
// ======================================================

process SUBSAMPLE_ONT {

    publishDir "${params.outdir}/ont/${params.run_id}/${sample_id}/subsample", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("subsampled.fastq.gz")

    script:
    """
    seqtk sample -s100 ${reads} 5000 | gzip > subsampled.fastq.gz
    """
}

process RAVEN {

    publishDir "${params.outdir}/ont/${params.run_id}/${sample_id}/assembly", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("assembly.fasta")

    script:
    """
    raven -t ${task.cpus} ${reads} > assembly.fasta
    """
}

process MEDAKA {

    publishDir "${params.outdir}/ont/${params.run_id}/${sample_id}/medaka", mode: 'copy'

    input:
    tuple val(sample_id), path(assembly), path(reads)

    output:
    tuple val(sample_id), path("consensus.fasta")

    script:
    """
    medaka_consensus \
        -i ${reads} \
        -d ${assembly} \
        -o medaka_out \
        -t ${task.cpus} \
        -b 40

    cp medaka_out/consensus.fasta .
    """
}

// ======================================================
// ONT FAN-OUT TOOLS
// ======================================================

process NANOPLOT {

    publishDir "${params.outdir}/ont/${params.run_id}/${sample_id}/nanoplot", mode: 'copy'
    time '30m'  // 30 minute timeout
    errorStrategy 'ignore'  // Don't fail pipeline if NanoPlot hangs

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("nanoplot/nanoplot/NanoPlot-report.html"), optional: true

    script:
    """
    echo "Starting NanoPlot for ${sample_id}"
    echo "Input file: ${reads}"
    echo "File size: \$(ls -lh ${reads})"
    
    # Check if file exists and is readable
    if [ ! -f "${reads}" ]; then
        echo "ERROR: Input file ${reads} does not exist"
        exit 1
    fi

    NanoPlot --fastq ${reads} -o nanoplot --threads ${task.cpus}
    
    echo "NanoPlot completed for ${sample_id}"
    """
}

process CHECKM {

    publishDir "${params.outdir}/ont/${params.run_id}/${sample_id}/checkm", mode: 'copy'

    input:
    tuple val(sample_id), path(fasta)

    output:
    tuple val(sample_id), path("checkm.tsv")

    script:
    """
    mkdir -p in
    cp ${fasta} in/

    checkm taxon_set genus Legionella legionella
    checkm analyze legionella in out -x fasta
    checkm qa legionella out -f checkm.tsv -o 2
    """
}

process PROKKA {

    publishDir "${params.outdir}/ont/${params.run_id}/${sample_id}/prokka", mode: 'copy'

    input:
    tuple val(sample_id), path(fasta)

    output:
    tuple val(sample_id), path("prokka/${sample_id}.gff"), emit: gff
    tuple val(sample_id), path("prokka/${sample_id}.faa"), emit: faa

    script:
    """
    prokka \
        --outdir prokka \
        --prefix ${sample_id} \
        --cpus ${task.cpus} \
        ${fasta}
    """
}

process AMRFINDER {

    publishDir "${params.outdir}/ont/${params.run_id}/${sample_id}/amrfinder", mode: 'copy'

    input:
    tuple val(sample_id), path(gff), path(faa)

    output:
    tuple val(sample_id), path("${sample_id}_amr.tsv")

    script:
    """
    amrfinder \
        -g ${gff} \
        -p ${faa} \
        -a prokka \
        --plus \
        -o ${sample_id}_amr.tsv
    """
}

process EL_GATO_ONT {

    publishDir "${params.outdir}/ont/${params.run_id}/${sample_id}/el_gato", mode: 'copy'

    input:
    tuple val(sample_id), path(fasta)

    output:
    tuple val(sample_id), path("report.json")

    script:
    """
    el_gato.py --assembly ${fasta} --out report.json
    """
}

// ======================================================
// DASHBOARD LAUNCH
// Passes the per-run results directory so the dashboard
// can auto-discover all samples within it.
// ======================================================

process LAUNCH_DASHBOARD {

    cache false

    input:
    val ready        // just a trigger — actual paths are resolved from outdir

    output:
    stdout

    script:
    def results_dir = "${launchDir}/${params.outdir}/ont/${params.run_id}"
    """
    streamlit run ${params.dashboard} \
        -- \
        --results_dir ${results_dir} \
        > ${launchDir}/streamlit.log 2>&1 &

    echo "Streamlit PID: \$!"
    echo "Dashboard : http://localhost:8501"
    echo "Results   : ${results_dir}"
    echo "Log       : ${launchDir}/streamlit.log"
    sleep 5
    """
}

// ======================================================
// WORKFLOW ROUTER
// ======================================================

workflow {

    // --------------------------------------------------
    // ILLUMINA PE FASTQ
    // --------------------------------------------------
    if ( params.read1 && params.read2 ) {

        Channel
            .of( tuple(params.sample_id, file(params.read1), file(params.read2)) )
            .set { illumina_ch }

        EL_GATO_FASTQ(illumina_ch)
    }

    // --------------------------------------------------
    // ILLUMINA / ASSEMBLY FASTA
    // --------------------------------------------------
    else if ( params.fasta ) {

        Channel
            .of( file(params.fasta) )
            .set { fasta_ch }

        EL_GATO_FASTA(fasta_ch)
    }

    // --------------------------------------------------
    // ONT — MULTI-SAMPLE VIA SAMPLESHEET
    // --------------------------------------------------
    else if ( params.samplesheet ) {

        log.info "ONT multi-sample pipeline | run_id: ${params.run_id}"

        sample_ch = Channel
            .fromPath(params.samplesheet, checkIfExists: true)
            .splitCsv(header: true)
            .map { row ->
                assert row.sample_id : "Missing sample_id in samplesheet row: ${row}"
                assert row.read1     : "Missing read1 in samplesheet row: ${row}"
                tuple(row.sample_id, file(row.read1, checkIfExists: true))
            }

        // ---- Linear assembly ----
        subsampled = SUBSAMPLE_ONT(sample_ch)
        assembled  = RAVEN(subsampled)
        polished   = MEDAKA( assembled.join(subsampled) )

        // ---- Fan-out from subsampled reads ----
        nanoplot_out  = NANOPLOT(subsampled)

        // ---- Fan-out from polished assembly ----
        checkm_out    = CHECKM(polished)
        prokka_out    = PROKKA(polished)
        elgato_out    = EL_GATO_ONT(polished)

        // ---- AMRFinder needs GFF + FAA from Prokka ----
        amr_out = AMRFINDER( prokka_out.gff.join(prokka_out.faa) )

        // ---- Wait for all samples to finish, then launch dashboard ----
        // collect() ensures we wait for every sample before triggering
        trigger_ch = checkm_out
            .mix( elgato_out, amr_out )  // Temporarily removed nanoplot_out
            .map { sid, f -> f.toString() }
            .collect()
            .map { "done" }

        LAUNCH_DASHBOARD(trigger_ch)
    }

    else {
        error """
        No valid input provided. Choose one of:
          --samplesheet  samplesheet.csv   (ONT, multi-sample)
          --read1 R1.fq --read2 R2.fq      (Illumina PE)
          --fasta assembly.fa              (assembly / Illumina FASTA)
        """
    }
}
