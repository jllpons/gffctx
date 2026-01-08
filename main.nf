include { samplesheetToList; paramsSummaryLog } from 'plugin/nf-schema'

include { GFFCTX } from './workflows/gffctx.nf'

help_message = """
gffctx.nf  Annotate genomic intervals with genomic feature context
==================================================================

Usage:
    nextflow run main.nf --input samples.csv [options]

Required Arguments:
    --input <csv>               : CSV file with columns 'annotation' and 'intervals' specifying GFF files for multiple samples

Optional Arguments:

    --help                      : Print help message and exit
    --version                   : Print version and exit
"""

init_summary = """
G F F C T X   P I P E L I N E   v${params.manifest.version}
======================================
Run as                  : ${workflow.commandLine}
Started at              : ${workflow.start}
Config files            : ${workflow.configFiles}
Container images        : ${workflow.containerEngine}:${workflow.container}
"""

// DESC: Validate input arguments and initialize pipeline, printing a small summary
// ARGS: None, uses variables defined at the beginning of the script
// OUTS: `$init_summary` as log message at `INFO` level
//       `$help_message` as stdout if `--help` flag is set
//       `$version` as stdout if `--version` flag is set
//       Proper error message and exit code if required arguments are missing
// RETS: None
def validateParams() {
// `--help` and `--version` flags
    if (params.help) {
        println help_message
        System.exit(0)
    }
    if (params.version) {
        println "${params.manifest.name} v${params.manifest.version}"
        System.exit(0)
    }

// Check required arguments
    if (params.input == null) {
        println help_message
        log.error "Missing required argument: --input"
        System.exit(1)
    }
    if (!file(params.input).exists()) {
        log.error "File not found: ${params.input}"
        System.exit(1)
    }
    if (params.outdir == null) {
        println help_message
        log.error "Missing required argument: --outdir"
        System.exit(1)
    }
}



// DESC: Parse the samplesheet and convert it to a list
// ARGS: None, uses variables defined at the beginning of the script
// OUTS: Channel with the samplesheet converted to a list
// RETS: Channel with the samplesheet converted to a list
def parseSamplesheet() {
    def ch_samplesheet = Channel
        .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))

    return ch_samplesheet
}

// DESC: Display completion message based on workflow status
// ARGS: None, uses variables defined at the beginning of the script
// OUTS: Completion message at `INFO` or `ERROR` level
// RETS: None
def completionMsg() {

    if (workflow.success) {
        if (workflow.stats.ignoredCount == 0) {
            log.info "Pipeline completed successfully!"
        }
        else {
            log.info "Pipeline completed successully, but with errored processes"
        }
    }
    else {
        log.error "Pipeline completed with errors"
    }

}

// Main workflow
workflow {

    main:

    // Validate input parameters
    validateParams()

    // Initialization Summary - Everything looks good so far
    log.info init_summary
    // Print summary of supplied parameters
    log.info paramsSummaryLog(workflow)

    // Parse and validate samplesheet.csv from `--input`
    ch_samplesheet = parseSamplesheet()


    GFFCTX(
        ch_samplesheet
    )


    // Display any error encountered during the workflow
    workflow.onComplete {
        completionMsg()
    }
}
