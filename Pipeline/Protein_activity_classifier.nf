#!/usr/bin/env nextflow

// Generate log filename based on current date
logFile = "${workflow.launchDir}/nextflow-${workflow.startTime.format('yyyyMMdd')}.log"

// Logging setup
log.info "Starting Workflow..."
log.info "Logging to: $logFile"

// Redirect Nextflow logs to the log file
workflow.log = file(logFile)

params.uniprot_code = ''
params.active_ligands = ''
params.inactive_ligands = ''
params.query_coverage_threshold = ''
params.identity_threshold = ''
params.gap_open_penalty = ''
params.gap_extend_penalty = ''
params.seed = ''

process runAligner {
    beforeScript 'chmod +x ./Aligner#.py'
    errorStrategy 'retry'
    maxRetries 3

    script:
    """
    ./Aligner#.py \\
        ${params.uniprot_code ? "--uniprot_code=${params.uniprot_code}" : ""} \\
        ${params.active_ligands ? "--active_ligands=${params.active_ligands}" : ""} \\
        ${params.inactive_ligands ? "--inactive_ligands=${params.inactive_ligands}" : ""} \\
        ${params.query_coverage_threshold ? "--query_coverage_threshold=${params.query_coverage_threshold}" : ""} \\
        ${params.identity_threshold ? "--identity_threshold=${params.identity_threshold}" : ""} \\
        ${params.gap_open_penalty ? "--gap_open_penalty=${params.gap_open_penalty}" : ""} \\
        ${params.gap_extend_penalty ? "--gap_extend_penalty=${params.gap_extend_penalty}" : ""}
    """
}

process runTensor {
    beforeScript 'chmod +x ./Tensor#.py'
    errorStrategy 'retry'
    maxRetries 3
    
    input:
    path 'temp_dir', from: 'temp_dir'
    
    script:
    """
    ./Tensor#.py
    """
}

process runRFModel {
    beforeScript 'chmod +x ./RF_model.py'
    errorStrategy 'retry'
    maxRetries 3
    
    input:
    path 'temp_dir', from: 'temp_dir'
    
    script:
    """
    ./RF_model.py ${params.seed ? "--seed=${params.seed}" : ""}
    """
}

workflow.onComplete {
    log.info "Workflow completed."
}
