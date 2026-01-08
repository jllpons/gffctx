import groovy.json.JsonBuilder
// DESC: Generate QC metadata in JSON format
// ARGS: None
// OUTS: Generates a qc_metadata.json file in the `${params.outdir}/qc/` directory
// RETS: None
def generateQcMetadata() {
    def metadata = [
        // Identify
        pipeline_name: workflow.manifest.name,
        session_id: workflow.sessionId.toString(),

        // Versions
        pipeline_version: workflow.manifest.version,
        nextflow_version: nextflow.version.toString(),

        // Runtime
        run_name: workflow.runName,
        started_at: workflow.start.format("yyyy-MM-dd'T'HH:mm:ss"),

        // User and execution
        launched_by: System.getProperty("user.name"),
        command_line: workflow.commandLine,

        // Directories
        launch_dir: workflow.launchDir.toString(),
        work_dir: workflow.workDir.toString(),
        project_dir: workflow.projectDir.toString(),

        // Configuration
        config_files: workflow.configFiles.collect { it.toString() },
        container_engine: workflow.containerEngine ?: 'none',
        profile: workflow.profile
    ]

    def json = new JsonBuilder(metadata)

    file("${params.outdir}").mkdirs()
    def metadata_file = file("${params.outdir}/qc_metadata.json")
    metadata_file.text = json.toPrettyString()
}
