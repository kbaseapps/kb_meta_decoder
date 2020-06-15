/*
A KBase module: kb_meta_decoder
*/

module kb_meta_decoder {
    typedef structure {
        string report_name;
        string report_ref;
    } ReportResults;

    typedef structure {
	string workspace_name;
	string workspace_id;
        string assembly_ref;
        string reads_ref;
    } MapReadsParams;

    typedef structure {
	string workspace_name;
	string workspace_id;
        string assembly_ref;
        string reads_ref;  /* should be: string mapped_reads_ref; */
	string output_vcf;
	int min_mapping_quality;
	int min_depth;
    } CallVariantsParams;

    typedef structure {
	string workspace_name;
	string workspace_id;
        string assembly_ref;
        string reads_ref;  /* should be: string mapped_reads_ref; */
    } CalcPopStatsParams;

    /**
      Map reads to a reference assembly.  Should save BAM-like object.
    */
    funcdef map_reads_to_reference(MapReadsParams params) returns (ReportResults output) authentication required;

    /**
      Call variants in a reference assembly.  Should be based on mapped reads (BAM file).
    */
    funcdef call_variants(CallVariantsParams params) returns (ReportResults output) authentication required;

    /**
      Calculates population statistics
    */
    funcdef calculate_population_statistics(CalcPopStatsParams params) returns (ReportResults output) authentication required;
};
