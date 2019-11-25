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
        string mapped_reads_ref;
    } CallVariantsParams;

    /**
      Map reads to a reference assembly
    */
    funcdef map_reads_to_reference(MapReadsParams params) returns (ReportResults output) authentication required;

    /**
      Call variants in a reference assembly, based on mapped reads
    */
    funcdef call_variants(CallVariantsParams params) returns (ReportResults output) authentication required;
};
