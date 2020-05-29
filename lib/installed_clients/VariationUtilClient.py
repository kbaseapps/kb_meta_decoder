# -*- coding: utf-8 -*-
############################################################
#
# Autogenerated by the KBase type compiler -
# any changes made here will be overwritten
#
############################################################

from __future__ import print_function
# the following is a hack to get the baseclient to import whether we're in a
# package or not. This makes pep8 unhappy hence the annotations.
try:
    # baseclient and this client are in a package
    from .baseclient import BaseClient as _BaseClient  # @UnusedImport
except ImportError:
    # no they aren't
    from baseclient import BaseClient as _BaseClient  # @Reimport


class VariationUtil(object):

    def __init__(
            self, url=None, timeout=30 * 60, user_id=None,
            password=None, token=None, ignore_authrc=False,
            trust_all_ssl_certificates=False,
            auth_svc='https://ci.kbase.us/services/auth/api/legacy/KBase/Sessions/Login',
            service_ver='dev',
            async_job_check_time_ms=100, async_job_check_time_scale_percent=150, 
            async_job_check_max_time_ms=300000):
        if url is None:
            raise ValueError('A url is required')
        self._service_ver = service_ver
        self._client = _BaseClient(
            url, timeout=timeout, user_id=user_id, password=password,
            token=token, ignore_authrc=ignore_authrc,
            trust_all_ssl_certificates=trust_all_ssl_certificates,
            auth_svc=auth_svc,
            async_job_check_time_ms=async_job_check_time_ms,
            async_job_check_time_scale_percent=async_job_check_time_scale_percent,
            async_job_check_max_time_ms=async_job_check_max_time_ms)

    def save_variation_from_vcf(self, params, context=None):
        """
        Save a variation (and trait?) object to Kbase given a reference genome, object output name,
        Variant Call Format (VCF) file, and sample attribute file.
        :param params: instance of type "save_variation_input" (## funcdef
           save_variation_from_vcf ## required input params:
           genome_or_assembly_ref: KBaseGenomes.Genome or
           KBaseGenomeAnnotations.Assembly object reference *** variation
           input data *** vcf_staging_file_path: path to location data
           associated with samples variation_object_name: output name for
           KBase variation object *** sample input data ***
           sample_attribute_ref: x/y/z reference to kbase sample attribute
           optional params: NA output report: report_name report_ref HTML
           visualization: Manhattan plot *** Visualization *** plot_maf:
           generate histogram of minor allele frequencies plot_hwe: generate
           histogram of Hardy-Weinberg Equilibrium p-values) -> structure:
           parameter "workspace_name" of String, parameter
           "genome_or_assembly_ref" of type "obj_ref" (An X/Y/Z style
           reference), parameter "vcf_staging_file_path" of type "filepath"
           (KBase file path to staging files), parameter
           "variation_object_name" of String, parameter
           "sample_attribute_ref" of type "obj_ref" (An X/Y/Z style reference)
        :returns: instance of type "save_variation_output" -> structure:
           parameter "variation_ref" of String, parameter "report_name" of
           String, parameter "report_ref" of String
        """
        return self._client.run_job('VariationUtil.save_variation_from_vcf',
                                    [params], self._service_ver, context)

    def export_variation_as_vcf(self, params, context=None):
        """
        Export KBase variation object as Variant Call Format (VCF) file
        :param params: instance of type "export_variation_input" (## funcdef
           export_variation_as_vcf ## required input params: Variation object
           reference optional params: NA output report: Shock id pointing to
           exported vcf file) -> structure: parameter "input_var_ref" of type
           "obj_ref" (An X/Y/Z style reference)
        :returns: instance of type "export_variation_output" -> structure:
           parameter "shock_id" of String
        """
        return self._client.run_job('VariationUtil.export_variation_as_vcf',
                                    [params], self._service_ver, context)

    def get_variation_as_vcf(self, params, context=None):
        """
        Given a reference to a variation object, and output name: return a Variant Call Format (VCF)
        file path and name.
        :param params: instance of type "get_variation_input" (## funcdef
           get_variation_as_vcf ## required input params: Variation object
           reference output file name optional params: NA output report: path
           to returned vcf name of variation object) -> structure: parameter
           "variation_ref" of type "obj_ref" (An X/Y/Z style reference),
           parameter "filename" of String
        :returns: instance of type "get_variation_output" -> structure:
           parameter "path" of type "filepath" (KBase file path to staging
           files), parameter "variation_name" of String
        """
        return self._client.run_job('VariationUtil.get_variation_as_vcf',
                                    [params], self._service_ver, context)

    def status(self, context=None):
        return self._client.run_job('VariationUtil.status',
                                    [], self._service_ver, context)
