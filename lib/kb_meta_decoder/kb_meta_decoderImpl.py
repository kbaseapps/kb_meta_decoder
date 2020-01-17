# -*- coding: utf-8 -*-
#BEGIN_HEADER
import sys
import logging
import os
import uuid
import subprocess
from pprint import pprint, pformat

from installed_clients.AssemblyUtilClient import AssemblyUtil as AUClient
from installed_clients.DataFileUtilClient import DataFileUtil as DFUClient
from installed_clients.ReadsUtilsClient import ReadsUtils
from installed_clients.KBaseReportClient import KBaseReport
#END_HEADER


class kb_meta_decoder:
    '''
    Module Name:
    kb_meta_decoder

    Module Description:
    A KBase module: kb_meta_decoder
    '''

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "0.0.1"
    GIT_URL = "git@github.com:kbaseapps/kb_meta_decoder.git"
    GIT_COMMIT_HASH = "b5ac2935324fdbde2415cc82099f915df158be24"

    #BEGIN_CLASS_HEADER

    def log(self, target, message):
        if target is not None:
            target.append(message)
        print(message)
        sys.stdout.flush()

    # get the contigs from the genome as FASTA
    def download_assembly(self, token, assembly_ref):
        SERVICE_VER = 'release'

        try:
            auClient = AUClient(self.callback_url, token=token, service_ver=SERVICE_VER)
        except Exception as e:
            raise ValueError('Unable to instantiate auClient with callback_url: '+ self.callback_url +' ERROR: ' + str(e))
        try:
            dfuClient = DFUClient(self.callback_url, token=token, service_ver=SERVICE_VER)
        except Exception as e:
            raise ValueError('Unable to instantiate dfuClient with callback_url: '+ self.callback_url +' ERROR: ' + str(e))

        contig_file = auClient.get_assembly_as_fasta({'ref':assembly_ref}).get('path')
        sys.stdout.flush()   # don't remember why this matters
        contig_file_path = dfuClient.unpack_file({'file_path': contig_file})['file_path']
        return contig_file_path

    # get the reads as FASTQ
    def download_reads(self, token, reads_ref):
        try:
            readsUtils_Client = ReadsUtils (url=self.callback_url, token=token)  # SDK local                   

            readsLibrary = readsUtils_Client.download_reads ({'read_libraries': [reads_ref],
                                                                 'interleaved': 'true'                                                              
            })
            reads_file_path = readsLibrary['files'][reads_ref]['files']['fwd']
        except Exception as e:
            raise ValueError('Unable to get reads library object from workspace: (' + reads_ref +")\n" + str(e))

        return reads_file_path

    # run samtools to map reads
    def map_reads(self, console, input_contigs, input_reads):
        try:
            # first index the contigs
            self.log(console,"Indexing contigs.\n");
            cmdstring = "/bwa/bwa index "+input_contigs

            cmdProcess = subprocess.Popen(cmdstring, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
            for line in cmdProcess.stdout:
                print(line)
            cmdProcess.wait()
            print('return code: ' + str(cmdProcess.returncode))
            if cmdProcess.returncode != 0:
                raise ValueError('Error running bwa index, return code: ' +
                                 str(cmdProcess.returncode) + '\n')

            # store output
            bam_file_path = os.path.join(self.scratch,"bam_file_"+str(uuid.uuid4())+".bam")

            # then map the reads
            self.log(console,"Mapping reads to contigs.\n");
            cmdstring = "/bwa/bwa mem "+input_contigs+" "+input_reads+"|samtools view -S -b >"+bam_file_path
            cmdProcess = subprocess.Popen(cmdstring, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
            for line in cmdProcess.stdout:
                print(line)
            cmdProcess.wait()
            print('return code: ' + str(cmdProcess.returncode)+ '\n')
            if cmdProcess.returncode != 0:
                raise ValueError('Error running samtools, return code: ' +
                                 str(cmdProcess.returncode) + '\n')

            if not os.path.isfile(bam_file_path) \
                or os.path.getsize (bam_file_path) == 0:
                raise ValueError('Error generating samtools output\n')

        except Exception as e:
            raise ValueError('Unable to map reads\n' + str(e))

        return bam_file_path

    # run samtools to get summary statistics
    def get_bam_stats(self, console, bam_file_path):
        try:
            # then map the reads
            self.log(console,"Getting bam stats.\n");
            cmdstring = "samtools flagstat "+bam_file_path
            cmdProcess = subprocess.Popen(cmdstring, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
            for line in cmdProcess.stdout:
                self.log(console,str(line))
            cmdProcess.wait()
            print('return code: ' + str(cmdProcess.returncode)+ '\n')
            if cmdProcess.returncode != 0:
                raise ValueError('Error running samtools, return code: ' +
                                 str(cmdProcess.returncode) + '\n')

        except Exception as e:
            raise ValueError('Unable to calculate bam stats\n' + str(e))

        return

    # sort bam
    def sort_bam(self, console, bam_file_path):
        try:
            # store output
            sorted_bam_file_path = os.path.join(self.scratch,"sorted_bam_file_"+str(uuid.uuid4())+".bam")

            # run sort
            self.log(console,"Sorting mapped reads.\n");
            cmdstring = "samtools sort "+bam_file_path+" -o "+sorted_bam_file_path

            cmdProcess = subprocess.Popen(cmdstring, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
            for line in cmdProcess.stdout:
                print(line)
            cmdProcess.wait()
            print('return code: ' + str(cmdProcess.returncode)+ '\n')
            if cmdProcess.returncode != 0:
                raise ValueError('Error running samtools sort, return code: ' +
                                 str(cmdProcess.returncode) + '\n')

        except Exception as e:
            raise ValueError('Unable to sort bam\n' + str(e))

        return sorted_bam_file_path

    # index bam
    def index_bam(self, console, bam_file_path):
        try:
            # run index
            self.log(console,"Indexing mapped reads.\n");
            cmdstring = "samtools index "+bam_file_path

            cmdProcess = subprocess.Popen(cmdstring, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
            for line in cmdProcess.stdout:
                print(line)
            cmdProcess.wait()
            print('return code: ' + str(cmdProcess.returncode)+ '\n')
            if cmdProcess.returncode != 0:
                raise ValueError('Error running samtools index, return code: ' +
                                 str(cmdProcess.returncode) + '\n')

        except Exception as e:
            raise ValueError('Unable to index bam\n' + str(e))

        return

    # call variants
    def call_variants_bcftools(self, console, contigs_file_path, sorted_bam_file_path):
        try:
            # store output
            vcf_file_path = os.path.join(self.scratch,"vcf_"+str(uuid.uuid4())+".vcf")

            # run mpileup
            self.log(console,"Calling variants.\n");
            cmdstring = "bcftools mpileup -Ou -f "+contigs_file_path+" "+sorted_bam_file_path+" | bcftools call -mv > "+vcf_file_path

            cmdProcess = subprocess.Popen(cmdstring, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
            for line in cmdProcess.stdout:
                print(line)
            cmdProcess.wait()
            print('return code: ' + str(cmdProcess.returncode)+ '\n')
            if cmdProcess.returncode != 0:
                raise ValueError('Error running bcftools, return code: ' +
                                 str(cmdProcess.returncode) + '\n')

        except Exception as e:
            raise ValueError('Unable to call variants\n' + str(e))

        return vcf_file_path

    # run bcftools to get summary statistics
    def get_vcf_stats(self, console, vcf_file_path):
        try:
            self.log(console,"Getting vcf stats.\n");
            cmdstring = "bcftools stats "+vcf_file_path
            cmdProcess = subprocess.Popen(cmdstring, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
            for line in cmdProcess.stdout:
                self.log(console,str(line))
            cmdProcess.wait()
            print('return code: ' + str(cmdProcess.returncode)+ '\n')
            if cmdProcess.returncode != 0:
                raise ValueError('Error running bcftools, return code: ' +
                                 str(cmdProcess.returncode) + '\n')

        except Exception as e:
            raise ValueError('Unable to calculate vcf stats\n' + str(e))

        return

    # set up paths for meta_decoder
    def setup_paths(self, console):
        try:
            self.log(console,"Setting up meta_decoder paths.\n");
            cmdstring = "cd /meta_decoder && git pull && ln -s /usr/bin/python3 bin/python && mkdir input_dir && cd input_dir && ln -s /kb/module/work/tmp/*.fastq . && ln -s /kb/module/work/*.fa . && cd .. && mkdir output_dir"
            cmdProcess = subprocess.Popen(cmdstring, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
            for line in cmdProcess.stdout:
                print(line)
            cmdProcess.wait()
            print('return code: ' + str(cmdProcess.returncode)+ '\n')
            if cmdProcess.returncode != 0:
                raise ValueError('Error running setup, return code: ' +
                                 str(cmdProcess.returncode) + '\n')

        except Exception as e:
            raise ValueError('Unable to set up paths\n' + str(e))

        return

    # run meta_decoder
    def run_meta_decoder(self, console):
        try:
            self.log(console,"Running meta_decoder.\n");
            cmdstring = "cd /meta_decoder && ./bin/python meta_decoder.py -i input_dir -inf .fastq --r input_dir --rf .fa --s 1 --o output_dir --t 1 --bwa /bwa/bwa && sh 0.sh"

            cmdProcess = subprocess.Popen(cmdstring, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
            for line in cmdProcess.stdout:
                print(line)
            cmdProcess.wait()
            print('return code: ' + str(cmdProcess.returncode)+ '\n')
            if cmdProcess.returncode != 0:
                raise ValueError('Error running meta_decoder, return code: ' +
                                 str(cmdProcess.returncode) + '\n')

        except Exception as e:
            raise ValueError('Unable to run meta_decoder\n' + str(e))

        return

    # make html
    def make_html(self, console):
        try:
            self.log(console,"Making HTML.\n");
            cmdstring = "cd /meta_decoder && ./bin/python meta_decoder.py --o output_dir --html T"

            cmdProcess = subprocess.Popen(cmdstring, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
            for line in cmdProcess.stdout:
                print(line)
            cmdProcess.wait()
            print('return code: ' + str(cmdProcess.returncode)+ '\n')
            if cmdProcess.returncode != 0:
                raise ValueError('Error running meta_decoder for html, return code: ' +
                                 str(cmdProcess.returncode) + '\n')

        except Exception as e:
            raise ValueError('Unable to run meta_decoder for html\n' + str(e))

        return


    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.scratch = os.path.abspath(config['scratch'])
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.shared_folder = config['scratch']
        logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
                            level=logging.INFO)
        #END_CONSTRUCTOR
        pass


    def map_reads_to_reference(self, ctx, params):
        """
        Map reads to a reference assembly.  Should save BAM-like object.
        :param params: instance of type "MapReadsParams" -> structure:
           parameter "workspace_name" of String, parameter "workspace_id" of
           String, parameter "assembly_ref" of String, parameter "reads_ref"
           of String
        :returns: instance of type "ReportResults" -> structure: parameter
           "report_name" of String, parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN map_reads_to_reference
        console = []
        self.log(console, 'Running map_reads_to_reference with parameters: ')
        self.log(console, "\n"+pformat(params))

        token = ctx['token']
        env = os.environ.copy()
        env['KB_AUTH_TOKEN'] = token

        # param checks
        required_params = ['workspace_name',
                           'workspace_id',
                           'assembly_ref',
                           'reads_ref'
                          ]
        for required_param in required_params:
            if required_param not in params or params[required_param] == None:
                raise ValueError ("Must define required param: '"+required_param+"'")
            
        # load provenance
        provenance = [{}]
        if 'provenance' in ctx:
            provenance = ctx['provenance']
        provenance[0]['input_ws_objects']=[str(params['assembly_ref']),str(params['reads_ref'])]

        # get the contigs from the genome as FASTA
        contigs_file_path = self.download_assembly(token, params['assembly_ref'])

        # get the reads as FASTQ
        reads_file_path = self.download_reads(token, params['reads_ref'])

        # run samtools
        print("got contigs as "+contigs_file_path)
        print("got reads as "+reads_file_path)
        bam_file_path = self.map_reads(console, contigs_file_path, reads_file_path)
        print("got bam output "+bam_file_path)

        # get bam stats
        self.get_bam_stats(console,bam_file_path)

        # build report
        #
        reportName = 'kb_map_reads_report_'+str(uuid.uuid4())

        reportObj = {'objects_created': [],
                     'message': "\n".join(console),
                     'direct_html': None,
                     'direct_html_link_index': None,
                     'file_links': [],
                     'html_links': [],
                     'html_window_height': 220,
                     'workspace_name': params['workspace_name'],
                     'report_object_name': reportName
                     }

        # text report
        try:
            reportObj['message'] = "\n".join(console)
            msg = "\n".join(console)
        except:
            raise ValueError ("no report generated")

        # save report object
        #
        SERVICE_VER = 'release'
        report = KBaseReport(self.callback_url, token=ctx['token'], service_ver=SERVICE_VER)
        report_info = report.create_extended_report(reportObj)

        output = { 'report_name': report_info['name'], 'report_ref': report_info['ref'] }

        #END map_reads_to_reference

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method map_reads_to_reference return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]

    def call_variants(self, ctx, params):
        """
        Call variants in a reference assembly.  Should be based on mapped reads (BAM file), and save VCF-like object.
        :param params: instance of type "CallVariantsParams" -> structure:
           parameter "workspace_name" of String, parameter "workspace_id" of
           String, parameter "assembly_ref" of String, parameter "reads_ref"
           of String
        :returns: instance of type "ReportResults" -> structure: parameter
           "report_name" of String, parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN call_variants
        console = []
        self.log(console, 'Running call_variants with parameters: ')
        self.log(console, "\n"+pformat(params))

        token = ctx['token']
        env = os.environ.copy()
        env['KB_AUTH_TOKEN'] = token

        # param checks
        required_params = ['workspace_name',
                           'workspace_id',
                           'assembly_ref',
                           'reads_ref'
                          ]
        for required_param in required_params:
            if required_param not in params or params[required_param] == None:
                raise ValueError ("Must define required param: '"+required_param+"'")
            
        # load provenance
        provenance = [{}]
        if 'provenance' in ctx:
            provenance = ctx['provenance']
        provenance[0]['input_ws_objects']=[str(params['assembly_ref']),str(params['reads_ref'])]

        # get the contigs from the genome as FASTA
        contigs_file_path = self.download_assembly(token, params['assembly_ref'])

        # get the reads as FASTQ
        reads_file_path = self.download_reads(token, params['reads_ref'])

        # run samtools
        print("got contigs as "+contigs_file_path)
        print("got reads as "+reads_file_path)
        bam_file_path = self.map_reads(console, contigs_file_path, reads_file_path)
        print("got bam output "+bam_file_path)

        # sort bam
        sorted_bam_file_path = self.sort_bam(console, bam_file_path)
        print("got sorted bam output "+sorted_bam_file_path)

        # index sorted bam
        self.index_bam(console,sorted_bam_file_path)
        print("indexed bam")

        # call variants using bcftools
        vcf_file_path = self.call_variants_bcftools(console, contigs_file_path, sorted_bam_file_path)
        print("got vcf output "+vcf_file_path)

        # get vcf stats
        self.get_vcf_stats(console,vcf_file_path)

        # build report
        #
        reportName = 'kb_call_variants_report_'+str(uuid.uuid4())

        reportObj = {'objects_created': [],
                     'message': "\n".join(console),
                     'direct_html': None,
                     'direct_html_link_index': None,
                     'file_links': [],
                     'html_links': [],
                     'html_window_height': 220,
                     'workspace_name': params['workspace_name'],
                     'report_object_name': reportName
                     }

        # text report
        try:
            reportObj['message'] = "\n".join(console)
            msg = "\n".join(console)
        except:
            raise ValueError ("no report generated")

        # save report object
        #
        SERVICE_VER = 'release'
        report = KBaseReport(self.callback_url, token=ctx['token'], service_ver=SERVICE_VER)
        report_info = report.create_extended_report(reportObj)

        output = { 'report_name': report_info['name'], 'report_ref': report_info['ref'] }

        #END call_variants

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method call_variants return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]

    def calculate_population_statistics(self, ctx, params):
        """
        Calculates population statistics
        :param params: instance of type "CalcPopStatsParams" -> structure:
           parameter "workspace_name" of String, parameter "workspace_id" of
           String, parameter "assembly_ref" of String, parameter "reads_ref"
           of String
        :returns: instance of type "ReportResults" -> structure: parameter
           "report_name" of String, parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN calculate_population_statistics
        console = []
        self.log(console, 'Running call_variants with parameters: ')
        self.log(console, "\n"+pformat(params))

        token = ctx['token']
        env = os.environ.copy()
        env['KB_AUTH_TOKEN'] = token

        # param checks
        required_params = ['workspace_name',
                           'workspace_id',
                           'assembly_ref',
                           'reads_ref'
                          ]
        for required_param in required_params:
            if required_param not in params or params[required_param] == None:
                raise ValueError ("Must define required param: '"+required_param+"'")
            
        # load provenance
        provenance = [{}]
        if 'provenance' in ctx:
            provenance = ctx['provenance']
        provenance[0]['input_ws_objects']=[str(params['assembly_ref']),str(params['reads_ref'])]

        # get the contigs from the genome as FASTA
        contigs_file_path = self.download_assembly(token, params['assembly_ref'])

        # get the reads as FASTQ
        reads_file_path = self.download_reads(token, params['reads_ref'])

        # set up meta_decoder paths
        self.setup_paths(console)

        # run meta_decoder
        self.run_meta_decoder(console)

        # make HTML
        self.make_html(console)

        # load in an example html file
        html_file_path = "/meta_decoder/output_dir/"+os.path.basename(reads_file_path)+"_"+os.path.basename(contigs_file_path).replace(".fa_assembly.fa",".fa")+".flt.vcf.Tajima.D.html"
        with open(html_file_path, "r") as html_file:
            html_output=html_file.readlines()

        # build report
        reportName = 'kb_calc_pop_stats_report_'+str(uuid.uuid4())

        reportObj = {'objects_created': [],
                     'message': "",
                     'direct_html': "\n".join(html_output),
                     'direct_html_link_index': None,
                     'file_links': [],
                     'html_links': [],
                     'html_window_height': 220,
                     'workspace_name': params['workspace_name'],
                     'report_object_name': reportName
                     }

        # text report
        try:
            reportObj['message'] = "\n".join(console)
            msg = "\n".join(console)
        except:
            raise ValueError ("no report generated")

        # save report object
        #
        SERVICE_VER = 'release'
        report = KBaseReport(self.callback_url, token=ctx['token'], service_ver=SERVICE_VER)
        report_info = report.create_extended_report(reportObj)

        output = { 'report_name': report_info['name'], 'report_ref': report_info['ref'] }

        #END calculate_population_statistics

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method calculate_population_statistics return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]
    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK",
                     'message': "",
                     'version': self.VERSION,
                     'git_url': self.GIT_URL,
                     'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
