# -*- coding: utf-8 -*-
#BEGIN_HEADER
import sys
import logging
import os
import uuid
import subprocess
import copy
import re
from pprint import pprint, pformat
from datetime import datetime

from installed_clients.WorkspaceClient import Workspace
from installed_clients.AssemblyUtilClient import AssemblyUtil as AUClient
from installed_clients.DataFileUtilClient import DataFileUtil as DFUClient
from installed_clients.KBaseReportClient import KBaseReport as ReportClient
from installed_clients.ReadsUtilsClient import ReadsUtils as RUClient
from installed_clients.VariationUtilClient import VariationUtil as VUClient
from installed_clients.KBParallelClient import KBParallel
from SetAPI.SetAPIServiceClient import SetAPI
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
    VERSION = "0.0.2"
    GIT_URL = "git@github.com:kbaseapps/kb_meta_decoder.git"
    GIT_COMMIT_HASH = "b4e83d791155b12efdf049d6021e1446b4bb0a71"

    #BEGIN_CLASS_HEADER
    SERVICE_VER = 'release'

    def log(self, target, message):
        if target is not None:
            target.append(message)
        print(message)
        sys.stdout.flush()

    # get the contigs from the genome as FASTA
    def download_assembly(self, token, assembly_ref):
        try:
            auClient = AUClient(self.callback_url, token=token, service_ver=self.SERVICE_VER)
        except Exception as e:
            raise ValueError('Unable to instantiate auClient with callback_url: '+ self.callback_url +' ERROR: ' + str(e))
        try:
            dfuClient = DFUClient(self.callback_url, token=token, service_ver=self.SERVICE_VER)
        except Exception as e:
            raise ValueError('Unable to instantiate dfuClient with callback_url: '+ self.callback_url +' ERROR: ' + str(e))

        contig_file = auClient.get_assembly_as_fasta({'ref':assembly_ref}).get('path')
        sys.stdout.flush()   # don't remember why this matters
        contig_file_path = dfuClient.unpack_file({'file_path': contig_file})['file_path']
        return contig_file_path

    # get the reads as FASTQ
    def download_reads(self, token, reads_ref):
        try:
            ruClient = RUClient(url=self.callback_url, token=token)
            readsLibrary = ruClient.download_reads ({'read_libraries': [reads_ref],
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
            cmdstring = "/usr/local/bin/bwa index "+input_contigs
            self.log(console,"command: "+cmdstring);
            cmdProcess = subprocess.Popen(cmdstring, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
            for line in cmdProcess.stdout:
                print(line.decode("utf-8").rstrip())
            cmdProcess.wait()
            print('return code: ' + str(cmdProcess.returncode))
            if cmdProcess.returncode != 0:
                raise ValueError('Error running bwa index, return code: ' +
                                 str(cmdProcess.returncode) + '\n')

            # store output
            bam_file_path = os.path.join(self.scratch,"bam_file_"+str(uuid.uuid4())+".bam")

            # then map the reads
            self.log(console,"Mapping reads to contigs.\n");
            cmdstring = "/usr/local/bin/bwa mem "+input_contigs+" "+input_reads+"|samtools view -S -b >"+bam_file_path
            self.log(console,"command: "+cmdstring);
            cmdProcess = subprocess.Popen(cmdstring, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
            for line in cmdProcess.stdout:
                print(line.decode("utf-8").rstrip())
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
            self.log(console,"command: "+cmdstring);
            cmdProcess = subprocess.Popen(cmdstring, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
            for line in cmdProcess.stdout:
                self.log(console,line.decode("utf-8").rstrip())
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
            self.log(console,"command: "+cmdstring);
            cmdProcess = subprocess.Popen(cmdstring, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
            for line in cmdProcess.stdout:
                print(line.decode("utf-8").rstrip())
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
            self.log(console,"command: "+cmdstring);
            cmdProcess = subprocess.Popen(cmdstring, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
            for line in cmdProcess.stdout:
                print(line.decode("utf-8").rstrip())
            cmdProcess.wait()
            print('return code: ' + str(cmdProcess.returncode)+ '\n')
            if cmdProcess.returncode != 0:
                raise ValueError('Error running samtools index, return code: ' +
                                 str(cmdProcess.returncode) + '\n')

        except Exception as e:
            raise ValueError('Unable to index bam\n' + str(e))

        return

    # call variants
    def call_variants_bcftools(self, console, contigs_file_path, sorted_bam_file_path, min_mapping_quality, min_depth):
        try:
            # store raw and filtered output
            raw_vcf_file_path = os.path.join(self.scratch,"vcf_"+str(uuid.uuid4())+".raw.vcf")
            vcf_file_path = os.path.join(self.scratch,"vcf_"+str(uuid.uuid4())+".flt.vcf")

            # run mpileup
            self.log(console,"Calling variants.\n");
            cmdstring = "bcftools mpileup -a FMT/AD -B -d3000 -q"+str(min_mapping_quality)+" -Ou -f "+contigs_file_path+" "+sorted_bam_file_path+" | bcftools call --ploidy 1 -m -A > "+raw_vcf_file_path
            self.log(console,"command: "+cmdstring);
            cmdProcess = subprocess.Popen(cmdstring, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
            for line in cmdProcess.stdout:
                print(line.decode("utf-8").rstrip())
            cmdProcess.wait()
            print('return code: ' + str(cmdProcess.returncode)+ '\n')
            if cmdProcess.returncode != 0:
                raise ValueError('Error running bcftools, return code: ' +
                                 str(cmdProcess.returncode) + '\n')

            # filter and call variants
            self.log(console,"Filtering variants.\n");
            cmdstring = "bcftools filter -s LowQual -e 'DP>"+str(min_depth)+"' "+raw_vcf_file_path+" | bcftools view -v snps > "+vcf_file_path
            self.log(console,"command: "+cmdstring);
            cmdProcess = subprocess.Popen(cmdstring, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
            for line in cmdProcess.stdout:
                print(line.decode("utf-8").rstrip())
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
            self.log(console,"command: "+cmdstring);
            cmdProcess = subprocess.Popen(cmdstring, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
            for line in cmdProcess.stdout:
                self.log(console,line.decode("utf-8").rstrip())
            cmdProcess.wait()
            print('return code: ' + str(cmdProcess.returncode)+ '\n')
            if cmdProcess.returncode != 0:
                raise ValueError('Error running bcftools, return code: ' +
                                 str(cmdProcess.returncode) + '\n')

        except Exception as e:
            raise ValueError('Unable to calculate vcf stats\n' + str(e))

        return

    # set up paths for meta_decoder
    def setup_paths(self, console, reads_file_path, contigs_file_path):
        self.log(console,"Setting up meta_decoder paths.\n");

        timestamp = int((datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds() * 1000)
        output_dir = os.path.join(self.scratch, 'output_' + str(timestamp))
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        try:
            cmdstring = "cd /meta_decoder && ln -s /usr/bin/python3 bin/python && mkdir input_dir && cd input_dir && ln -s "+reads_file_path+" . && ln -s "+contigs_file_path+" . && cd .. && ln -s "+output_dir+" output_dir"
            self.log(console,"command: "+cmdstring);
            cmdProcess = subprocess.Popen(cmdstring, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
            for line in cmdProcess.stdout:
                print(line.decode("utf-8").rstrip())
            cmdProcess.wait()
            print('return code: ' + str(cmdProcess.returncode)+ '\n')
            if cmdProcess.returncode != 0:
                raise ValueError('Error running setup, return code: ' +
                                 str(cmdProcess.returncode) + '\n')

        except Exception as e:
            raise ValueError('Unable to set up paths\n' + str(e))

        return output_dir

    # run meta_decoder
    def run_meta_decoder(self, console):
        try:
            self.log(console,"Running meta_decoder.\n");
            cmdstring = "cd /meta_decoder && ./bin/python meta_decoder.py -i input_dir -inf .fastq --r input_dir --rf .fa --s 1 --o output_dir --t 1 --bwa /usr/local/bin/bwa && /bin/sh sub_metadecoder/0.sh"
            self.log(console,"command: "+cmdstring);
            cmdProcess = subprocess.Popen(cmdstring, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
            for line in cmdProcess.stdout:
                print(line.decode("utf-8").rstrip())
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
            cmdstring = "cd /meta_decoder && ./bin/python meta_decoder.py --o output_dir --html T && sh meta.decoder.visual.sh"
            self.log(console,"command: "+cmdstring);
            cmdProcess = subprocess.Popen(cmdstring, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
            for line in cmdProcess.stdout:
                print(line.decode("utf-8").rstrip())
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
        self.workspace_url = config['workspace-url']
        self.service_wizard_url = config['srv-wiz-url']
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.shared_folder = config['scratch']
        logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
                            level=logging.INFO)

        if not os.path.exists(self.scratch):
            os.makedirs(self.scratch)
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

        output_files = []
        try:
            dfuClient = DFUClient(self.callback_url, token=token, service_ver=self.SERVICE_VER)
        except Exception as e:
            raise ValueError('Unable to instantiate dfuClient with callback_url: '+ self.callback_url +' ERROR: ' + str(e))
        dfu_output = dfuClient.file_to_shock({'file_path': bam_file_path,
                                              'make_handle': 0})
        output_files.append({'shock_id': dfu_output['shock_id'],
                             'name': os.path.basename(bam_file_path),
                             'label': 'BAM file',
                             'description': 'BAM file'})

        # get bam stats
        self.get_bam_stats(console,bam_file_path)

        # build report
        reportName = 'kb_map_reads_report_'+str(uuid.uuid4())

        reportObj = {'objects_created': [],
                     'message': "\n".join(console),
                     'direct_html': None,
                     'direct_html_link_index': None,
                     'file_links': output_files,
                     'html_links': [],
                     'workspace_name': params['workspace_name'],
                     'report_object_name': reportName
                     }

        # save report object
        try:
            reportClient = ReportClient(self.callback_url, token=ctx['token'], service_ver=self.SERVICE_VER)
            report_info = reportClient.create_extended_report(reportObj)
        except:
            raise ValueError ("no report generated")

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
        Call variants in a reference assembly.  Should be based on mapped reads (BAM file).
        :param params: instance of type "CallVariantsParams" -> structure:
           parameter "workspace_name" of String, parameter "workspace_id" of
           String, parameter "assembly_ref" of String, parameter "reads_refs"
           of list of String, parameter "min_mapping_quality" of Long,
           parameter "min_depth" of Long
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
                           'reads_refs',
                           'min_mapping_quality',
                           'min_depth'
                          ]
        for required_param in required_params:
            if required_param not in params or params[required_param] == None:
                raise ValueError ("Must define required param: '"+required_param+"'")
            
        # load provenance
        provenance = [{}]
        if 'provenance' in ctx:
            provenance = ctx['provenance']
        provenance[0]['input_ws_objects']=[str(params['assembly_ref'])]
        for reads_ref in params['reads_refs']:
            provenance[0]['input_ws_objects'].append(str(reads_ref))

        try:
            self.log(console,'instantiating KBParallel')
            parallelClient = KBParallel(self.callback_url)
        except:
            raise ValueError("unable to instantiate KBParallelClient")
        try:
            wsClient = Workspace(self.workspace_url, token=token)
        except Exception as e:
            raise ValueError("unable to instantiate wsClient. "+str(e))

        # unpack sets, if any
        all_reads_refs = []
        # object info
        [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11)  # object_info tuple
        for reads_ref in params['reads_refs']:
            try:
                reads_obj_info = wsClient.get_object_info_new({'objects':[{'ref':reads_ref}]})[0]
                reads_obj_type = reads_obj_info[TYPE_I]
                reads_obj_type = re.sub ('-[0-9]+\.[0-9]+$', "", reads_obj_type)  # remove trailing version
                if reads_obj_type == 'KBaseSets.ReadsSet':
                    # unpack it
                    try:
                        setAPI_Client = SetAPI(url=self.service_wizard_url, token=ctx['token'])
                        self.log(console,'getting reads set '+reads_ref)
                        readsSet_obj = setAPI_Client.get_reads_set_v1({'ref':reads_ref,'include_item_info':1})

                    except Exception as e:
                        raise ValueError('SetAPI FAILURE: Unable to get read library set object: (' + reads_ref+ ')\n' + str(e))
                    for readsLibrary_obj in readsSet_obj['data']['items']:
                        all_reads_refs.append(readsLibrary_obj['ref'])
                else:
                    # use other reads objects "as is"
                    all_reads_refs.append(reads_ref)
            except Exception as e:
                raise ValueError('Unable to get read library object: (' + str(reads_ref) +')' + str(e))

        # RUN call variants in parallel
        report_text = ''
        parallel_tasks = []
        for reads_ref in all_reads_refs:
            single_params = copy.deepcopy(params)
            del single_params['reads_refs']
            single_params['reads_ref'] = reads_ref

            parallel_tasks.append( {'module_name': 'kb_meta_decoder',
                                    'function_name': 'call_variants_single',
                                    'version': 'dev',
                                    'parameters': single_params
                                    } )
            self.log(console, "QUEUING calling variants for reads: " + reads_ref)
            self.log(console, "\n" + pformat(single_params))
            self.log(console, str(datetime.now()))

        # run them in batch
        batch_run_params = {'tasks': parallel_tasks,
                            'runner': 'parallel',
                            'concurrent_local_tasks': 1,
                            'concurrent_njsw_tasks': 0,
                            'max_retries': 2}

        self.log(console, "RUNNING all variant calling in parallel.")
        self.log(console, "\n" + pformat(batch_run_params))
        self.log(console, str(datetime.now()))
        kbparallel_results = parallelClient.run_batch(batch_run_params)

        for fd in kbparallel_results['results']:
           if fd['result_package']['error'] is not None:
               raise ValueError('kb_parallel failed to complete without throwing an error on at least one of the nodes.')

        ### combine individual reports
        output_files = []
        for fd in kbparallel_results['results']:
            # self.log(console, "FD = " + pformat(fd))
            this_report_ref = fd['result_package']['result'][0]['report_ref']

            try:
                this_report_obj = wsClient.get_objects([{'ref': this_report_ref}])[0]['data']
            except:
                raise ValueError("unable to fetch report: " + this_report_ref)
            # self.log(console, "REPORT = "+pformat(this_report_obj))
            report_text += this_report_obj['text_message']
            report_text += "\n\n"
            reads_pattern = '\'reads_ref\': \'(\d+/\d+/\d+)\''
            match = re.search(reads_pattern, this_report_obj['text_message']) 
            if match:
                reads_ref = match.group(1)
            else:
                reads_ref = "unknown"
            for link in this_report_obj['file_links']:
                output_files.append({'shock_id': link['URL'].split('/')[-1],
                                     'name': link['name'],
                                     'label': link['label']+' for '+reads_ref,
                                     'description': link['description']+' for '+reads_ref})

        ### build and save the combined report
        reportName = 'kb_call_variants_report_'+str(uuid.uuid4())
        reportObj = {'objects_created': [],
                     'message': report_text,
                     'direct_html': None,
                     'direct_html_link_index': None,
                     'file_links': output_files,
                     'html_links': [],
                     'workspace_name': params['workspace_name'],
                     'report_object_name': reportName
                     }

        SERVICE_VER = 'release'
        # save report object
        try:
            reportClient = ReportClient(self.callback_url, token=ctx['token'], service_ver=self.SERVICE_VER)
            report_info = reportClient.create_extended_report(reportObj)
        except:
            raise ValueError ("no report generated")

        ### construct the output to send back
        output = {'report_name': report_info['name'], 'report_ref': report_info['ref']}

        #END call_variants

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method call_variants return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]

    def call_variants_single(self, ctx, params):
        """
        Call variants in a reference assembly.  Should be based on mapped reads (BAM file).
        :param params: instance of type "CallVariantsSingleParams" ->
           structure: parameter "workspace_name" of String, parameter
           "workspace_id" of String, parameter "assembly_ref" of String,
           parameter "reads_ref" of String, parameter "min_mapping_quality"
           of Long, parameter "min_depth" of Long
        :returns: instance of type "ReportResults" -> structure: parameter
           "report_name" of String, parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN call_variants_single
        console = []
        self.log(console, 'Running call_variants_single with parameters: ')
        self.log(console, "\n"+pformat(params))

        token = ctx['token']
        env = os.environ.copy()
        env['KB_AUTH_TOKEN'] = token

        # param checks
        required_params = ['workspace_name',
                           'workspace_id',
                           'assembly_ref',
                           'reads_ref',
                           'min_mapping_quality',
                           'min_depth'
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
        vcf_file_path = self.call_variants_bcftools(console, contigs_file_path, sorted_bam_file_path, params['min_mapping_quality'], params['min_depth'])
        print("got vcf output "+vcf_file_path)

        # save VCF file
        try:
            print("NOT using VUClient to save VCF - need object type release in KBase")
            # vuClient = VUClient(self.callback_url)
            # vcf_result = vuClient.save_variation_from_vcf({
            #    'workspace_name': params['workspace_name'],
            #    'genome_or_assembly_ref': params['assembly_ref'],
            #    'vcf_staging_file_path': vcf_file_path,
            #    'variation_object_name': params['output_vcf']
            #})
            #vcf_object = vcf_result['variation_ref']
            #print("SUCCESS using VUClient to save VCF")
        except Exception as e:
            print("vuclient didn't work.")
            cmdstring = "cat /kb/module/work/tmp/*vcf.errors_summary*"
            self.log(console,"command: "+cmdstring);
            cmdProcess = subprocess.Popen(cmdstring, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
            for line in cmdProcess.stdout:
                print(line.decode("utf-8").rstrip())
            cmdProcess.wait()
            raise

        # save bam and vcf files for report
        output_files = []
        try:
            dfuClient = DFUClient(self.callback_url, token=token, service_ver=self.SERVICE_VER)
        except Exception as e:
            raise ValueError('Unable to instantiate dfuClient with callback_url: '+ self.callback_url +' ERROR: ' + str(e))
        dfu_output = dfuClient.file_to_shock({'file_path': sorted_bam_file_path,
                                              'make_handle': 0})
        output_files.append({'shock_id': dfu_output['shock_id'],
                             'name': os.path.basename(sorted_bam_file_path),
                             'label': 'BAM file (sorted)',
                             'description': 'BAM file (sorted)'})
        dfu_output = dfuClient.file_to_shock({'file_path': vcf_file_path,
                                              'make_handle': 0})
        output_files.append({'shock_id': dfu_output['shock_id'],
                             'name': os.path.basename(vcf_file_path),
                             'label': 'VCF file',
                             'description': 'VCF file'})

        # get vcf stats
        self.get_vcf_stats(console,vcf_file_path)

        # build report
        reportName = 'kb_call_variants_single_report_'+str(uuid.uuid4())

        # don't claim the VCF object, or it will overwrite
        # Ranjan's linked report.
        reportObj = {'objects_created': [], # {'ref':vcf_object}],
                     'message': "\n".join(console),
                     'direct_html': None,
                     'direct_html_link_index': None,
                     'file_links': output_files,
                     'html_links': [],
                     'workspace_name': params['workspace_name'],
                     'report_object_name': reportName
                     }

        # save report object
        try:
            reportClient = ReportClient(self.callback_url, token=ctx['token'], service_ver=self.SERVICE_VER)
            report_info = reportClient.create_extended_report(reportObj)
        except:
            raise ValueError ("no report generated")

        output = { 'report_name': report_info['name'], 'report_ref': report_info['ref'] }

        #END call_variants_single

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method call_variants_single return value ' +
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
        output_dir = self.setup_paths(console, reads_file_path, contigs_file_path)

        # run meta_decoder
        self.run_meta_decoder(console)

        # make HTML
        self.make_html(console)

        suffixes = {".sites.pi": "Nucleotide_diversity per site",
                    ".windowed.pi": "Nucleotide_diversity per 1000bp",
                    ".frq": "Allele frequency for each site",
                    ".frq.count": "Raw allele counts for each site",
                    ".TsTv": "Transition / Transversion ratio in bins of size 1000bp",
                    ".TsTv.summary": "Simple summary of all Transitions and Transversions",
                    ".TsTv.count": "Transition / Transversion ratio as a function of alternative allele count",
                    ".het": "Measure of heterozygosity on a per-individual basis",
                    ".hwe": "p-value for each site from a Hardy-Weinberg Equilibrium test",
                    ".Tajima.D": "Tajima D statistic in bins with size of 1000bp",
                    ".relatedness": "Relatedness statistic of unadjusted Ajk statistic<br>Expectation of Ajk is zero for individuals within a populations, and one for an individual with themselves",
                    ".snpden": "Number and density of SNPs in bins of size of 1000bp",
                    ".indel.hist": "Histogram file of the length of all indels (including SNPs)"}

        # need another dfu client
        try:
            dfuClient = DFUClient(self.callback_url, token=token, service_ver=self.SERVICE_VER)
        except Exception as e:
            raise ValueError('Unable to instantiate dfuClient with callback_url: '+ self.callback_url +' ERROR: ' + str(e))

        # remove BAM files, or output zip gets too big; save some.
        output_files = []

        data_file_path = os.path.join(output_dir,os.path.basename(reads_file_path)+"_"+os.path.basename(contigs_file_path)+".bam")
        if os.path.isfile(data_file_path):
            os.remove(data_file_path)
        else:
            print('missing expected file '+data_file_path+'\n')
        data_file_path = os.path.join(output_dir,os.path.basename(reads_file_path)+"_"+os.path.basename(contigs_file_path)+".sorted.bam")
        if os.path.isfile(data_file_path):
            dfu_output = dfuClient.file_to_shock({'file_path': data_file_path,
                                                  'make_handle': 0})
            os.remove(data_file_path)
            output_files.append({'shock_id': dfu_output['shock_id'],
                                 'name': os.path.basename(data_file_path),
                                 'label': 'BAM file (sorted)',
                                 'description': 'BAM file (sorted)'})
        else:
            print('missing expected file '+data_file_path+'\n')

        data_file_path = os.path.join(output_dir,os.path.basename(reads_file_path)+"_"+os.path.basename(contigs_file_path)+".sorted.bam.bai")
        if os.path.isfile(data_file_path):
            os.remove(data_file_path)
        else:
            print('missing expected file '+data_file_path+'\n')
        data_file_path = os.path.join(output_dir,os.path.basename(reads_file_path)+"_"+os.path.basename(contigs_file_path)+".raw.vcf")
        if os.path.isfile(data_file_path):
            dfu_output = dfuClient.file_to_shock({'file_path': data_file_path,
                                              'make_handle': 0})
            os.remove(data_file_path)
            output_files.append({'shock_id': dfu_output['shock_id'],
                                 'name': os.path.basename(data_file_path),
                                 'label': 'VCF file (raw)',
                                 'description': 'VCF file (raw)'})
        else:
            print('missing expected file '+data_file_path+'\n')

        data_file_path = os.path.join(output_dir,os.path.basename(reads_file_path)+"_"+os.path.basename(contigs_file_path)+".flt.vcf")
        if os.path.isfile(data_file_path):
            dfu_output = dfuClient.file_to_shock({'file_path': data_file_path,
                                                  'make_handle': 0})
            os.remove(data_file_path)
            output_files.append({'shock_id': dfu_output['shock_id'],
                                 'name': os.path.basename(data_file_path),
                                 'label': 'VCF file (filtered)',
                                 'description': 'VCF file (filtered)'})
        else:
            print('missing expected file '+data_file_path+'\n')

        # make index/explanation of HTML output files
        # and load output
        output_html = []
        html_message = "<p><b>MetaDecoder output</b>:<p>Meta_decoder maps a metagenome to reference genomes and infer the intraspecies strain-level diversity in a sample. It finds maximum likelihood estimates for the strain genotypes and the strain abundances. Based on the strain genotypes and abundances, it estimates can compares the strain-level diversity across metagenome samples by calculating the nucleotide diversity, allele frequency, Transition / Transversion ratio, SNP density, and other statistic measurements.<p><ul>"
        for suffix, description in suffixes.items():
            data_file_base = os.path.basename(reads_file_path)+"_"+os.path.basename(contigs_file_path)+".flt.vcf"+suffix
            data_file_path = os.path.join(output_dir,data_file_base)
            html_file_base = data_file_base+".html"
            html_file_path = os.path.join(output_dir,html_file_base)
            html_message+='<li><a href="'+html_file_base+'">'+suffix+' file</a>: '+description+"\n"
            if os.path.isfile(data_file_path):
                print("Loading "+data_file_path+"\n")
                try:
                    dfu_output = dfuClient.file_to_shock({'file_path': data_file_path,
                                                          'make_handle': 0})
                    os.remove(data_file_path)
                    output_files.append({'shock_id': dfu_output['shock_id'],
                                         'name': 'meta_decoder_output'+suffix,
                                         'label': suffix+' file'})
                except:
                    print("failed to load "+data_file_path+"\n")
            if os.path.isfile(html_file_path):
                print("Indexing "+html_file_path+"\n")
                try:
                    new_buf = ['<body><a href="kb_meta_decoder_report_index.html">Back to index</a><p>']
                    with open (html_file_path, 'r') as html_handle:
                        for line in html_handle.readlines():
                            new_buf.append(line.lstrip())
                    new_buf.append('</body>')
                    with open (html_file_path, 'w') as html_handle:
                        for line in new_buf:
                            html_handle.write(line)
                except:
                    print("failed to index "+html_file_path+"\n")

        html_message+="</ul>\n"
        html_file_path = os.path.join(output_dir, 'kb_meta_decoder_report_index.html')
        with open(html_file_path, 'w') as html_handle:
            html_handle.write(html_message)
        try:
            dfu_output = dfuClient.file_to_shock({'file_path': html_file_path,
                                                  'pack': 'zip',
                                                  'make_handle': 0})
            output_html.append({'shock_id': dfu_output['shock_id'],
                                'name': 'kb_meta_decoder_report_index.html',
                                'label': 'All meta_decoder reports'})
        except:
            raise ValueError('Logging exception loading html_report to shock')

        # build report
        reportName = 'kb_calc_pop_stats_report_'+str(uuid.uuid4())

        reportObj = {'objects_created': [],
                     'direct_html': None,
                     'direct_html_link_index': 0,
                     'file_links': output_files,
                     'html_links': output_html,
                     'workspace_name': params['workspace_name'],
                     'report_object_name': reportName
                     }

        # save report object
        try:
            reportClient = ReportClient(self.callback_url, token=ctx['token'], service_ver=self.SERVICE_VER)
            report_info = reportClient.create_extended_report(reportObj)
        except:
            raise ValueError ("no report generated")

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
