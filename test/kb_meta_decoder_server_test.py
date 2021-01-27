# -*- coding: utf-8 -*-
import os
import time
import unittest
import sys
from configparser import ConfigParser
from pprint import pprint
import requests

from kb_meta_decoder.kb_meta_decoderImpl import kb_meta_decoder
from kb_meta_decoder.kb_meta_decoderServer import MethodContext
from kb_meta_decoder.authclient import KBaseAuth as _KBaseAuth

from installed_clients.WorkspaceClient import Workspace
from installed_clients.AssemblyUtilClient import AssemblyUtil
from installed_clients.AbstractHandleClient import AbstractHandle as HandleService

class kb_meta_decoderTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.token = os.environ.get('KB_AUTH_TOKEN', None)
        config_file = os.environ.get('KB_DEPLOYMENT_CONFIG', None)
        cls.cfg = {}
        config = ConfigParser()
        config.read(config_file)
        for nameval in config.items('kb_meta_decoder'):
            cls.cfg[nameval[0]] = nameval[1]
        # Getting username from Auth profile for token
        authServiceUrl = cls.cfg['auth-service-url']
        auth_client = _KBaseAuth(authServiceUrl)
        user_id = auth_client.get_user(cls.token)
        # WARNING: don't call any logging methods on the context object,
        # it'll result in a NoneType error
        cls.ctx = MethodContext(None)
        cls.ctx.update({'token': cls.token,
                        'user_id': user_id,
                        'provenance': [
                            {'service': 'kb_meta_decoder',
                             'method': 'please_never_use_it_in_production',
                             'method_params': []
                             }],
                        'authenticated': 1})
        cls.wsURL = cls.cfg['workspace-url']
        cls.shockURL = cls.cfg['shock-url']
        cls.hs = HandleService(url=cls.cfg['handle-service-url'],
                               token=cls.token)
        cls.wsClient = Workspace(cls.wsURL)
        cls.serviceImpl = kb_meta_decoder(cls.cfg)
        cls.scratch = cls.cfg['scratch']
        cls.callback_url = os.environ['SDK_CALLBACK_URL']
        suffix = int(time.time() * 1000)
        cls.wsName = "test_meta_decoder_" + str(suffix)
        cls.wsinfo = cls.wsClient.create_workspace({'workspace': cls.wsName})
        print('created workspace ' + cls.getWsName())
        cls.wsID = cls.getWsID()
        print('wsID is '+str(cls.wsID))
        cls.au = AssemblyUtil(cls.callback_url, token=cls.token, service_ver='release')
        cls.setupTestData()

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')
        if hasattr(cls, 'nodes_to_delete'):
            for node in cls.nodes_to_delete:
                cls.delete_shock_node(node)
        if hasattr(cls, 'handles_to_delete'):
            if cls.handles_to_delete:
                cls.hs.delete_handles(cls.hs.hids_to_handles(cls.handles_to_delete))
                print('Deleted handles ' + str(cls.handles_to_delete))

    @classmethod
    def getWsName(cls):
        return cls.wsinfo[1]

    @classmethod
    def getWsID(cls):
        return cls.wsinfo[0]

    def getImpl(self):
        return self.serviceImpl

    @classmethod
    def delete_shock_node(cls, node_id):
        header = {'Authorization': 'Oauth {0}'.format(cls.token)}
        requests.delete(cls.shockURL + '/node/' + node_id, headers=header,
                        allow_redirects=True)
        print('Deleted shock node ' + node_id)

    # Helper script borrowed from the transform service, logger removed              
    @classmethod
    def upload_file_to_shock(cls, file_path):
        """                                                                          
        Use HTTP multi-part POST to save a file to a SHOCK instance.                 
        """

        header = dict()
        header["Authorization"] = "Oauth {0}".format(cls.token)

        if file_path is None:
            raise Exception("No file given for upload to SHOCK!")

        with open(os.path.abspath(file_path), 'rb') as dataFile:
            files = {'upload': dataFile}
            print('POSTing data')
            response = requests.post(
                cls.shockURL + '/node', headers=header, files=files,
                stream=True, allow_redirects=True)
            print('got response')

        if not response.ok:
            response.raise_for_status()

        result = response.json()

        if result['error']:
            raise Exception(result['error'][0])
        else:
            return result["data"]

    @classmethod
    def upload_file_to_shock_and_get_handle(cls, test_file):
        '''                                                                          
        Uploads the file in test_file to shock and returns the node and a            
        handle to the node.                                                          
        '''
        print('loading file to shock: ' + test_file)
        node = cls.upload_file_to_shock(test_file)
        pprint(node)
        cls.nodes_to_delete.append(node['id'])

        print('creating handle for shock id ' + node['id'])
        handle_id = cls.hs.persist_handle({'id': node['id'],
                                           'type': 'shock',
                                           'url': cls.shockURL
                                         })
        cls.handles_to_delete.append(handle_id)

        md5 = node['file']['checksum']['md5']
        return node['id'], handle_id, md5, node['file']['size']

    @classmethod
    def upload_reads(cls, wsobjname, object_body, fwd_reads,
                     rev_reads=None, single_end=False, sequencing_tech='Illumina',
                     single_genome='1'):

        ob = dict(object_body)  # copy                                               
        ob['sequencing_tech'] = sequencing_tech
        ob['wsname'] = cls.getWsName()
        ob['name'] = wsobjname
        if single_end or rev_reads:
            ob['interleaved'] = 0
        else:
            ob['interleaved'] = 1
        print('\n===============staging data for object ' + wsobjname +
              '================')
        print('uploading forward reads file ' + fwd_reads['file'])
        fwd_id, fwd_handle_id, fwd_md5, fwd_size = \
            cls.upload_file_to_shock_and_get_handle(fwd_reads['file'])

        ob['fwd_id'] = fwd_id
        rev_id = None
        rev_handle_id = None
        if rev_reads:
            print('uploading reverse reads file ' + rev_reads['file'])
            rev_id, rev_handle_id, rev_md5, rev_size = \
                cls.upload_file_to_shock_and_get_handle(rev_reads['file'])
            ob['rev_id'] = rev_id
        obj_ref = cls.readUtilsImpl.upload_reads(ob)
        objdata = cls.wsClient.get_object_info_new({
            'objects': [{'ref': obj_ref['obj_ref']}]
            })[0]
        cls.staged[wsobjname] = {'info': objdata,
                                 'ref': cls.make_ref(objdata),
                                 'fwd_node_id': fwd_id,
                                 'rev_node_id': rev_id,
                                 'fwd_handle_id': fwd_handle_id,
                                 'rev_handle_id': rev_handle_id
                                }

    @classmethod
    def upload_assembly(cls, a_name, fa_file_path):
        wsname = cls.getWsName()
        print('\n===============staging data for object ' + a_name +
              '================')
        cls.au.save_assembly_from_fasta(
            {'file': {'path': fa_file_path},
             'workspace_name': wsname,
             'assembly_name': a_name})
        cls.staged[a_name] = {'info': objdata,
                              'ref': cls.make_ref(objdata)}

    @classmethod
    def setupTestData(cls):
        print('Shock url ' + cls.shockURL)
        # print('WS url ' + cls.wsClient.url)                                                       
        # print('Handle service url ' + cls.hs.url)                                                 
        print('staging data')

        mapped_reads = {'file': 'data/mapped_reads.fastq',
                        'name': 'mapped_reads.fastq',
                        'type': 'fastq'}
        cls.upload_reads('mapped_reads', {'single_genome': 0}, mapped_reads)
        # cls.upload_assembly('assembly_nohits', 'data/genome_fragment_nohits.fa');
        cls.upload_assembly('assembly_withhits', 'data/genome_fragment_withhits.fa');

        print('Data staged.')


    @classmethod
    def make_ref(self, object_info):
        return str(object_info[6]) + '/' + str(object_info[0]) + \
            '/' + str(object_info[4])


    @unittest.skip("call variant test takes 10 min")
    def test_call_variants(self):
        # Prepare test objects in workspace if needed using
        # self.getWsClient().save_objects({'workspace': self.getWsName(),
        #                                  'objects': []})
        #
        # Run your method by
        # ret = self.getImpl().your_method(self.getContext(), parameters...)
        #
        # Check returned data with
        # self.assertEqual(ret[...], ...) or other unittest methods

        # these test objects are in appdev:
        test_ws_name = "jmc:narrative_1576867921697"
        test_ws_id = "35222"
        test_assembly = "35222/2/1"
        test_reads = "35222/3/1"
        # ret = self.serviceImpl.calculate_population_statistics(self.ctx, {'workspace_name': test_ws_name,
        ret = self.serviceImpl.call_variants_single(self.ctx, {'workspace_name': test_ws_name,
                                                        'workspace_id': test_ws_id,
                                                        'assembly_ref' : test_assembly,
                                                        'reads_ref' : test_reads,
                                                        'min_mapping_quality' : '30',
                                                        'min_depth' : '50'})

    @unittest.skip("parallel version is working")
    def test_call_variants_readsset(self):
        # Prepare test objects in workspace if needed using
        # self.getWsClient().save_objects({'workspace': self.getWsName(),
        #                                  'objects': []})
        #
        # Run your method by
        # ret = self.getImpl().your_method(self.getContext(), parameters...)
        #
        # Check returned data with
        # self.assertEqual(ret[...], ...) or other unittest methods

        # these test objects are in appdev:
        test_ws_name = "jmc:narrative_1576867921697"
        test_ws_id = "35222"
        test_assembly = "35222/2/1"
        test_reads = "35222/93/1"
        # ret = self.serviceImpl.calculate_population_statistics(self.ctx, {'workspace_name': test_ws_name,
        ret = self.serviceImpl.call_variants(self.ctx, {'workspace_name': test_ws_name,
                                                        'workspace_id': test_ws_id,
                                                        'assembly_ref' : test_assembly,
                                                        'reads_refs' : [test_reads],
                                                        'min_mapping_quality' : '30',
                                                        'min_depth' : '50'})
    @unittest.skip("parallel version is working")
    def test_call_variants_readsset_multi(self):
        # Prepare test objects in workspace if needed using
        # self.getWsClient().save_objects({'workspace': self.getWsName(),
        #                                  'objects': []})
        #
        # Run your method by
        # ret = self.getImpl().your_method(self.getContext(), parameters...)
        #
        # Check returned data with
        # self.assertEqual(ret[...], ...) or other unittest methods

        # these test objects are in appdev:
        test_ws_name = "jmc:narrative_1576867921697"
        test_ws_id = "35222"
        test_assembly = "35222/2/1"
        test_reads = ["35222/93/1", "35222/103/1"]
        # ret = self.serviceImpl.calculate_population_statistics(self.ctx, {'workspace_name': test_ws_name,
        ret = self.serviceImpl.call_variants(self.ctx, {'workspace_name': test_ws_name,
                                                        'workspace_id': test_ws_id,
                                                        'assembly_ref' : test_assembly,
                                                        'reads_refs' : test_reads,
                                                        'min_mapping_quality' : '30',
                                                        'min_depth' : '50'})

    @unittest.skip("current focus on calling variants")
    def test_calc_pop_stat(self):
        # these test objects are in appdev:
        test_ws_name = "jmc:narrative_1576867921697"
        test_ws_id = "35222"
        test_assembly = "35222/2/1"
        test_reads = "35222/3/1"
        ret = self.serviceImpl.calculate_population_statistics(self.ctx, {'workspace_name': test_ws_name,
                                                        'workspace_id': test_ws_id,
                                                        'assembly_ref' : test_assembly,
                                                        'reads_ref' : test_reads})


    def test_call_variants_small(self):
        # Prepare test objects in workspace if needed using
        # self.getWsClient().save_objects({'workspace': self.getWsName(),
        #                                  'objects': []})
        #
        # Run your method by
        # ret = self.getImpl().your_method(self.getContext(), parameters...)
        #
        # Check returned data with
        # self.assertEqual(ret[...], ...) or other unittest methods

        test_ws_name = self.wsName
        test_ws_id = self.wsID
        test_assembly = self.staged['assembly_withhits']['ref']
        test_reads = [self.staged['mapped_reads']['ref']]
        # ret = self.serviceImpl.calculate_population_statistics(self.ctx, {'workspace_name': test_ws_name,
        ret = self.serviceImpl.call_variants(self.ctx, {'workspace_name': test_ws_name,
                                                        'workspace_id': test_ws_id,
                                                        'assembly_ref' : test_assembly,
                                                        'reads_refs' : test_reads,
                                                        'min_mapping_quality' : '30',
                                                        'min_depth' : '50'})
