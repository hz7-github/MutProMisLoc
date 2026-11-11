import json
import tempfile
import unittest
from pathlib import Path
from datetime import datetime, timezone
from unittest.mock import patch
from scripts._internal.identity_recovery import References

class RecoveryPolicyTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        h = self.root / 'hgnc.json'
        h.write_text(json.dumps({'response': {'docs': [{'hgnc_id':'HGNC:1','status':'Approved','symbol':'TEST','uniprot_ids':['P12345']}]}}))
        self.refs = References(h, self.root/'old', self.root/'new', False)
        self.entry = {'primaryAccession':'P12345','secondaryAccessions':['Q12345'],'entryType':'UniProtKB reviewed (Swiss-Prot)','organism':{'taxonId':9606},'genes':[{'geneName':{'value':'TEST'}}]}
    def tearDown(self):
        self.refs.close()
        self.tmp.cleanup()
    def work(self, *accessions):
        return {'row':{'Gene Name':'TEST','Accession_new':None},'originals':[{'record':{'uniprot_id':a}} for a in accessions]}
    def test_conflict(self):
        other = dict(self.entry, primaryAccession='Q99999')
        with patch.object(self.refs, 'get_entry', side_effect=lambda a: (self.entry if a=='P12345' else other,None)):
            self.assertEqual(self.refs.identify(self.work('P12345','Q99999'))[1], 'SOURCE_ACCESSIONS_CONFLICT')
    def test_secondary_mapping(self):
        with patch.object(self.refs,'get_entry',return_value=(self.entry,None)):
            self.assertEqual(self.refs.identify(self.work('P12345','Q12345'))[1],'IDENTITY_CONFIRMED')
    def test_unresolved_source(self):
        with patch.object(self.refs,'get_entry',return_value=(None,'OFFLINE_CACHE_MISS')):
            self.assertEqual(self.refs.identify(self.work('P12345'))[1],'SOURCE_ACCESSION_UNRESOLVED')
    def test_unreviewed_switch(self):
        self.entry['entryType']='UniProtKB unreviewed (TrEMBL)'
        with patch.object(self.refs,'get_entry',return_value=(self.entry,None)):
            self.assertEqual(self.refs.identify(self.work('P12345'))[1],'UNREVIEWED_NOT_ALLOWED')
            self.refs.include_unreviewed=True
            self.assertEqual(self.refs.identify(self.work('P12345'))[1],'IDENTITY_CONFIRMED')
    def test_cache_dates(self):
        self.assertFalse(self.refs.fresh({'fetched_at':'2000-01-01T00:00:00+00:00'}))
        self.assertFalse(self.refs.fresh({}))
        self.assertFalse(self.refs.fresh({'fetched_at':'bad'}))
        self.assertTrue(self.refs.fresh({'fetched_at':datetime.now(timezone.utc).isoformat()}))
    def test_expired_bulk_cache(self):
        import hashlib
        old=self.root/'old';old.mkdir()
        request={'url':'https://rest.uniprot.org/uniprotkb/P12345.json','params':{}}
        key=hashlib.sha256(json.dumps(request,sort_keys=True).encode()).hexdigest()
        (old/(key+'.json')).write_text(json.dumps({'schema':'uniprot-request-cache-v1','request':request,'status':'OK','payload':self.entry,'fetched_at':'2000-01-01T00:00:00+00:00'}))
        refs=References(self.root/'hgnc.json',old,self.root/'new',False)
        try:self.assertIsNone(refs.get_entry('P12345')[0])
        finally:refs.close()

if __name__ == '__main__':
    unittest.main()
