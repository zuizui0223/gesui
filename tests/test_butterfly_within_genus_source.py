"""Checks that prevent confusing genera, abundance and discordant recordings."""
import importlib.util
import sys
import unittest
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'analysis/20261122_butterfly_within_genus'))
try:
    from analyse import build_pair_strata, eligible_scores
except ImportError:
    build_pair_strata = eligible_scores = None


class SourceAssemblyTests(unittest.TestCase):
    def require_builder(self):
        self.assertIsNotNone(build_pair_strata, 'new discordant-pair builder is not implemented')

    def fixtures(self):
        cells = pd.DataFrame([
            {'h': h, 'h6': region, 'period': 'p1', 'plant_n_incidence': n, 'foo': 0.0}
            for h, region, n in [('1','r1',1),('2','r1',2),('3','r1',3),('4','r1',4)]
        ])
        scores = pd.DataFrame([
            {'butterflyKey': 'a', 'accepted_name':'Alpha one', 'genus':'Alpha', 'host_nutrient_mean':2.0},
            {'butterflyKey': 'b', 'accepted_name':'Alpha two', 'genus':'Alpha', 'host_nutrient_mean':7.0},
            {'butterflyKey': 'c', 'accepted_name':'Beta one', 'genus':'Beta', 'host_nutrient_mean':9.0},
        ])
        counts = pd.DataFrame([
            {'h':h,'period':'p1','butterflyKey':s,'count':n}
            for h,s,n in [('1','a',1),('2','b',1),('3','a',1),('3','b',1),('4','c',1)]
        ])
        return cells, counts, scores

    def test_only_same_genus_and_discordant_cells_are_used(self):
        self.require_builder()
        strata = build_pair_strata(*self.fixtures(), covariates=['foo'])
        self.assertEqual(len(strata),1)
        self.assertEqual(strata[0]['species_a'],'a')
        self.assertEqual(strata[0]['species_b'],'b')
        self.assertEqual(strata[0]['h'],['1','2'])
        self.assertEqual(strata[0]['y'],[1,0])

    def test_record_multiplicity_does_not_change_primary(self):
        self.require_builder()
        cells, counts, scores = self.fixtures()
        expected = build_pair_strata(cells,counts,scores,covariates=['foo'])
        counts.loc[counts.butterflyKey=='a','count']=100000
        self.assertEqual(expected,build_pair_strata(cells,counts,scores,covariates=['foo']))

    def test_constant_pair_outcome_is_not_informative(self):
        self.require_builder()
        cells, counts, scores = self.fixtures()
        counts = counts[~((counts.h=='2') & (counts.butterflyKey=='b'))]
        self.assertEqual(build_pair_strata(cells,counts,scores,covariates=['foo']),[])

    def test_tied_pair_retained_when_genus_has_other_trait_variation(self):
        self.require_builder()
        cells, counts, scores = self.fixtures()
        scores.loc[scores.butterflyKey=='b','host_nutrient_mean']=2.0
        scores.loc[scores.butterflyKey=='c',['accepted_name','genus']]=['Alpha three','Alpha']
        strata=build_pair_strata(cells,counts,scores,covariates=['foo'])
        ab=[s for s in strata if s['species_a']=='a' and s['species_b']=='b']
        self.assertEqual(len(ab),1)
        self.assertEqual(ab[0]['trait_difference'],0.0)

    def test_missing_nonaccepted_and_alias_conflict_are_not_scored(self):
        self.assertIsNotNone(eligible_scores, 'new score eligibility is not implemented')
        rows=[]
        for key,status,value,conflict in [('1','ACCEPTED',2.0,False),('2','SYNONYM',5.0,False),('3','ACCEPTED',None,False),('4','ACCEPTED',6.0,True)]:
            rows.append(dict(butterflyKey=key,accepted_name='Alpha '+key,taxon_status=status,host_nutrient_mean=value,actual_used_alias_conflict=conflict,trait_coverage=1.0))
        self.assertEqual(eligible_scores(pd.DataFrame(rows)).butterflyKey.tolist(),['1'])


if __name__=='__main__':
    unittest.main()
