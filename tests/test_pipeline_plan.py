"""Command-plan tests; these do not execute Augur or infer phylogenies."""
from pathlib import Path
import copy
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from run_nextstrain import load_configuration, build_steps, check_settings


class PipelinePlanTests(unittest.TestCase):
    def setUp(self):
        self.config,self.settings=load_configuration(ROOT/'config/analysis.json','NA')
        self.steps=build_steps(self.config,self.settings,'NA',4)
        self.commands={step.name:step.argv for step in self.steps}

    def test_unconfirmed_parameters_block_execution(self):
        with self.assertRaisesRegex(ValueError,'parameters_confirmed'):
            check_settings(self.settings)

    def test_confirmed_na_configuration_checks(self):
        self.settings['parameters_confirmed']=True
        check_settings(self.settings)

    def test_ha_model_must_be_supplied(self):
        _,ha=load_configuration(ROOT/'config/analysis.json','HA')
        ha['parameters_confirmed']=True
        with self.assertRaisesRegex(ValueError,'substitution_model'):
            check_settings(ha)

    def test_reference_alignment_uses_existing_tip(self):
        self.assertIn('--reference-name',self.commands['align'])
        self.assertNotIn('--reference-sequence',self.commands['align'])

    def test_na_iqtree_model_and_seed(self):
        cmd=self.commands['tree']
        self.assertEqual(cmd[cmd.index('--substitution-model')+1],'GTR')
        self.assertIn('--tree-builder-args=-seed 666962',cmd)

    def test_iqtree_seed_not_assumed_for_treetime(self):
        self.assertNotIn('--seed',self.commands['refine'])
        self.assertNotIn('--seed',self.commands['ancestral'])

    def test_no_unsupported_seed_flag_for_traits(self):
        self.assertNotIn('--seed',self.commands['traits'])

    def test_refine_requests_per_site_divergence(self):
        cmd=self.commands['refine']
        self.assertEqual(cmd[cmd.index('--divergence-units')+1],'mutations-per-site')

    def test_inclusion_files_form_one_argument_group(self):
        self.settings['include_file']='private/include.txt'
        steps=build_steps(self.config,self.settings,'NA',4)
        filt=next(s.argv for s in steps if s.name=='filter')
        self.assertEqual(filt.count('--include'),1)
        self.assertEqual(len(filt[filt.index('--include')+1:]),2)

    def test_clinical_browser_export_blocked(self):
        self.settings['parameters_confirmed']=True
        self.settings['itol_columns']=['country','age']
        with self.assertRaisesRegex(ValueError,'iTOL columns'):
            check_settings(self.settings)

    def test_divergence_export_does_not_require_valid_temporal_export(self):
        itol=next(s for s in self.steps if s.name=='itol')
        self.assertTrue(any(p.name=='NA_divergence.nwk' for p in itol.outputs))
        self.assertFalse(any(p.name=='NA_time_years.nwk' for p in itol.outputs))

    def test_ha_patristic_step_uses_original_script(self):
        config,ha=load_configuration(ROOT/'config/analysis.json','HA')
        steps=build_steps(config,ha,'HA',4)
        step=next(s for s in steps if s.name=='patristic')
        self.assertEqual(Path(step.argv[1]).name,'patristic_distance.py')
        self.assertEqual(step.argv[step.argv.index('--reference')+1],'EPI_ISL_983345')


if __name__=='__main__':
    unittest.main()
