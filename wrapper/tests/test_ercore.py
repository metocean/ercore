from unittest import TestCase
from nose.tools import assert_true
import os
import datetime as dt

import yaml

from msl_actions import stdout_logger

from msl_actions.modeling.ercore_wrap import ERCoreWrapper

class TestErcoreCase(TestCase):
    """Tests for ercore model wrapper"""

    def setUp(self):
        with open('model.ercore.nz.yaml') as yamlfile:
            self.config = yaml.load(yamlfile)
        cycle = dt.datetime.utcnow()-dt.timedelta(days=2)
        self.cycle_dt = dt.datetime(cycle.year,cycle.month,cycle.day,0)

    def test_instance(self):
        wrapper = ERCoreWrapper(**self.config)
        assert_true(wrapper)      

    def test_input_query(self):
        wrapper = ERCoreWrapper(**self.config)
        wrapper.set_cycle(dt.datetime(2000,1,1,0))
        query = wrapper._get_movers_query()
        assert_true(isinstance(query, dict))

    def test_uds_query(self):
        wrapper = ERCoreWrapper(**self.config)
        wrapper.set_cycle(self.cycle_dt)
        wrapper.get_input_data()
        assert_true(os.path.isfile(wrapper.inputfile))

    def test_ercore_config(self):
        wrapper = ERCoreWrapper(**self.config)
        wrapper.set_cycle(self.cycle_dt)
        wrapper.get_input_data()
        config = wrapper.get_ercore_config()
        with open('/tmp/ercore.yaml', 'w') as erconfig:
            yaml.dump(config, erconfig)

    def test_getsites(self):
        wrapper = ERCoreWrapper(**self.config)
        wrapper.set_cycle(self.cycle_dt)
        sites = wrapper.get_sites()
        assert_true(sites)

    def test_preprocess(self):
        wrapper = ERCoreWrapper(**self.config)
        wrapper.set_cycle(self.cycle_dt)
        wrapper.preprocess()
        assert_true(wrapper.config_files)

    def test_run_model(self):
        wrapper = ERCoreWrapper(**self.config)
        wrapper.set_cycle(self.cycle_dt)
        wrapper.run()
        assert_true(os.listdir(wrapper.outdir))