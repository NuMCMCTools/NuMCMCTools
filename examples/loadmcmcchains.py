#!/usr/bin/env python3

import os
from numcmctools import MCMCSamples

script_dir = os.path.dirname(os.path.abspath(__file__))
root_file_path = os.path.join(script_dir, "testchaindata.root")

chain = MCMCSamples(root_file_path, "mcmc")
