#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : ztf_forced_phot/setup.py
# License           : MIT
# Author            : Adam Miller
# Date              : Unspecified
# Last Modified Date: 20.11.2021
# Last Modified By  : v. Brinnel

from setuptools import setup # type: ignore[import]

setup(
	name = 'bts_phot',
	version = '0.1',
	description = 'Generates flux calibrated photometry from the ZTF forced photometry service provided by IPAC',
	author = 'Adam Miller',
	maintainer = 'Adam Miller',
	packages = ['bts_phot'],
	include_package_data = True,
	package_data = {'': ['cal_data/zp_thresholds_quadID.txt']},
	install_requires = ['pandas>=1.3.0', 'supersmoother'],
	python_requires = '>=3.6'
)
