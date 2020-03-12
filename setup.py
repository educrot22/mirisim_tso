#!/usr/bin/env python

"""
Setup file for installing mirisim_tso
"""

import io
import os
import re

from setuptools import setup


def read(*names, **kwargs):
    with io.open(
        os.path.join(os.path.dirname(__file__), *names),
        encoding=kwargs.get("encoding", "utf8")
    ) as fp:
        return fp.read()


def find_version(*file_paths):
    version_file = read(*file_paths)
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]",
                              version_file, re.M)
    if version_match:
        return version_match.group(1)
    raise RuntimeError("Unable to find version string.")


setup(
    name="mirisim_tso",
    version=find_version("mirisim_tso", "version.py"),
    description="mirisim_tso",
    url="TODO",
    author="Marine Martin lagarde",
    author_email="miri@xxx.yyy",
    packages=['mirisim_tso'],
    package_dir={'mirisim_tso': 'mirisim_tso'},
    package_data={},
    scripts=[],
    data_files=[('', ['README.adoc'])]
)
