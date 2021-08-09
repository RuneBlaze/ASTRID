import os
import re
import sys
import sysconfig
import platform
import subprocess
from pathlib import Path

from distutils.version import LooseVersion
from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext
from setuptools.command.test import test as TestCommand


class BazelExtension(Extension):
    def __init__(self, name):
        Extension.__init__(self, name, sources=[])


class BazelBuild(build_ext):
    def run(self):
        subprocess.run(['bazel', 'build', '//src:asterid.so', '-c', 'opt'])
        for ext in self.extensions:
            self.move_output(ext)

    def move_output(self, ext):
        # build_temp = Path(self.build_temp).resolve()
        dest_path = Path(self.get_ext_fullpath(ext.name)).resolve()
        dest_directory = dest_path.parents[0]
        dest_directory.mkdir(parents=True, exist_ok=True)
        self.copy_file("bazel-bin/src/asterid.so", dest_path)
        
        
ext_modules = [
  BazelExtension('asterid'),
]

setup(
    name="asterid",
    version="0.1.6.post1",
    author="runeblaze",
    author_email="runeblaze@excite.co.jp", 
  packages=find_packages(),
  ext_modules=ext_modules,
  cmdclass=dict(build_ext=BazelBuild),
  zip_safe=False,
  python_requires=">=3.6",
)