from distutils.core import setup
from setuptools.command.test import test as TestCommand
import sys
from geo2fastq.config import VERSION

DESCRIPTION = """
Download GEO sequencing experiments and process to map and create track hubs.
"""

class PyTest(TestCommand):
    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_args = []
        self.test_suite = True
    def run_tests(self):
        #import here, cause outside the eggs aren't loaded
        import pytest
        errno = pytest.main(self.test_args)
        sys.exit(errno)


setup (name = 'geo2fastq',
	version = VERSION,
	description = DESCRIPTION,
	author='Simon van Heeringen',
	author_email='s.vanheeringen@ncmls.ru.nl',
	license='MIT',
	packages=[
		'geo2fastq'
	],
	scripts=[
		"scripts/geo2fastq",
		"scripts/geo2trackhub",
	],
	
    data_files=[
        ('config', ["config/geo2fastq.yaml"]),
    ],
    #include_package_data=True,
    tests_require=['pytest'],
    cmdclass = {'test': PyTest},
)
