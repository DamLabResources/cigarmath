#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

with open('README.md') as readme_file:
    readme = readme_file.read()

with open('HISTORY.md') as history_file:
    history = history_file.read()

requirements = [ ]

test_requirements = ['pytest>=3', ]

setup(
    author="Will Dampier",
    author_email='wnd22@drexel.edu',
    python_requires='>=3.6',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    description="Performs common operations on cigar strings",
    install_requires=requirements,
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='cigarmath',
    name='cigarmath',
    packages=find_packages(include=['cigarmath', 'cigarmath.*']),
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/judowill/cigarmath',
    version='0.1.0',
    zip_safe=False,
)
