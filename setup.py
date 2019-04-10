#!/usr/bin/env python
# -*- coding: utf-8 -*-


try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


with open('README.md') as readme_file:
    readme = readme_file.read()

with open('HISTORY.md') as history_file:
    history = history_file.read().replace('.. :changelog:', '')

with open('requirements.txt') as requirements_file:
    requirements = requirements_file.read()

test_requirements = [
    'pytest', 'coverage', "flake8"
]

setup(
    name='kh-tools',
    version='0.1.0',
    description="Kmer hashing tools contains data cleaning and visualization code for analyzing kmer-hashing similarity matrices",
    long_description=readme + '\n\n' + history,
    author="Olga Botvinnik",
    author_email='olga.botvinnik@czbiohub.org',
    url='https://github.com/czbiohub/kh-tools',
    packages=[
        'kh-tools',
    ],
    package_dir={'kh-tools':
                 'kh-tools'},
    include_package_data=True,
    install_requires=requirements,
    license="MIT",
    zip_safe=False,
    keywords='kh-tools',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ],
    entry_points={
        'console_scripts': [
            'kh-tools = kh-tools.commandline:cli'
        ]
    },
    test_suite='tests',
    tests_require=test_requirements
)
