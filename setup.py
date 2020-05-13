#!/usr/bin/env python
# -*- coding: utf-8 -*-


try:
    from setuptools import setup, find_packages
except ImportError:
    from distutils.core import setup, find_packages


with open("README.md") as readme_file:
    readme = readme_file.read()

with open("HISTORY.md") as history_file:
    history = history_file.read().replace(".. :changelog:", "")

with open("requirements.txt") as requirements_file:
    requirements = requirements_file.read()

test_requirements = ["pytest", "coverage", "flake8"]

setup(
    name="sencha",
    description="Sencha is a Python package for directly translating RNA-seq reads into coding protein sequence.",
    long_description=readme + "\n\n" + history,
    long_description_content_type="text/markdown",
    author="Olga Botvinnik",
    author_email="olga.botvinnik@czbiohub.org",
    url="https://github.com/czbiohub/sencha",
    packages=find_packages(
        exclude=["tests", "*.tests", "*.tests.*", "tests.*", "test_*"]
    ),
    include_package_data=True,
    install_requires=requirements,
    license="MIT",
    zip_safe=False,
    keywords="sencha",
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
    ],
    setup_requires=[
        "setuptools>=38.6.0",
        "setuptools_scm",
        "setuptools_scm_git_archive",
    ],
    entry_points={"console_scripts": ["sencha = sencha.commandline:cli"]},
    use_scm_version={"write_to": "sencha/version.py"},
    test_suite="tests",
    tests_require=test_requirements,
)
