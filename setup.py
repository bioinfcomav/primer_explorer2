import setuptools
from pathlib import Path
from setuptools import find_packages
from primer_explorer._version import version

with open("README.md", "r") as fh:
    long_description = fh.read()

requirements = ['matplotlib', 'openpyxl', 'pyranges', 'pandas']
scripts = [str(f) for f in Path('./bin').glob('*.py')]

setuptools.setup(
    name="primer_explorer2",  # Replace with your own username
    version=version,
    author="P.Ziarsolo",
    author_email="pziarsolo@gmail.com",
    description="A small library to help K-seq users to select primers.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/pziarsolo/primer_explorer2",
    packages=find_packages(),
    package_dir={"primer_explorer.primer3": "primer_explorer/primer3"},
    install_requires=requirements,
    package_data={"primer_explorer.primer3": ["primer3_config/*.ds", "primer3_config/*.dh",
                                              "primer3_config/interpretations/*.ds",
                                              "primer3_config/interpretations/*.dh",
                                              "bin/*"]},
    scripts=scripts,
    license="GNU General Public License v3.0",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
