import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="mhcshrubs",
    version="0.2.0",
    author="Christiaan H. van Dorp",
    author_email="chvandorp@gmail.com",
    description="Find associations between MHC genotype and disease traits using MHC similarities",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/chvandorp/MHCshrubs",
    packages=setuptools.find_packages(),
    include_package_data=True,
    package_data={
        '' : ['stan/*.stan', 'jags/*.bug', 'resources/*']
    },
    entry_points={
        'console_scripts': [
          'mhcshrubs = mhcshrubs.main:main'
        ]
    },
    install_requires=[
        'numpy',
        'scipy',
        'PyQt5',
        'ete3',
        'networkx',
        'tqdm',
        'pystan',
        'matplotlib'
    ],
    classifiers=(
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ),
)
