import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="mhcshrubs",
    version="1.0.0",
    author="Christiaan H. van Dorp",
    author_email="chvandorp@gmail.com",
    description="Find associations between MHC genotype and disease traits using MHC similarities",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/chvandorp/MHCshrubs",
    packages=setuptools.find_packages(),
    classifiers=(
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ),
)

