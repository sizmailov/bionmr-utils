import setuptools

with open("README.rst", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="bionmr-utils",
    version="0.0.2",
    author="Sergei Izmailov",
    author_email="sergei.a.izmailov@gmail.com",
    description="Aggregate package of internal bionmr utils",
    long_description=long_description,
    long_description_content_type="text/x-rst",
    url="https://github.com/sizmailov/bionmr-utils",
    include_package_data = True,
    packages=setuptools.find_packages(),
    classifiers=(
        "Programming Language :: Python :: 3",
        "Operating System :: POSIX :: Linux",
    ),
    install_requires=[
        "pyxmolpp2>=0.7.0,<1.0",
        "pyxmolpp2-stubs",
        "tqdm",
        "pandas",
    ],
    package_data = {
        "bionmr_utils": [ 'data/rename_tables/*.csv' ]
    }
)