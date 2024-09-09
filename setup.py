from setuptools import setup

# TODO: the version from a file
__version__ = "0.0.0"

setup(
    name="planter",
    description="a panoply of assemblies",
    packages=["planter"],
    install_requires=["snakemake", "biopython", "build"],
    entry_points={"console_scripts": ["planter = planter.run:main"]},
    include_package_data=True,
    zip_safe=True,
    python_requires=">=3.6",
)
