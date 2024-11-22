from setuptools import setup, find_packages

setup(
    name="planter",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "duckdb",
        "pandas",
        "biopython",
        # add other dependencies here
    ],
    entry_points={
        "console_scripts": [
            # add CLI tools if needed
        ],
    },
)
