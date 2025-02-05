from setuptools import setup, find_packages

setup(
    name="celegans_geno_converter",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[

        "pandas",

        ],
    entry_points={
        "console_scripts": [
            "celegans_geno_converter=celegans_geno_converter.main:main",
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3", 
        "License :: OSI Approved :: MIT License", 
        "Operating System :: OS Independent", 
    ], 
)

