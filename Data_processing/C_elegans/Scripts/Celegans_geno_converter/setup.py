from setuptools import setup, find_packages

setup(
    name="Celegans_geno_converter",
    version="0.1.0",
    packages=find_packages(),
    install_requires=["pandas"],
    entry_points={
        "console_scripts": [
            "Celegans_geno_converter=Celegans_geno_converter.main:main",
        ],
    },
)

