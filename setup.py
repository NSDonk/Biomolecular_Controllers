from setuptools import setup, find_packages

setup(
    name="biomolecular_controllers",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "tellurium",
        "numpy",
        "scipy",
    ],
)
