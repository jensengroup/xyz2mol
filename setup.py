from setuptools import setup

setup(
    name="xyz2mol",
    version="0.1.1",
    description="Convert Cartesian coordinates to one or more molecular graphs",
    url="https://github.com/jensengroup/xyz2mol",
    py_modules=["xyz2mol"],
    entry_points={"console_scripts": ["xyz2mol=xyz2mol:main"]},
)