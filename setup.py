from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='geospatialroutine',
    version='0.1',
    author="Saif Aati",
    author_email="saif@caltech.edu, saifaati@gmail.com",
    description="Geospatial routines",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/SaifAati/Geospatial-Routine.git',
    python_requires='>=3.5',
     classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
    packages=['geoRoutine', 'geoRoutine.Filter', 'geoRoutine.Plotting',
              'geoRoutine.VectorField','geoRoutine.Registration'],
)
