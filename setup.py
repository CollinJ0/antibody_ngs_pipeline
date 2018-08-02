import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name='antibody_ngs_pipeline',
    version='1.0.0dev',
    #packages=['antibody_ngs_pipeline',],
    description='Bulk antibody sequence preprocessing, annotaion with abstar, upload to MongoDB and S3',
    author="CollinJ0",
    url= 'https://www.github.com/CollinJ0/antibody_ngs_pipeline',
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=setuptools.find_packages(),
    scripts=['bin/antibody_ngs_pipeline'],
    classifiers=(
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ),
)
