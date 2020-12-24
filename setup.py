import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name='pyfgwas',
    version='0.0.4',
    author='Anthony Aylward',
    author_email='aaylward@eng.ucsd.edu',
    description='Wrapper for fGWAS',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/anthony-aylward/pyfgwas.git',
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent"
    ],
    install_requires=[
        'matplotlib',
        'pandas',
        'seaborn'
    ],
    entry_points={
        'console_scripts': ['pyfgwas=pyfgwas.pyfgwas:main',]
    },
    include_package_data=True
)
