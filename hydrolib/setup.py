from setuptools import setup, find_packages

setup(
    name='hydrolib',
    version='1.0.0',
    description='Python library for hydrologic analysis including USGS data retrieval and Bulletin 17C flood frequency analysis',
    author='HydroLib',
    packages=find_packages(),
    python_requires='>=3.8',
    install_requires=[
        'numpy>=1.20.0',
        'pandas>=1.3.0',
        'matplotlib>=3.4.0',
        'scipy>=1.7.0',
        'requests>=2.25.0',
    ],
    extras_require={
        'report': ['python-docx>=0.8.11'],
    },
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Hydrology',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
    ],
)
