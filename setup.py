from distutils.core import setup

setup(
    name='FlowPy',
    version='0.1',
    author='Scott White',
    packages=['flowpy', 'flowpy.models'],
    description='Python library for analyzing Flow Cytometry Standard (FCS) files',
    requires=[
        'matplotlib',
        'numpy',
        'flowio',
        'flowutils'
    ]
)
