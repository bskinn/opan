from setuptools import setup

setup(
    name='opan',
    version='0.4rc1',
    provides=['opan'],
    requires=['numpy (>=1.10)', 'scipy (>=0.12)', 'h5py (>=2.4)'],
    zip_safe=False,
    packages=['opan', 'opan.test', 'opan.utils', 'opan.vpt2'],
    package_data={'opan': ['test/resource/test.trj',
                           'test/resource/inertia/*.hess',
                           'test/resource/inertia/*.xyz',
                           'test/resource/orca/test_orca*']},
    url='https://www.github.com/bskinn/opan',
    license='The MIT License',
    author='Brian Skinn',
    author_email='bskinn@alum.mit.edu',
    description='Open Anharmonic -- A computational chemistry toolkit',
    classifiers=['License :: OSI Approved :: MIT License',
                 'Natural Language :: English',
                 'Intended Audience :: Science/Research',
                 'Operating System :: OS Independent',
                 'Programming Language :: Python :: 3 :: Only',
                 'Programming Language :: Python :: 3.4',
                 'Programming Language :: Python :: 3.5',
                 'Development Status :: 3 - Alpha',
                 'Topic :: Scientific/Engineering :: Chemistry']
)
