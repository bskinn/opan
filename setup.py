from setuptools import setup

setup(
    name='opan',
    version='0.4.0',
    packages=['opan', 'opan.test', 'opan.utils', 'opan.vpt2'],
    package_data={'opan': ['test/resource/test.trj',
                           'test/resource/inertia/*.hess',
                           'test/resource/inertia/*.xyz',
                           'test/resource/orca/test_orca*']},
    url='https://www.github.com/bskinn/opan',
    license='The MIT License',
    author='Brian Skinn',
    author_email='bskinn@alum.mit.edu',
    description='Open Anharmonic',
    classifiers=['License :: OSI Approved :: MIT License',
                 'Natural Language :: English',
                 'Environment :: Console',
                 'Intended Audience :: Science/Research',
                 'Operating System :: OS Independent',
                 'Programming Language :: Python :: 3 :: Only',
                 'Topic :: Scientific/Engineering :: Chemistry']
)
