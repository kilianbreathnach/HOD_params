from distutils.core import setup

setup(
    name='bossvpf',
    version='0.1.0',
    description='Measure the VPF in BOSS data',
    author='Kilian Walsh',
    author_email='kilian@nyu.edu',
    license='LICENSE.txt',
    packages=['bossvpf'],
    scripts=['scripts/run-expt.py'],
    description='Measure the VPF in BOSS.',
    long_description=open('README.txt').read(),
    install_requires=[
        "scipy >= 0.8.0",
    ],
)
