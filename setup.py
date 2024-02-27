from setuptools import setup, find_packages


setup(
    name='HGVStoFASTA',
    version='0.4',
    license='GPL-3.0',
    author="Taner Karagol",
    author_email='taner.karagol@gmail.com',
    url='https://github.com/karagol-taner/HGVS-missense-variants-to-FASTA',
    keywords='HGVS, FASTA',
    install_requires=[
          'biopython',
      ],
    packages=['HGVStoFASTA'],
    package_dir={'HGVStoFASTA': 'HGVStoFASTA'},
    entry_points={'console_scripts': ['HGVStoFASTA = HGVStoFASTA.__main__:main']},
    include_package_data=True,

)