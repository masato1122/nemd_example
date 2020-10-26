from distutils.core import setup

packages_nemd_example = [
        'nemd_example',
        'nemd_example.nemd',
        'nemd_example.tools'
        ]

setup(name='nemd_example',
      version='0.1',
      description='example for NEMD simulatin',
      author='Masato Ohnishi',
      author_email='ohnishi@photon.t.u-tokyo.ac.jp',
      packages=packages_nemd_example,
      requires=['numpy', 'ase', 'pymatgen', 'pubchempy'],
      url='https://github.com/masato1122/nemd_example.git',
      license='MIT',
      provides=['nemd_example']
      )

