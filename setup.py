from setuptools import setup
import io

# read the contents of your README file
from os import path
this_directory = path.abspath(path.dirname(__file__))
with io.open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
  name = 'kinisot',
  packages = ['kinisot'],
  version = '2.0.1',
  description = 'Compute Kinetic Isotope Effects for compchem outputs',
  author = 'Paton Research Group',
  author_email = 'patonlab.colostate@edu',
  url = 'https://github.com/patonlab/kinisot',
  download_url = 'https://github.com/patonlab/Kinisot/archive/refs/tags/v2.0.1.zip',
  keywords = ['compchem', 'kie', 'dft', 'vibrational-scaling-factors', 'kinisot', 'frequencies' 'quiver'],
  classifiers = [],
  install_requires=['numpy'],
  python_requires='>=3.0',
  include_package_data=True,
)
