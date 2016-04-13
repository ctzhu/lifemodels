#from setuptools import setup

'''A setuptools based setup module, per pypa/sampleproject
'''

# Always prefer setuptools over distutils
from setuptools import find_packages
from setuptools import setup as old_setup
from numpy.distutils.misc_util import Configuration
from numpy.distutils.core import setup
# To use a consistent encoding
from codecs import open
from os import path
from os.path import join

here = path.abspath(path.dirname(__file__))
# Get the long description from the relevant file
with open(path.join(here, 'DESCRIPTION.rst'), encoding='utf-8') as f:
    long_description = f.read()


def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('lifemodels',parent_package, top_path)
    config.add_extension('mortdist', 
                         [join('lifemodels', item) for item in ['mortdist.pyf', 'mortdist.f90']])
    config.add_extension('mortdist_cpsh', 
                         [join('lifemodels', item) for item in ['mortdist_cpsh.pyf', 'mortdist_cpsh.f90']])
    return config

#config = Configuration('lifemodels', parent_package='', top_path=None)
#config.add_extension('mortdist', 
#                     [join('mortdist', item) for item in ['mortdist.pyf', 'mortdist.f90']],)
#config.add_extension('mortdist_cpsh', 
#                     [join('mortdist_cpsh', item) for item in ['mortdist_cpsh.pyf', 'mortdist_cpsh.f90']],)

#name='lifemodels' skipped as it is included in `config`
setup(version='2.02',
      description='MLE and GLM for a number of lifespan distributions commonly seen in aging research',
      long_description=long_description,
      url='http://github.com/ctzhu/lifemodels',
      author='Chen-Tseh Zhu',
      author_email='lei.ctzhu@gmail.com',
      license='MIT',

      # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
      classifiers=[
      # How mature is this project? Common values are
      # 3 - Alpha
      # 4 - Beta
      # 5 - Production/Stable
      'Development Status :: 3 - Alpha',
      # Indicate who your project is intended for
      'Intended Audience :: Developers',
      'Topic :: Software Development :: Build Tools',
      # Pick your license as you wish (should match "license" above)
      'License :: OSI Approved :: MIT License',
      # Specify the Python versions you support here. In particular, ensure
      # that you indicate whether you support Python 2, Python 3 or both.
      'Programming Language :: Python :: 2',
      'Programming Language :: Python :: 2.6',
      'Programming Language :: Python :: 2.7',
      #'Programming Language :: Python :: 3',
      #'Programming Language :: Python :: 3.2',
      #'Programming Language :: Python :: 3.3',
      #'Programming Language :: Python :: 3.4',
      ],

      keywords='MLE lifespan distribution',

      packages=['lifemodels'],
      install_requires=['matplotlib', 'numpy', 'scipy', 'pandas', 'patsy', 'pypdf'],
      zip_safe=False,
      **configuration(top_path='').todict())
