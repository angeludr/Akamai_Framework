from setuptools import setup, find_packages

NAME = 'drp'
VERSION = '1.0.0'

entry_points = {
    'console_scripts': [
        "mosaic = drp.mosaic_deimos:main"
    ]}

setup(name=NAME,
      provides=NAME,
      version=VERSION,
      packages=find_packages(),
      package_data={'kcwidrp': ['configs/*.cfg', 'data/*',
                                'data/extin/*', 'data/stds/*']},
      scripts=[],
      entry_points=entry_points,
      install_requires=['astropy~=4.0',
                        'numpy~=1.20',
                        'bokeh~=2.0.0',
                        'psutil~=5.7.0',
                        'pytest~=5.4.1',
                        # 'requests~=2.23.0',
                        'pandas~=1.0.3',
                        'matplotlib~=3.1.3',
                        'setuptools~=46.1.1',
                        'selenium',
                        # 'geckodriver',
                        'phantomjs',
                        'keckdrpframework'],
      python_requires="~=3.7"
      )