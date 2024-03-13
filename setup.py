from setuptools import setup,find_packages

setup(name='gft',
	version='0.1',
	description='Module to compute and display the General Fourier Family Transform (GFT) and inverse.',
	url='http://github.com/robb-brown/GFT',
	author='Robert A. Brown',
	author_email='robb@shadowlabresearch.com',
	packages=find_packages(),
	install_requires=[
		'numpy',
		'scipy',
	],
	zip_safe=False)

