from setuptools import setup, find_packages

setup(name='pyiq',
      version='0.0.2',
      description='Python implement of the iq package, which is designed for protein quantification in mass spectrometry-based proteomics',
      license='MIT',
      author='Ji Hongchao',
      author_email='ji.hongchao@foxmail.com',
      url='https://github.com/hcji/pyiq',
      packages=find_packages(),
	  classifiers=[
	  'Development Status :: 4 - Beta',
	  'Programming Language :: Python :: 3.6',
	  'Programming Language :: Python :: 3.7',
      'Programming Language :: Python :: 3.8'
	  ]
     )