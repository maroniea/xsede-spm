from setuptools import setup

setup(name='capsolpy',
      version='0.1',
      description='Capacitance solver for an SPM tip',
      url='https://github.com/maroniea/xsede-spm',
      author='Emily Maroni and Ryan Dwyer',
      author_email='maroniea2023@mountunion.edu, dwyerry@mountunion.edu',
      license='MIT',
      packages=['capsol'],
      install_requires=[
          'tqdm',
      ],
      zip_safe=False)