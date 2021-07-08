from setuptools import setup

setup(
    name='capice_updatr',
    version='0.0.1',
    packages=[''],
    url='',
    license='LGPL-3.0',
    author='Robert Sietsma',
    author_email='work.robertsietsma@gmail.com',
    description='CAPICE updating service',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: LGPL-3.0',
        'Programming Language :: Python :: 3.9'
    ],
    python_requires='3.9',
    install_requires=[
        'matplotlib==3.4.2',
        'numpy==1.21.0',
        'pandas==1.3.0',
        'scipy==1.7.0',
        'xgboost==1.4.2'
    ]
)
