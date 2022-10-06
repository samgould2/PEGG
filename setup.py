from setuptools import setup

setup(

    name = 'pegg',
    author = 'Samuel Gould',
    author_email = 'samgould@mit.edu',
    url = 'https://github.com/samgould2/PEGG',
    version = '1.0.2',
    description = 'Prime Editing Guide Generator',
    py_modules = ["pegg"],
    package_dir = {'': 'src'},

    install_requires = ["Bio>=1.4.0",
        "matplotlib>=3.5.1",
        "mock>=4.0.3",
        "numpy>=1.21.5",
        "pandas>=1.4.2",
        "seaborn>=0.11.2",
        "setuptools>=61.2.0",
        "Sphinx>=4.4.0"
    ],

    classifiers=[
            "Programming Language :: Python :: 3",
            "Programming Language :: Python :: 3.9",
            "License :: OSI Approved :: MIT License",
            "Operating System :: OS Independent"
        ]
)