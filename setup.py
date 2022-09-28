from setuptools import setup

setup(

    name = 'pegg',
    author = 'Samuel Gould',
    author_email = 'samgould@mit.edu',
    url = 'https://github.com/samgould2/PEGG',
    version = '0.0.1',
    description = 'Prime Editing Guide Generator',
    py_modules = ["pegg"],
    package_dir = {'': 'src'},

    classifiers=[
            "Programming Language :: Python :: 3",
            "Programming Language :: Python :: 3.9",
            "License :: OSI Approved :: MIT License",
            "Operating System :: OS Independent"
        ]
)