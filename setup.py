import codecs
import os
import re
import versioneer  # https://github.com/warner/python-versioneer
from setuptools import setup, find_packages

NAME = "ToFImaging"
META_PATH = os.path.join("ToFImaging", "__init__.py")
KEYWORDS = ["class", "attribute", "boilerplate"]
CLASSIFIERS = [
    "Development Status :: 5 - Production/Stable",
    "Intended Audience :: Developers",
    "Natural Language :: English",
    "License :: OSI Approved :: GPL License",
    "Operating System :: OS Independent",
    "Programming Language :: Python",
    "Programming Language :: Python :: 3",
    "Programming Language :: Python :: 3.5",
    "Programming Language :: Python :: 3.6",
    "Programming Language :: Python :: 3.7",
    "Programming Language :: Python :: 3.8",
    "Programming Language :: Python :: Implementation :: CPython",
    "Programming Language :: Python :: Implementation :: PyPy",
    "Topic :: Software Development :: Libraries :: Python Modules",
]

THIS_DIR = os.path.abspath(os.path.dirname(__file__))


def read_requirements_from_file(filepath):
    '''Read a list of requirements from the given file and split into a
    list of strings. It is assumed that the file is a flat
    list with one requirement per line.
    :param filepath: Path to the file to read
    :return: A list of strings containing the requirements
    '''
    with open(filepath, 'rU') as req_file:
        return req_file.readlines()


def read(*parts):
    """
    Build an absolute path from *parts* and and return the contents of the
    resulting file.  Assume UTF-8 encoding.
    """
    with codecs.open(os.path.join(THIS_DIR, *parts), "rb", "utf-8") as f:
        return f.read()


META_FILE = read(META_PATH)


def find_meta(meta):
    """
    Extract __*meta*__ from META_FILE.
    """
    # print (r"^__{meta}__ = ['\"]([^'\"]*)['\"]".format(meta=meta))
    # print (META_FILE)
    # print (re.M)
    meta_match = re.search(
        r"^__{meta}__ = ['\"]([^'\"]*)['\"]".format(meta=meta),
        META_FILE, re.M
    )
    if meta_match:
        return meta_match.group(1)
    raise RuntimeError("Unable to find __{meta}__ string.".format(meta=meta))


install_requires = read_requirements_from_file(os.path.join(THIS_DIR, 'requirements.txt'))
# test_requires = read_requirements_from_file(os.path.join(THIS_DIR, 'requirements_dev.txt'))


if __name__ == "__main__":
    scripts = [
    ]
    setup(
        name=NAME,
        description=find_meta("description"),
        license=find_meta("license"),
        url=find_meta("url"),
        version=versioneer.get_version(),
        cmdclass=versioneer.get_cmdclass(),
        author=find_meta("author"),
        author_email=find_meta("email"),
        maintainer=find_meta("author"),
        maintainer_email=find_meta("email"),
        keywords=KEYWORDS,
        long_description=read("README.md"),
        packages=find_packages(exclude=['test', 'test.*']),
        zip_safe=False,
        classifiers=CLASSIFIERS,
        install_requires=install_requires,
        # tests_require=test_requires,
        package_dir={},
        package_data={},
        scripts=scripts,
        setup_requires=['pytest-runner'],
    )