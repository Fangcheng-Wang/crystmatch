import setuptools
from crystmatch import __name__, __version__, __author__, __email__, __description__, __url__

with open("README.md", "r", encoding="utf-8") as fhand:
    long_description = fhand.read()

setuptools.setup(
    name=__name__,
    version=__version__,
    author=__author__,
    author_email=__email__,
    description=(__description__),
    long_description=long_description,
    long_description_content_type="text/markdown",
    url=__url__,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
    install_requires=["numpy", "scipy", "matplotlib", "spglib"],
    packages=setuptools.find_packages(),
    python_requires=">=3.9",
    entry_points={
        "console_scripts": [
            "crystmatch = crystmatch:main",
        ]
    }
)