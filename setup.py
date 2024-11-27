import setuptools

with open("README.md", "r", encoding="utf-8") as fhand:
    long_description = fhand.read()

setuptools.setup(
    name="crystmatch",
    version="1.0.5",
    author="Fang-Cheng Wang",
    author_email="wfc@pku.edu.cn",
    description=("Package for enumerating crystal-structure matches in solid-solid phase transitions."),
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://www.phy.pku.edu.cn/xzli/RESEARCH.htm",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=["numpy", "scipy", "matplotlib", "spglib"],
    packages=setuptools.find_packages(),
    python_requires=">=3.9",
    entry_points={
        "console_scripts": [
            "crystmatch = crystmatch.cli:main",
        ]
    }
)