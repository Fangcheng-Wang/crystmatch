import os
import re
import setuptools

def parse_init_file():
    init_path = os.path.join(os.path.dirname(__file__), 'crystmatch', '__init__.py')
    variables = {}
    pattern = re.compile(r'^__(\w+)__\s*=\s*[\'"]([^\'"]+)[\'"]')
    with open(init_path, 'r', encoding='utf-8') as f:
        for line in f:
            match = pattern.match(line.strip())
            if match:
                key, value = match.groups()
                variables[key] = value
    return variables

metadata = parse_init_file()

with open("README.md", "r", encoding="utf-8") as fhand:
    long_description = fhand.read()

setuptools.setup(
    name=metadata["name"],
    version=metadata["version"],
    author=metadata["author"],
    author_email=metadata["email"],
    description=(metadata["description"]),
    long_description=long_description,
    long_description_content_type="text/markdown",
    url=metadata["url"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
    install_requires=["numpy", "scipy", "matplotlib", "spglib", "tqdm"],
    packages=setuptools.find_packages(),
    python_requires=">=3.9",
    entry_points={
        "console_scripts": [
            "crystmatch = crystmatch:main",
        ]
    }
)