"""Install AffiTox docking pipeline package (import path: docking_benchmark2)."""

from pathlib import Path

from setuptools import find_packages, setup

readme = Path(__file__).parent / "README.md"
long_description = readme.read_text(encoding="utf-8") if readme.exists() else ""

setup(
    name="affitox",
    version="1.0.0",
    description="AffiTox: toxicity-oriented docking and binding-affinity benchmark",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    python_requires=">=3.9",
    install_requires=[
        "pandas>=1.3.0",
        "numpy<2.0",
        "matplotlib>=3.3.0",
        "seaborn>=0.11.0",
        "pyyaml>=5.4.0",
        "scipy>=1.7.0",
        "meeko>=0.5.0",
        "biopython>=1.79",
    ],
    extras_require={
        "analysis": [
            "pandas>=1.3.0",
            "numpy<2.0",
            "matplotlib>=3.3.0",
            "seaborn>=0.11.0",
            "scipy>=1.7.0",
            "pytest>=7.0",
            "rdkit>=2022.3.0",
        ],
        "rdkit": [],
    },
    entry_points={
        "console_scripts": [
            "affitox-pipeline=docking_benchmark2.cli.run_benchmark:main",
            "toxdock-pipeline=docking_benchmark2.cli.run_benchmark:main",
        ],
    },
)
