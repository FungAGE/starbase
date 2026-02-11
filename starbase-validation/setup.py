"""Starbase Validation - Standalone ETL and quality gate for Django → SQLite publication."""

from setuptools import find_packages, setup

setup(
    name="starbase-validation",
    version="0.1.0",
    description="Validation and ETL app: extract from Django MySQL, validate quality, load versioned SQLite",
    author="Starbase",
    python_requires=">=3.9",
    packages=find_packages(),
    include_package_data=True,
    package_data={"starbase_validation": ["config/*.yaml"]},
    install_requires=[
        "pyyaml>=6.0",
        "sqlalchemy>=2.0",
        "pymysql>=1.0",
        "pandas>=2.0",
    ],
    extras_require={
        "dev": ["pytest>=7.0", "pytest-cov>=4.0"],
    },
    entry_points={
        "console_scripts": [
            "starbase-validate=starbase_validation.pipeline:main",
        ],
    },
)
