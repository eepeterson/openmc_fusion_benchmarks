
#!/usr/bin/env python
import setuptools
import glob

setuptools.setup(
    name="benchmark_tools",
    version="0.1.0-dev",
    scripts=glob.glob('scripts/*'),
    include_package_data=True,
    description="",
    packages=setuptools.find_packages()
)
