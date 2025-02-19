from setuptools import setup, Extension
import pybind11
import os

ext_modules = [
    Extension(
        "matOpsPy.matOpsPy",  # This is the module name. It will be accessible as matOpsPy in Python.
        sources=["bindings/bindings.cpp"],  # Path to your C++ source file
        include_dirs=[
            pybind11.get_include(),           # Include pybind11 headers
            os.path.join(os.path.dirname(__file__), "include")  # Include your own headers
        ],
        language="c++",
        extra_compile_args=[
            "-O3",
            "-Wall",
            "-std=c++11",
            "-fPIC",
            "-fvisibility=hidden"
        ],  # Optimization flag (adjust as needed)
    ),
]

setup(
    name="matOpsPy",
    version="0.3.0",
    author="Mayank Vats",
    author_email="dev-theorist.e5xna@simplelogin.com",
    description="Python bindings for matOps library",
    packages=["matOpsPy"],  # Make sure the .so file is inside matOpsPy/
    install_requires=[
        'pybind11>=2.13.6',
    ],
    setup_requires=[
        'pybind11>=2.13.6',
        'wheel',
    ],
    include_package_data=True,
    keywords=['Linear Algebra', 'Matrices'],
    ext_modules=ext_modules,
    zip_safe=False,
    classifiers=[
        'Development Status :: 3 - Alpha',
        # "3 - Alpha", "4 - Beta" or "5 - Production/Stable" as the current state of your package
        'Intended Audience :: Developers',
        'Topic :: Software Development :: Build Tools',
        'License :: OSI Approved :: MIT License',
    ],
)