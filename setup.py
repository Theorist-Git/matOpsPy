from setuptools import setup, Extension
import pybind11

# ext_modules = [
#     Extension(
#         'matOpsPy',                       # module name
#         ['bindings/bindings.cpp'],        # source file(s)
#         include_dirs=[
#             pybind11.get_include(),
#             'include',                    # path to your header-only lib
#         ],
#         language='c++'
#     ),
# ]

setup(
    name="matOpsPy",
    version="0.1.0",
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
    zip_safe=False,
    classifiers=[
        'Development Status :: 3 - Alpha',
        # "3 - Alpha", "4 - Beta" or "5 - Production/Stable" as the current state of your package
        'Intended Audience :: Developers',
        'Topic :: Software Development :: Build Tools',
        'License :: OSI Approved :: MIT License',
    ],
)