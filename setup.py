from setuptools import setup

setup(
    name="pyMatOps",
    version="0.1.0",
    author="Mayank Vats",
    author_email="dev-theorist.e5xna@simplelogin.com",
    description="Python bindings for matOps library",
    packages=["pyMatOps"],  # Make sure the .so file is inside pyMatOps/
    include_package_data=True,
    zip_safe=False,
)