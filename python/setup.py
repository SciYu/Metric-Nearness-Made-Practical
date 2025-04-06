from setuptools import setup, Extension
import pybind11

ext_modules = [
    Extension(
        "matrix_optimization",
        ["matrix_optimization.cpp"],
        include_dirs=[pybind11.get_include()],
        language="c++"
    ),
]

setup(
    name="matrix_optimization",
    ext_modules=ext_modules,
    zip_safe=False,
)
