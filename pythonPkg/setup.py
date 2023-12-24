# A good chunk of this is from
# https://github.com/pypa/setuptools/blob/4b64119ab29b1418dc4975ea4212ab38f112c2ab/setuptools/_distutils/ccompiler.py#L964
# Linux:
#  CC=cc python3 setup.py build
#  pip install -e .

from setuptools import setup
from setuptools import Extension

descriptionStr="\
   Neeldeman Wunsch, Waterman Smith, and Hirschberg\
   pairwise alingers written in C for python. See the\
   docstrings for more details.\
";

srcFilesStr = [
      "pyAlnSeq.c"
]; # Source files to compile library

compileFlags = [
   "-O3",
   "-DDELINSSNP",
   "-static",
   "-Wno-unused-function"
]
setup(
   name = "alnSeq",
   version = "20231222",
   description = descriptionStr,
   author = "https://github.com/jeremyButtler/alnSeq",
   ext_modules=[
      Extension(
         "alnSeq",
         srcFilesStr,
         extra_compile_args = compileFlags
   )]
);
