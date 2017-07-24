from distutils.core import setup, Extension

# define the extension module
prob_module = Extension('prob_module', sources=['prob_module.c'])

# run the setup
setup(ext_modules=[prob_module])