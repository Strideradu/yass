from distutils.core import setup, Extension

# define the extension module
prob_module = Extension('prob_module',
                        sources=['prob_module.c', 'proba.c', 'global_var.c', 'tuple.c', 'kword.c', 'util.c', 'prdyn.c'])

# run the setup
setup(ext_modules=[prob_module])
