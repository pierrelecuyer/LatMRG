#! /usr/bin/env python
# encoding: utf-8

top = '.'
out = 'build'

from waflib import Utils
from shutil import copyfile
import os
import stat

import imp
#import sys
def waftool(name):
    return imp.load_module('waf_' + name, *imp.find_module(name, ['./waftools', './latmrg/waftools']))

version = waftool('version')
compiler = waftool('compiler')
deps = waftool('deps')

def options(ctx):
    ctx.load('compiler_c compiler_cxx gnu_dirs waf_unit_test')
    ctx.add_option('--link-static', action='store_true', help='statically link with dependencies')
    ctx.add_option('--build-docs', action='store_true', default=False, help='build documentation')
    ctx.add_option('--ntl', action='store', help='prefix under which NTL is installed')
    ctx.add_option('--gmp', action='store', help='prefix under which GMP is installed')
    ctx.add_option('--latticetester', action='store', help='prefix under which latticetester is installed')

def configure(ctx):
    print('→ prefix is ' + ctx.options.prefix)
    ctx.env.prefix = ctx.options.prefix
    build_platform = Utils.unversioned_sys_platform()
    ctx.msg("Build platform", build_platform)

    cxxflags     = ['-O3', '-fext-numeric-literals']
    CXX_EXTENSIONS = YES
    # -fext-numeric-literals
    
    ctx.load('compiler_c compiler_cxx gnu_dirs waf_unit_test')
    ctx.check(features='cxx', cxxflags='-std=c++14')
    ctx.env.append_unique('CXXFLAGS', ['-std=c++14', '-O3', '-Wall'])
    ctx.check(features='c', cflags='-std=c99')
    ctx.env.append_unique('CFLAGS', ['-std=c99', '-Wall'])
    ctx.env.INCLUDES_TEST      = ['examples'] #/usr/include



    # suppress Boost ublas warnings
    compiler.add_cxx_flag_if_supported(ctx, '-Wno-unused-local-typedefs')
    compiler.add_cxx_flag_if_supported(ctx, '-Wno-unused-function')
    compiler.add_cxx_flag_if_supported(ctx, '-Wnon-virtual-dtor')
    compiler.add_cxx_flag_if_supported(ctx, '-Wshorten-64-to-32') # clang

    if ctx.options.link_static:
        #flags = ['-static', '-static-libgcc', '-static-libstdc++']
        flags = ['-static-libgcc', '-static-libstdc++']
        if ctx.check(features='cxx cxxprogram',
                linkflags=flags,
                mandatory=False):
            ctx.env.append_unique('LINKFLAGS', flags)

    # options
    if ctx.options.ntl:
        deps.add_deps_path(ctx, 'NTL', ctx.options.ntl)
    if ctx.options.gmp:
        deps.add_deps_path(ctx, 'GMP', ctx.options.gmp)

    ctx_check = deps.shared_or_static(ctx, ctx.check)

    # NTL
    # if ctx.options.ntltypes:  # ntlttypes are now mandatory   
    ctx_check(features='cxx cxxprogram', header_name='NTL/vector.h')
    ctx_check(features='cxx cxxprogram', lib='ntl', uselib_store='NTL')

    # GMP
    ctx_check(features='cxx cxxprogram', header_name='gmp.h')
    ctx_check(features='cxx cxxprogram', lib='gmp', uselib_store='GMP')

    # Lattice Tester
    # ctx_check(features='cxx cxxprogram', header_name='../latticetester/include/latticetester/IntLattice.h')
    # ctx_check(features='cxx cxxprogram', header_name='../latticetester/include/latticetester/IntLattice.h')
    # ctx_check(features='cxx cxxprogram', header_name='/usr/local/lib/latticetester/IntLattice.h')
    # ctx_check(features='cxx cxxprogram', lib='latticetester', uselib_store='latticetester')

    # Doxygen
    if ctx.options.build_docs:
        ctx.env.BUILD_DOCS = True
        if not ctx.find_program('doxygen', var='DOXYGEN', mandatory=False):
            ctx.fatal('Doxygen is required for building documentation.\n' +
                      'Get it from http://www.stack.nl/~dimitri/doxygen/')

    ctx.version_file('latmrg')
    version_tag = ctx.set_version('latmrg')
    ctx.define('LATMRG_VERSION', version_tag)
    ctx.msg("Setting LatMRG version", version_tag)

    # definitions of DEBUG and NDEBUG flags: not used in this version
    # if not hasattr(ctx.options, 'nested') or not ctx.options.nested:
    #     # build variants
    #     env = ctx.env.derive()
    #     env.detach()

    #     # release (default)
    #     ctx.env.append_unique('CXXFLAGS', ['-O3'])
    #     ctx.define('NDEBUG', 1)

    #     ctx.setenv('debug', env)
    #     ctx.env.append_unique('CXXFLAGS', ['-g'])
    #     ctx.define('DEBUG', 1)
    

def distclean(ctx):

    verfile = ctx.path.find_node('VERSION')
    if verfile:
        verfile.delete()
    from waflib import Scripting
    Scripting.distclean(ctx)
    

def post(ctx):
    if ctx.cmd == 'install': 
        if not os.path.exists(ctx.env.prefix + '/bin/examples/data'):
           os.makedirs(ctx.env.prefix + '/bin/examples/data')
        copyfile('data/factorm', ctx.env.prefix + '/bin/examples/data/yafu')
        st = os.stat(ctx.env.prefix + '/bin/examples/data/yafu')
        os.chmod(ctx.env.prefix + '/bin/examples/data/yafu', st.st_mode | stat.S_IEXEC)


def build(ctx):

    ctx.add_group('group1')
    ctx.add_group('group2')

    if ctx.variant:
        print("Building variant `%s'" % (ctx.variant,))

    ctx.set_group('group1')
    ctx.recurse('src')
    ctx.recurse('data')    

    ctx.set_group('group2')
    ctx.recurse('examples')

    ctx.add_group('group6')
    ctx.set_group('group6')
    if not hasattr(ctx.options, 'nested') or not ctx.options.nested:
        #ctx.recurse('progs')
        if ctx.env.BUILD_DOCS:
            ctx.recurse('doc') 

    ctx.recurse('data')   
    
    ctx.add_post_fun(post)
   
    

    # ctx.recurse('analysis') # probably broken, not used in the current version


# # build variants not used in this version

# from waflib.Build import BuildContext
# class debug(BuildContext):
#     cmd = 'debug'
#     variant = 'debug'
