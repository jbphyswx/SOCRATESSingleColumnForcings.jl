# A simple script for updating the manifest
# files in all of our environments.
#
using Pkg: Pkg
using PkgDevTools: PkgDevTools

root = dirname(@__DIR__)

PkgDevTools.update_deps(root; auto_all = true)


