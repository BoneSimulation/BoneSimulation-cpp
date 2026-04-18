#!/bin/bash
set -e

echo "Installing required packages for Arch Linux..."

sudo pacman -Syu --noconfirm \
    base-devel \
    cmake \
    eigen \
    git \
    vtk \
    cgal \
    boost \
    libtiff \
    openmp \
    hpdf \
    gdb \
    valgrind \
    clang

echo "All dependencies installed."