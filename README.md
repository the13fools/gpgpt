Use this starter project to accelerate your geometry processing research, yesterday. 

# 🤝 libigl 🤝 TinyAD 🤝 Polyscope 🤝 Physics thread 🤝 foliations 🤝

Welcome to the repository for the wonderful, wonderous, GPGPT!

# The General Purpose Geometry Processing Toolkit!  

This is an example project which minimally mixes a maximal number of independently useful libraries.  Put ur rings together captin planet.  

You might be confusing us with (gptoolbox)[https://github.com/alecjacobson/gptoolbox], we don't support that yet though.  

In this repo, we mix [libigl](https://github.com/libigl/libigl/) and
[TinyAD](https://github.com/patr-schm/TinyAD) and [polyscope](polyscope.run) and [a physics thread](https://github.com/evouga/libigl-example-physics-project).

On top of this, we also include our own heavily gpt generated example library for performing a geometric optimization problem and logging it's state to file for later visualization.  

If you want higher fidelity, by a path tracer like the one in blender, using scripts like this (toolbox)[https://github.com/HTDerekLiu/BlenderToolbox].

animations can then be produced using some ffmpeg scripts we also bake in.  

TODO, a solver, maybe nasoq.  Enzyme, the autodiff tool.  what else?  

# Main Attraction #1 

Mint2D.  

This is a reference implementation of the theory in the companion draft submission.  

Here we implement, blah blah blah.  


# Main Attraction #2 

Mint3D.  In which we provide a solution for recovering smooth integrable symmetric frame fields subject to various boundary conditions and optimization objectives.  

In a small step from Mint2D, we implement our primary contribution in this work.  

---

To open a nut, one can act like Serre and strike with unrelenting force, but [I] would gradually soak it to soften it, and almost without effort eventually open it.

       - Grothendieck

# Main Attraction #3 

An open source reference implementation of cube cover.  This powerful technique is broadly useful but there is not currently an avaliable reference implementation for the community to use, this project addresses this gap in the open source geometry processing software landscape.  







--------

COMING SOON!  

GPTGP-1 $20/month no additional functionality.  

--------

This example repository tries to be extremely easy to install, batteries included.  Sorry for making you download stuff each time, if you know enough to get mad about this, you can help improve the make files ;).  

--------

This repository includes a lot of gpt generated code, but its a cyborg project with lots of other snippets mixed in.  It's partially an experiment of working with a large language model.  


## Compile

Compile this project using the standard cmake routine:

    mkdir build
    cd build
    cmake ..
    make

This should find and build the dependencies and create a `bin/##` binary.

## Run

From within the `build` directory just issue:

    ./example

A glfw app should launch displaying an animating Armadillo parametrization.

![](armadillo.gif)

_Derived from
[parametrization_libigl.cc](https://github.com/patr-schm/TinyAD-Examples/blob/main/apps/parametrization_libigl.cc)_

## SUPER ANNOYING DETAIL

install suite sparse and switch out the solver that tiny ad uses like this: 

```
#pragma once

// #include <Eigen/SparseCholesky>
#include <Eigen/CholmodSupport>

namespace TinyAD
{

template <
        typename PassiveT = double,
        typename SolverT =  Eigen::CholmodSupernodalLLT< Eigen::SparseMatrix<PassiveT>> >
struct LinearSolver
{
    SolverT solver;
    bool sparsity_pattern_dirty = true;
};

}

// Eigen::SimplicialLDLT<Eigen::SparseMatrix<PassiveT>>>

```

## Ordering

TinyAD operates most conveniently on nodal vector values (e.g., vertex positions in `V`). If `V` contains nodal vectors **per row**, then regardless of whether `V` is stored as column-major (e.g., `Eigen::MatrixXd`) or row-major (e.g., `Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>`) the internal order of tiny-ad will correspond to [non-standard](https://en.wikipedia.org/wiki/Vectorization_(mathematics)) **row-major vectorization**.

That is, if
```
V = [
  x₀ y₀ z₀
  x₁ y₁ z₁
  …
  xₙ yₙ zₙ
  ]
```

Then TinyAD will vectorize this into
```
x = [
  x₀
  y₀
  z₀
  x₁
  y₁
  z₁
  …
  xₙ
  yₙ
  zₙ
  ]
```

And use a corresponding ordering for gradients and Hessians.

This is ignorable if you're using the provided `func.x_from_data` and `func.x_to_data`. However, if you're mixing in your own gradients, Hessians, constraint projections, subspace bases, then take care!


# MAC BUILD DETAILS 

fresh install details: 
brew install cmake eigen suitesparse llvm libomp

------------

Oh boy, this was quite tricky to pin down.  

Thank you eerii .  Anonomous internet comments really came through here.  

https://github.com/glfw/glfw/issues/1743#issuecomment-1229177189

For now until we improve the cmake file we added the following to our ~/.zshrc

```
export CC=$(which clang)
export CXX=$(which g++-13)
export OpenMP_ROOT=$(brew --prefix)/opt/libomp
```

BTW, todo: how to enable gdb on mac https://stackoverflow.com/a/10441587



