# 🤝 libigl 🤝 TinyAD 🤝 Polyscope 🤝 Physics thread 🤝

A small example project mixing [libigl](https://github.com/libigl/libigl/) and
[TinyAD](https://github.com/patr-schm/TinyAD) and [polyscope](polyscope.run) and [a physics thread](https://github.com/evouga/libigl-example-physics-project).

Intended to generate animation traces which can be visualized by polyscope, or if you want higher fidelity, by a path tracer like the one in blender, using scripts like this (toolbox)[https://github.com/HTDerekLiu/BlenderToolbox].

animations can then be produced using some ffmpeg scripts we also bake in.  

TODO, a solver, maybe nasoq.  Enzyme, the autodiff tool.  what else?  

Use this starter project and accelerate your geometry processing research yesterday, with gpgpt.  

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

This should find and build the dependencies and create a `example` binary.

## Run

From within the `build` directory just issue:

    ./example

A glfw app should launch displaying an animating Armadillo parametrization.

![](armadillo.gif)

_Derived from
[parametrization_libigl.cc](https://github.com/patr-schm/TinyAD-Examples/blob/main/apps/parametrization_libigl.cc)_

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
