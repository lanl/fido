# Finite Difference Optimizer (fido)

An optimizer for locating stable finitie difference numerical boundary
schemes using nlopt, legion, and shoccs.  This optimizer is the code counterpart for the
boundary stencils given in the [uniform mesh](https://doi.org/10.1016/j.compfluid.2018.12.010) and [cut cell](https://doi.org/10.1016/j.jcp.2020.109794) papers .  Currently, the
documentation is woefully incomplete but there are several `.lua`
example files demonstrating how one could run the code.


# Building
The easiest way to build is with spack.  A spack repo with the recipe is supplied and can be activated via
```shell
    $ spack repo add /path/to/fido/spack-repo
```

As a quick check, a `spack repo list` should list
`/path/to/fido/spack-repo` as a package repository.  One should then
be able to simply install fido via `spack install fido`.  Rather than
directly installing the fido package, developers are encouraged to
setup a spack environment and only install the dependencies.  A
typical workflow for the dependencies  might be:

```shell
    $ spack env create fido
    $ spack env activate fido
    $ spack install --only dependencies fido
```

Followed by building fido from source:

```shell
    $ cd /path/to/fido
    $ cmake -Bbuild -H. -DCMAKE_BUILD_TYPE=...
    $ cd build && make
```


## Misc
Copyright assertion C19158
