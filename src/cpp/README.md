# paction_cpp (SPRUCE)

## Dependencies

* [CMake](http://www.cmake.org/) (>=3.1)
* [Boost](http://www.boost.org) (>= 1.38)
* [LEMON](http://lemon.cs.elte.hu/trac/lemon) graph library (>= 1.3)

## Compilation

    $ mkdir build
    $ cd build
    $ cmake ..
    $ make

## Usage

The executable is in the build directory. The syntax for usage is as follows.

```
./enumerate <CNA tree file> <CNA proportions file> <SNV tree file> <SNV proportions file> 0 <output prefix>
```

Example of using this code on the sample files is as follows.

```
./enumerate ../../../data/sample/overview_cna_tree.csv ../../../data/sample/overview_cna.csv ../../../data/sample/overview_snv_tree.csv ../../../data/sample/overview_snv.csv 0 ../../../data/sample/cpp
```
