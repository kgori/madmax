MADmax
======

Find regions of high coverage using Median Absolute Deviation (MAD)


build
=====

Requires [ldc2](https://github.com/ldc-developers/ldc/releases)

    dub build --compiler=ldc2 --build=release


usage
===

```madmax -b|--bamfile [PATH] -c|--mad-const [FLOAT] -w|--window-size [INT] -r|min-run [INT]```
```bamfile``` = path to BAM
```mad-const``` = output depths with a deviation from the median greater than ```mad-const * MAD```
```window-size``` = length of sliding window
```min-run``` = minimum number of consecutive positions to be considered a region of high depth

