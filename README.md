## stlplus

[![Build Status](https://travis-ci.org/hafen/stlplus.svg?branch=master)](https://travis-ci.org/hafen/stlplus)

![png](https://cloud.githubusercontent.com/assets/1275592/11673681/b7fce548-9dcf-11e5-8cd2-f3311b501ab9.png)

This package contains enhancements to the Seasonal Trend Decomposition using Loess (STL) implementation that comes with base R, `stl()`.

Here are some of the added features over `stl()`:

- Can handle NA values
- Higher order loess smoothing (more than just local constant and linear)
- Automated parameter choices for local quadratic
- Frequency component smoothing beyond seasonal and trend
- Plot methods for diagnostics

For (very) experimental inference, prediction, and variance reduction at endpoints, see the [operator](http://github.com/hafen/operator) package.

## References

- [Cleveland, R. B., Cleveland, W. S., McRae, J. E., & Terpenning, I. (1990). STL: A seasonal-trend decomposition procedure based on loess. *Journal of Official Statistics*, 6(1), 3-73.](http://cs.wellesley.edu/~cs315/Papers/stl%20statistical%20model.pdf)
- [Hafen, R. P. "Local regression models: Advancements, applications, and new methods." (2010).](http://ml.stat.purdue.edu/hafen/preprints/Hafen_thesis.pdf)

## Installation

```s
devtools::install_github("hafen/stlplus")
```

## License

This software is released under the BSD license.  Please read the [license](https://github.com/hafen/stlplus/blob/master/LICENSE.md) document.

