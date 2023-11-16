# arima2 (development version)

* Fixed minor spelling errors in documentation.

# arima2 3.1.0

* Fixed bug with plotting models with only MA or only AR coefficients. 

# arima2 3.0.5

* Added `max_inv_root` and `min_inv_root_dist` arguments to `arima`:
    * `max_inv_root` controls the maximum size of the inverted AR / MA polynomial roots. Defaults to 1. 
    * `min_inv_root_dist` controls the minimum allowed distance between AR and MA polynomial roots, in order to avoid nearly canceling roots. Defaults to 0. 

# arima2 3.0.4

* Added a check to ensure a random restart is considered an improvement only if it returns a model that has a proper variance matrix for estimated coefficients. 

# arima2 3.0.3

* Added `aicc` option to `aicTable`, and removed superfluous arguments. 

# arima2 3.0.2

* Added `...` argument to `aicTable` function. 

# arima2 3.0.1

* Added arXiv paper to package description
* Added parenthesis to function names in package description. 

# arima2 3.0.0

* Added a `NEWS.md` file to track changes to the package.
