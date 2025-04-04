# Routh-Hurwitz Criterion Analyzer

This application evaluates the stability of a system by analyzing the polynomial provided using the Routh-Hurwitz criterion. It computes whether the system is stable or unstable, reports the number of poles with non-negative real parts, displays these poles if applicable, and prints the Routh array.

## Usage

Run the application from the command line with the following syntax:

``` bash
$ python rhc.py "<poly>"
```

- `<poly>`: The polynomial expression to analyze (enclosed in quotes).

## Example

For a polynomial:

$$s^5 + s^4 + 10s^3 + 72s^2 + 152s + 240$$

run:

``` bash
$ python rhc.py "s^5 + s^4 + 10s^3 + 72s^2 + 152s + 240"
```

or:

``` bash
$ python rhc.py "s**5 + s**4 + 10s**3 + 72s**2 + 152s + 240"
```

## Requirements

- Python 3.x
- [Sympy](https://www.sympy.org/): A Python library for symbolic mathematics.
