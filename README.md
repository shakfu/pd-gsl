# puredata scientific library (psl)

Wraps some [gsl](https://www.gnu.org/software/gsl/) functions for puredata.

Tentatively using the name `psl` instead of `gsl` to prevent namespace collisions.

## Rationale

- Make some `gsl` functions accessible to puredata.
- See how far this can be done.
- Learn something new along the way.

## Usage

The external is a single statically linked puredata external with the name `[psl]` or `[gsl]` and which has a single inlet and a single outlet.

Basic usage is without arguments, where the object is sent messages to obtain a calculation:

```
[<func_name> <f1> <f2> .. <fN>(
```

In this case, `<func_name>` is an abbreviated name of a `gsl` function following by its arguments. The result of the calculation is sent out the outlet.

The result can be a single float value or a list of floats depending on the function.

Another usage variation is to include the name of a function in during object creation, such as `[gsl bessel_j0]`. In this case, the following rules apply:

- if the function has 1 argument, then the inlet feeds the function
- if the function has > than 1 argument then a `list` message of a length equal to its number of argument can be sent to the object as arguments:

```
[1.5 2.1 3.2(
```

Please see the file `help-psl.pd` for examples.


## To build


```
make
```

Note that the the static libraries included in this project are currently macos only. This is
a development conveniance during the early stage of this project. To make it work with other platforms just use the platform specific static libs instead.


## Development

Because quite a bit of the mapping work from `gsl` to `pd` external is quite repetitive, I have resorted to code-generation via a python script to make it more manageable.

This python3 script is in the `scripts` directory and has a dependency on the [mako](https://www.makotemplates.org) template library which is found to be quite suitable for code generation. The python3 `render.py` script uses a single `mako` template, `scripts/psl.c.mako`, and generates `psl.c` file directly in the project directory.

Effectively, this means that current development entails that changes are made to `scripts/psl.c.mako` rather than the generated `psl.c` file as per the following dev cycle:

1. modify `scripts/psl.c.mako`
2. run `scripts/render.py` to generate `psl.c`
3. make

typically running the following with each change:

```
./script/render.py && make
```

## TODO

It would be nice to have inlets which are mapped to the functions dynamically grow with the number of arguments of the funcion once it is entered as an argument.

- [ ] Revisit proxy inlets to check if it is possible to have the number of inlets dynamically grow with as per the selected functions (and its number of arguments).


