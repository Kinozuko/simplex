# Simplex Project

The following project was made with C to calculate  Simplex Method. This project was made by:

* [Victor Mendoza](https://gitlab.com/Kinozuko)
* Javier Mura
* Andres Prado

## Input format

```
n # Number of coefficients
m # Number of restrictions
c1 c2 .. .. .. .. cn # Coefficients of the objective function z (vector)
a1,1 a1,2........ a1,n # Constraints of the constraint matrix (matrix nxn)
a1,1 a1,2........ a1,n
.. .. .. .. .. .. ..
.. .. .. .. .. .. ..
.. .. .. .. .. .. ..
.. .. .. .. .. .. ..
am,1 am,2 .. .. .. am,n
r1 r2 .. .. .. .. rm # Restritions >, <, = (char)
b1 b2 .. .. .. .. bm # Vector coefficients to the right of the constraints (vector)
m / M # Indicates whether to maximize (M) or minimize (m) (char)

```
## Structure

```
.
├── data/
│   └── in
├── .gitignore
├── main.c
├── makefile
├── README.md

```

## Running Project

* make
* ./main < data/in