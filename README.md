# self_similar_public
A repository containing the code used to generate results given in "Self-similar solutions for resistive diffusion, Ohmic heating and Ettingshausen effects in plasmas of arbitrary β", found here: https://arxiv.org/abs/2111.05064

For any comments or queries, please contact gmmfarrow@gmail.com

For information on how to use these as test problems (as in our paper), please contact gmmfarrow@gmail.com

## Usage
For a casual user, the only file that should need to be edited at all is `input.py`. This changes the normalisation values, the boundary conditions and the logical switches for the different transport terms. There are also optional switches to e.g. save data or change the tolerances for the shooting method used. There is documentation within this file explaining what all of the different terms do

For a less casual user:
* `newton.py` contains the Newton-Raphson solver used in the shooting method. This is an approximate Newton-Raphson method that doesn't need an analytic derivative, rather approximating it using finite differences. Implementation was taken from Numerical Recipes.
* `transport_coeff.py` contains the coefficients for the fitting functions of the transport coefficients. Currently, using Braginskii and only Z = 1, 2 and 3 are included. Other coefficients could obviously be added.
* `unnorm.py` converts the input values for density, field and temperature into normalised values and calculates all of the dimensionless normalisation coefficients
* `gauss_pivot.py` contains a function to do Gauss-Jordan elimination with a pivoting method. Could use `scipy` or `numpy` for more efficient methods, but this works
* `equations.py` contains the actual functions that are to be evolved (e.g. the derivative functions). Any changes to e.g. the transport coefficients or the equations being solved need to be made here.
* `calculate.py` is effectively the "int main" (C++) or "program" (Fortran). Calls all of the other functions and steers the routine. Any changes to which variables should have boundary conditions enforced, or in what should be output should be made here

## Different input files

The input files for the cases produced in the paper have been provided (and named appropriately). The only input file that is documented is that of case A. The code will only run the file that is called `input.py`, so rename the files to ensure that the desired set of inputs are being used. Feel free to experiment with other inputs, but we're not guarantee that the code will run. As discussed in the paper, there are a lot of degrees of freedom and sometimes the solution will not converge.
