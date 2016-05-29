# rcgmin
(Non-linear) Conjugate Gradient Optimizer in R with Rasmussen and More-Thuente
Line Search.

This package combines an R translation of two Matlab routines. The first is
Carl Edward Rasmussen's [conjugate gradient minimization](http://learning.eng.cam.ac.uk/carl/code/minimize/),
which has its own line search method, which uses cubic and quadratic 
interpolation and extrapolation to find a step size which satisfies the Strong
Wolfe conditions. 

Additionally, it also offers an alternative line search method: a translation of
Dianne O'Leary's translation of the [More'-Thuente line search](https://www.cs.umd.edu/users/oleary/software/)
algorithm from [MINPACK](http://www.netlib.org/minpack/).

### Installing:
```R
# install.packages("devtools")
devtools::install_github("jlmelville/rcgmin")
```

### Documentation:
```R
# The optimization function
?conj_grad

# Details on writing a custom line search function
package?rcgmin
```

### Do We Really Need Another CG Minimization Package in R?

Not really. This package was born out of an interest in the relative performance 
of  different line search methods. If you can write a steepest descent routine, 
it's not much extra work to implement conjugate gradient minimization, but an 
efficient line search becomes much more important. This turns out to be a 
non-trivial endeavour. 

Backtracking line search initially seems like a good idea: start with a 
suitably large step size, then reduce the step size until the sufficient 
decrease condition is met. Unfortunately, as far as I can tell, there's no good 
way to find the starting step size with non-linear conjugate gradient methods.

Instead, you need to allow the step size to both grow and shrink as necessary.
Doing this efficiently involves cubic or quadratic interpolation and 
extrapolation, as well as adding safeguards to ensure the resulting step
size isn't too small or large. At this point, surely, dear reader, you are
beginning to consider the prospect of writing your own line search routine as
a series of decidedly un-fun headaches. Fortunately, the two line search
routines in this package already do the hard work, and seem to have been used 
fairly extensively in their Matlab (and/or MINPACK) forms. I recommend them 
highly over writing your own.

Unlike most other optimization routines, I have tried to document and expose
the interfaces needed to implement additional line search methods, which can be
supplied as an argument to the conjugate gradient function.

The code itself remains a fairly straight translation of the two independent 
Matlab original code. It therefore may lack in consistency in terms of variable
naming, efficiency or idiomatic R usage, and redundancy.

### License
BSD 2-Clause. The original matlab code of Rasmussen is licensed this way as
part of the 
[GPML Matlab package](http://www.gaussianprocess.org/gpml/code/matlab/doc/).

The More-Thuente' line search code was translated from the Matlab code of 
Dianne O'Leary, which was in turn translated from Fortran code in MINPACK, which
has this [license](http://www.netlib.org/minpack/disclaimer). Elsewhere, this
has been treated as [BSD-like](http://mail-archives.apache.org/mod_mbox/www-legal-discuss/200609.mbox/%3C2d12b2f00609101412t7a47d99akfdb46001561d0cc4@mail.gmail.com%3E), and 
there's a translation of the Matlab code in 
[Julia](https://github.com/JuliaOpt/Optim.jl), which is MIT licensed. So the 
BSD 2-Clause is probably ok. I would welcome any correction.
