# rcgmin
Conjugate Gradient Optimizer in R

This package combines an R translation of two Matlab routines: 
Carl Edward Rasmussen's [conjugate gradient minimization](http://learning.eng.cam.ac.uk/carl/code/minimize/)
and Dianne O'Leary's translation of the [More'-Thuente line search](https://www.cs.umd.edu/users/oleary/software/)
algorithm from [MINPACK](http://www.netlib.org/minpack/).

### Installing:
```R
# install.packages("devtools")
devtools::install_github("jlmelville/rcgmin")
```

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
