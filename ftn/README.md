This is the old fortran model code plus the fortran extensions for the new python model.

Manuall compilation of `ercore._flib_ercore` :

```
cd ftn/
f2py -c flib.pyf shoreline.f90 plume.f90  slipvel.f90 rnkpar.f90 interp.f90  output.f90
cd ../ercore/
ln -s ../ftn/_flib_ercore.so 
```
