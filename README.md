Simple C program to draw Julia sets for $`f_c(z)=z^2+c`$.

```shell
$ make
$ ./julia
```

The program also draws the orbit of $`0`$, the critical point of $`f_c(z)`$. As $`c`$ traces the boundary of the main cardioid ($`c=z-z^2`$) in the parameter plane, the fixed point traces the circle $`z=\frac{1}{2}e^{2\pi i\theta}`$ in the dynamic plane. Likewise, as $`c`$ traces the boundary of the period-$`2`$ bulb ($`c=-1+\frac{1}{4}e^{2\pi i\theta}`$), the fixed point pair traces the lemniscate $`|z(z+1)|=\frac{1}{4}`$, where the relation $`z^2+z+c+1=0`$ holds. Comment/uncomment the appropriate parameter setting to generate either of these cases. The parameter can also be scanned along any other desired path to see the correspoinding change in the Julia set and the orbit of the critical point.
