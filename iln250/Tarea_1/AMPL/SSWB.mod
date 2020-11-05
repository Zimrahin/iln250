#param n := 10000;
#Conjuntos
set I := {1..n}; #set

param w {I}  default Uniform(n/10, n/2);
param x {I} default trunc(Uniform(0, n+1));
param y {I} default trunc(Uniform(0, n+1));

var X;
var Y;

minimize z: sum {i in I} w[i] * sqrt((x[i]-X)^2 + (y[i] - Y)^2);
