param n := 10000;

#Conjuntos
set I := {1..n}; #set

param w{I}:= Uniform(n/10, n/2);
param x {I}:= trunc(Uniform(0, n+1));
param y {I}:= trunc(Uniform(0, n+1));

var X;
var Y;

minimize z: sum {i in I} w[i] * sqrt((x[i]-X)^2 + (y[i] - Y)^2);


#time() current time in seconds
#http://ampl.996311.n3.nabble.com/Random-generating-of-parameters-td6132.html
#usar _total_solve_system_time para obtener el tiempo total y de esta manera calcular el promedio
#usar solve_system_time para cada iteración y almacenar el máximo
