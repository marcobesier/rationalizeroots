
read "RationalizeRoots.mpl":

# Example from section 4.1
example_1 := RationalizeRoot(sqrt(1-x^2-y^2));

# Examples from section 4.2
example_2 := ParametrizePolynomial(u^2+x^2+y^2-1);

example_3 := ParametrizePolynomial(u^2-x-y-1);

example_4 := ParametrizePolynomial(u^2-x-y-1,Variables=[u,y]);

example_5 := ParametrizePolynomial(u^2+x^2+y^2-1,OutputVariables=[v,w]);

example_6 := ParametrizePolynomial(u^2+x^2-1,GeneralC=true);

example_7 := ParametrizePolynomial(u^2-x^3-x^2,GeneralT=true);

example_8 := RationalizeRoot(sqrt((1-x1-x2-x3)^2-4*x1*x2*x3));

example_9 := RationalizeRoot(sqrt((1-x1-x2-x3)^2-4*x1*x2*x3),ForceFDecomposition=true,FPolynomials=[x1,1-x1-x2-x3,x2*x3]);

example_10 := RationalizeRoot(sqrt((1-x1-x2-x3)^2-4*x1*x2*x3),ForceFDecomposition=true,FPolynomials=[1,1-x1-x2-x3,x1*x2*x3]);

# Example from section 2.4
example_11 := RationalizeRoot(sqrt(x^4+y^3));

# Examples from section 4.3
#
# rationalise simultaneously the roots sqrt(x+1) and sqrt(x+y+1)
# rationalise sqrt(x+y+1) first
example_12_1 := RationalizeRoot(sqrt(x+y+1),OutputVariables=[v,w]);
# second roots is now sqrt(x+1)=sqrt(v*(v+w))
example_12_2 := RationalizeRoot(sqrt(v*(v+w)));
# substitute
example_12_3 := normal(subs( op(1,example_12_2), example_12_1 ));

# rationalise simultaneously the roots sqrt(1-x^2) and sqrt(1-x^2-y^2)
example_13_1 := RationalizeRoot(sqrt(1-x^2-y^2),OutputVariables=[v,w]);
# express 1-x^2 in v,w, the denomiator is a square, we only need to rationalise the numerator
example_13_2 := RationalizeRoot(sqrt(numer(normal(subs(op(1,example_13_1),1-x^2)))));
# substitute
example_13_3 := normal(subs( op(1,example_13_2), example_13_1 ));

# Examples from section 4.4
#
# for example_14_1 the program does not find a rational parametrisation
example_14_1 := RationalizeRoot(sqrt( (x^4+x^4*y+x*y^2+x^2*y^2)/x^2 ));
# leaving out the factor 1/x^2 helps
example_14_2 := RationalizeRoot(sqrt( (x^4+x^4*y+x*y^2+x^2*y^2) ));

example_15_1 := RationalizeRoot(sqrt( (x^4+4*x^2*y^2+4)/4/x^2 ));
example_15_2 := RationalizeRoot(sqrt( (x^4+4*x^2*y^2+4) ));

