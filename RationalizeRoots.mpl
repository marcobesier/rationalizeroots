# ---------------------------------------------------------------------------------------------
#
#  RationalizeRoots.mpl: a Software Package for the Rationalization of Square Roots
#  Version 1.0.0
#
#  Copyright (C) 2019 Stefan Weinzierl
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
# ---------------------------------------------------------------------------------------------
#
# This package defines routines for the rationalization of roots.
#
# ---------------------------------------------------------------------------------------------
#
# This package defines
#
#  RationalizeRoot(my_root, options)
#  ParametrizePolynomial(poly, options)
#
#  RationalizeRootProjective(poly,var_lst,new_var_lst, options)
#  RationalizeRootPoint(poly,var_lst,point_lst,new_var_lst, options)
#  FindPoints(poly,var_lst)
#
#  ParametrizeFDecomposition(poly,var_lst,new_var_lst, options)
#  FindFDecomposition(poly, var_lst)
#  RationalizeFDecomposition(F1,F2,F3,var_lst,new_var_lst,u0, options)
#  CheckFDegrees(poly,var_lst,F_lst)
#  GetF1F3(poly,var_lst)
#  GetDegMax(F1,F2,F3,var_lst)
#  CheckSquare(poly,var_lst)
#
#  ConvertSubsLstToPoint(var_lst,point_subs_lst, options)
#  ConvertProjectiveCoordinatesToAffineCoordinates(projective_lst)
#  ConvertAffinePolynomialToProjectivePolynomial(poly,var_lst)
#  ConvertAffinePolynomialToProjectivePolynomialOfDegree(poly,var_lst,deg)
#  SymbolFactory(base_name,lst_indices)
#  NextMultiIndex( current_index, upper_limit )
#
# ---------------------------------------------------------------------------------------------
#
# This package depends on
#
#  combinat
#
# ---------------------------------------------------------------------------------------------

with(combinat):

# ---------------------------------------------------------------------------------------------
#
# Input: my_root is of the form r1*sqrt(r2), where r1,r2 are rational functions
#
# Options: 
#
#  Variables          = []             the variables are deduced from the expression, all symbols are treated as variables
#                     = [x_1,...,x_n]  only the listed variables are treated as variables, all other occurring symbols are treated as parameters
#  OutputVariables    = []             the default choice t_1, t_2, ..., for the new variables is used
#                     = [y_1,...,y_n]  these are used as new variables
#                                      if GeneralT=true and OutputVariables non-empty, then OutputVariables should contain (n+1) variables [y_0,y_1,...,y_n]
#  MultipleSolutions  = true           multiple solutions returned
#                     = false          the first solution found is returned
#  GeneralC           = true           solutions may depend on free parameters C_1, C_2, ...
#                     = false          the value zero is substituted for all occurring free parameters
#  GeneralT           = true           rational parametrisation is given in terms of homogeneous projective coordinates
#                     = false          rational parametrisation is given in terms of affine coordinates
#  AllPoints          = true           multiple solutions returned, if multiple points of multiplicity (d-1) were found
#                     = false          the first solution found is returned
#  AllCharts          = true           all charts in projective space are searched for solutions
#                     = false          only the default chart is searched for solutions
#  TryFDecomposition  = true           if no rational parametrisation is found with the standard method, try in addition the FDecomposition function
#                     = false          do not use the FDecomposition function
#  AllFDecomposition  = true           multiple solutions from FDecomposition returned
#                     = false          the first solution found is returned
#  FPolynomials       = []             use a heuristic algorithm to find possible F-decompositions
#                     = [F1,F2,F3]     use this F-decomposition, where sqrt(poly)=sqrt(F2^2-4*F1*F3)
#  ForceFDecomposition= true           do not use the standard algorithm, use only the F-decomposition algorithm
#                       false          use first the standard algorithm
#
# Output: A a list of rational parametrisations as substitution lists [ [x_1 = f_1(y_1,...,y_n), ..., x_n = f_n(y_1,...,y_n)], ... ].
#         If MultipleSolutions=false the list contains one substitution list.
#         If no rational parametrisation is found, the empty list is returned.
#
RationalizeRoot := proc( my_root, { Variables::list              := [], 
                                    OutputVariables::list        := [], 
                                    MultipleSolutions::boolean   := false, 
                                    GeneralC::boolean            := false, 
                                    GeneralT::boolean            := false, 
                                    AllPoints::boolean           := true, 
                                    AllCharts::boolean           := true, 
                                    TryFDecomposition::boolean   := true, 
                                    AllFDecomposition::boolean   := true,
                                    FPolynomials::list           := [],
                                    ForceFDecomposition::boolean := false })

local my_u,poly,my_vars, my_sol;

# generate symbol for the root
my_u := SymbolFactory('_u',0);

# if Variables is non-empty, add new variable,
# otherwise the variables are figured out later
my_vars := Variables;
if ( my_vars <> [] ) then
 my_vars := [ my_u, op(my_vars) ];
end if;

# get the (affine) polynomial
poly := numer(normal(my_u^2 - my_root^2));

my_sol := ParametrizePolynomial(poly,  ':-Variables'           = my_vars, 
                                       ':-OutputVariables'     = OutputVariables, 
                                       ':-MultipleSolutions'   = MultipleSolutions, 
                                       ':-GeneralC'            = GeneralC, 
                                       ':-GeneralT'            = GeneralT, 
                                       ':-AllPoints'           = AllPoints, 
                                       ':-AllCharts'           = AllCharts, 
                                       ':-TryFDecomposition'   = TryFDecomposition, 
                                       ':-AllFDecomposition'   = AllFDecomposition,
                                       ':-FPolynomials'        = FPolynomials,
                                       ':-ForceFDecomposition' = ForceFDecomposition,
                                       ':-IsProjective'        = false, 
                                       ':-OutputRoot'          = false );

return my_sol;

end proc:

# ---------------------------------------------------------------------------------------------
#
# Input: poly is a (multi-variate) polynomial
#
# Options: 
#
#  Variables          = []             the variables are deduced from the expression, all symbols are treated as variables
#                     = [x_1,...,x_n]  only the listed variables are treated as variables, all other occurring symbols are treated as parameters
#  OutputVariables    = []             the default choice t_1, t_2, ..., for the new variables is used
#                     = [y_1,...,y_m]  these are used as new variables
#  MultipleSolutions  = true           multiple solutions returned
#                     = false          the first solution found is returned
#  GeneralC           = true           solutions may depend on free parameters C_1, C_2, ...
#                     = false          the value zero is substituted for all occurring free parameters
#  GeneralT           = true           rational parametrisation is given in terms of homogeneous projective coordinates
#                     = false          rational parametrisation is given in terms of affine coordinates
#  AllPoints          = true           multiple solutions returned, if multiple points of multiplicity (d-1) were found
#                     = false          the first solution found is returned
#  AllCharts          = true           all charts in projective space are searched for solutions
#                     = false          only the default chart is searched for solutions
#  TryFDecomposition  = true           if no rational parametrisation is found with the standard method, try in addition the FDecomposition function
#                     = false          do not use the FDecomposition function
#  AllFDecomposition  = true           multiple solutions from FDecomposition returned
#                     = false          the first solution found is returned
#  FPolynomials       = []             use a heuristic algorithm to find possible F-decompositions
#                     = [F1,F2,F3]     use this F-decomposition, where poly=u^2-F2^2+4*F1*F3
#  ForceFDecomposition= true           do not use the standard algorithm, use only the F-decomposition algorithm
#                       false          use first the standard algorithm
#  IsProjective       = true           poly is a polynomial in projective space
#                     = false          poly is a polynomial in affine space
#  OutputRoot         = true           output contains parametrisations for all variables
#                     = false          first variable is the named root, suppress output for this variable
#
# Output: A list of rational parametrisations as substitution lists [ [x_1 = f_1(y_1,...,y_m), ..., x_n = f_n(y_1,...,y_m)], ... ].
#         If MultipleSolutions=false the list contains one substitution list.
#         If no rational parametrisation is found, the empty list is returned.
#
ParametrizePolynomial := proc( poly, { Variables::list              := [], 
                                       OutputVariables::list        := [], 
                                       MultipleSolutions::boolean   := false, 
                                       GeneralC::boolean            := false, 
                                       GeneralT::boolean            := false, 
                                       AllPoints::boolean           := true, 
                                       AllCharts::boolean           := true, 
                                       TryFDecomposition::boolean   := true, 
                                       AllFDecomposition::boolean   := true, 
                                       FPolynomials::list           := [],
                                       ForceFDecomposition::boolean := false,
                                       IsProjective::boolean        := false, 
                                       OutputRoot::boolean          := true })

local my_vars, my_poly, my_u, n_var, my_new_vars, i1, my_sol, my_u2;

my_vars := Variables;
if ( my_vars = [] ) then
 my_vars := indets(poly);
end if;

if ( IsProjective = true ) then

 my_poly := poly;

else # affine polynomial

 # generate symbol for projective variable
 my_u := SymbolFactory('_u',1);

 # add additional variable to the end
 my_vars := [ op(my_vars), my_u ];

 my_poly := ConvertAffinePolynomialToProjectivePolynomial(poly,my_vars);

end if;

n_var := nops(my_vars);

my_new_vars := OutputVariables;
if ( OutputVariables = [] ) then
 for i1 from 1 to (n_var-2) do
  my_new_vars := [ op(my_new_vars), SymbolFactory('t',i1) ];
 end do;

 # add one variable
 if ( GeneralT = true ) then
  my_new_vars := [ SymbolFactory('t',0), op(my_new_vars)];
 end if;

end if;

my_sol := [];
if ( ForceFDecomposition=false ) then
 my_sol := RationalizeRootProjective(my_poly,my_vars,my_new_vars,
                                     ':-MultipleSolutions' = MultipleSolutions, 
                                     ':-GeneralC'          = GeneralC, 
                                     ':-GeneralT'          = GeneralT, 
                                     ':-AllPoints'         = AllPoints, 
                                     ':-AllCharts'         = AllCharts, 
                                     ':-IsProjective'      = IsProjective, 
                                     ':-OutputRoot'        = OutputRoot );
end if;

# try FDecomposition 
if ( (my_sol = []) and ( TryFDecomposition= true ) and ( IsProjective = false ) ) then
 
 # generate additional symbol
 my_u2 := SymbolFactory('_u',2);
 my_vars := [ op(1..(nops(my_vars)-1),my_vars), my_u2, op(nops(my_vars),my_vars) ];

 my_sol := ParametrizeFDecomposition(poly,my_vars,my_new_vars,
                                     ':-MultipleSolutions'  = MultipleSolutions, 
                                     ':-GeneralC'           = GeneralC, 
                                     ':-GeneralT'           = GeneralT, 
                                     ':-AllPoints'          = AllPoints, 
                                     ':-AllCharts'          = AllCharts, 
                                     ':-AllFDecomposition'  = AllFDecomposition,
                                     ':-FPolynomials'       = FPolynomials,
                                     ':-OutputRoot'         = OutputRoot );

end if;

return my_sol;

end proc:

# ---------------------------------------------------------------------------------------------
#
# Input: poly is a homogeneous polynomial in the variables given by the list var_lst = [x1,x2,...,xn], i.e. in projective space.
#        new_var_lst = [t1,...,t_{n-2}] is a list of (n-2) variables parametrising the rationalisation.
#
# Options: 
#
#  MultipleSolutions = true           multiple solutions returned
#                    = false          the first solution found is returned
#  GeneralC          = true           solutions may depend on free parameters C_1, C_2, ...
#                    = false          the value zero is substituted for all occurring free parameters
#  GeneralT          = true           rational parametrisation is given in terms of homogeneous projective coordinates
#                    = false          rational parametrisation is given in terms of affine coordinates
#  AllPoints         = true           multiple solutions returned, if multiple points of multiplicity (d-1) were found
#                    = false          the first solution found is returned
#  AllCharts         = true           all charts in projective space are searched for solutions
#                    = false          only the default chart is searched for solutions
#  IsProjective      = true           poly is a polynomial in projective space
#                    = false          poly is a polynomial in affine space
#  OutputRoot        = true           output contains parametrisations for all variables
#                    = false          first variable is the named root, suppress output for this variable
#
# Output: A list of rational parametrisations as substitution lists [ [x_1 = f_1(t_1,...,t_{n-2}), ..., x_n = f_n(t_1,...,t_{n-2})], ... ].
#         If MultipleSolutions=false the list contains one substitution list.
#         If no rational parametrisation is found, the empty list is returned.
#
# if IsProjective    = true  parametrisations for x_1,...,x_n are returned (projective case)
#                    = false parametrisations for x_1,...,x_{n-1} are returned (x_n=1, affine case)
#
# if OutputRoot      = true  (and IsProjective = false) parametrisations for x_1,x_2,x_3 ..., x_{n-1} are returned
#                    = false (and IsProjective = false) parametrisations for     x_2,x-3 ..., x_{n-1} are returned
#
RationalizeRootProjective := proc( poly, var_lst, new_var_lst,
                                   { MultipleSolutions::boolean := false, 
                                     GeneralC::boolean          := false, 
                                     GeneralT::boolean          := false, 
                                     AllPoints::boolean         := true, 
                                     AllCharts::boolean         := true, 
                                     IsProjective::boolean      := false, 
                                     OutputRoot::boolean        := true })

local n_var,my_sol,i1,j1,k1,my_poly,my_var_lst,my_subs_points,my_point,my_temp_sol,my_temp_proj_sol,my_temp_affine_sol,my_affine_sol,my_proj_sol;

n_var := nops(var_lst);

my_sol := [];

# loop over charts
for i1 from n_var by -1 to 1 do
 if ( AllCharts or (i1=n_var) ) then

  my_poly := subs( op(i1,var_lst)=1, poly);
  my_var_lst := [];
  for j1 from 1 to n_var do
   if ( j1 <> i1 ) then
    my_var_lst := [op(my_var_lst), op(j1,var_lst) ];
   end if;
  end do;

  # find points of multiplicity (deg-1)
  my_subs_points := FindPoints(my_poly,my_var_lst);

  for j1 from 1 to nops(my_subs_points) do
   if ( AllPoints or (j1=1) ) then

    my_point := ConvertSubsLstToPoint(my_var_lst,op(j1,my_subs_points), ':-GeneralC' = GeneralC);
    if ( my_point <> [] ) then
     my_temp_sol := RationalizeRootPoint(my_poly,my_var_lst,my_point,new_var_lst, ':-GeneralT' = GeneralT);
     my_temp_proj_sol := [];
     for k1 from 1 to n_var do
      if ( k1=i1 ) then
       my_temp_proj_sol := [ op(my_temp_proj_sol), 1 ];
      elif (k1<i1) then
       my_temp_proj_sol := [ op(my_temp_proj_sol), op(k1,my_temp_sol) ];
      else
       my_temp_proj_sol := [ op(my_temp_proj_sol), op(k1-1,my_temp_sol) ];
      end if;
     end do;

     # output
     if ( IsProjective = false ) then # affine polynomial

      if ( normal(op(n_var,my_temp_proj_sol)) <> 0 ) then
       my_temp_affine_sol := ConvertProjectiveCoordinatesToAffineCoordinates(my_temp_proj_sol);

       my_affine_sol := [];
       for k1 from 1 to (n_var-1) do
        my_affine_sol := [ op(my_affine_sol), op(k1,var_lst) = op(k1,my_temp_affine_sol) ];
       end do;

       if ( OutputRoot ) then

        if ( MultipleSolutions = false ) then
         return [my_affine_sol];
        end if;
        my_sol := [ op(my_sol), my_affine_sol ];

       else

        if ( MultipleSolutions = false ) then
         return [[ op(2..nops(my_affine_sol),my_affine_sol) ]];
        end if;
        my_sol := [ op(my_sol), [ op(2..nops(my_affine_sol),my_affine_sol) ] ];

       end if; # ( OutputRoot )
      end if; # ( normal(op(n_var,my_temp_proj_sol)) <> 0 )

     else # ( IsProjective = true ), projective polynomial

       my_proj_sol := [];
       for k1 from 1 to n_var do
        my_proj_sol := [ op(my_proj_sol), op(k1,var_lst) = op(k1,my_temp_proj_sol) ];
       end do;

      if ( MultipleSolutions = false ) then
       return [my_proj_sol];
      end if;
      my_sol := [ op(my_sol), my_proj_sol ];

     end if; # IsProjective

    end if; # ( my_point <> [] )
   end if; # ( AllPoints or (j1=1) )
  end do; # j1, loop over points

 end if; # ( AllCharts or (i1=n_var) )
end do; # i1, loop over charts

# remove identical entries
my_sol := [op({op(my_sol)})];

return my_sol;

end proc:

# ---------------------------------------------------------------------------------------------
#
# Input: poly is a polynomial in the variables given by the list var_lst=[x1,...,xn].
#        point_lst = [a1,...,an] is a point of multiplicity of (deg-1).
#        new_var_lst = [t2,...,t_n] is a list of (n-1) variables parametrising the rationalisation.
#                    if  GeneralT = true, new_var_lst = [t1,t2,...,t_n] is a list of n homogeneous variables
#
# Options: 
#
#  GeneralT = true           rational parametrisation is given in terms of homogeneous projective coordinates
#           = false          rational parametrisation is given in terms of affine coordinates
#
# Output: [phi1,...,phin], a rational parametrisation of the hypersurface according to eq.(54) of arXiv:1809.10983
#
RationalizeRootPoint := proc(poly,var_lst,point_lst,new_var_lst, {GeneralT::boolean := false})

local n_var,deg,temp_var_lst,subs_lst,i1,ll,g,gd,gdm1,sol_lst;

n_var := nops(var_lst);
deg := degree( expand(poly), {op(var_lst)});

if ( GeneralT = true ) then
 temp_var_lst := new_var_lst;
else
 temp_var_lst := [ 1, op(new_var_lst) ];
end if;

subs_lst := [];
for i1 from 1 to n_var do
 subs_lst := [ op(subs_lst), op(i1,var_lst) = ll*op(i1,temp_var_lst) + op(i1,point_lst) ];
end do;

g    := collect(subs(subs_lst,poly),ll);
gd   := coeff(g,ll,deg);
gdm1 := coeff(g,ll,deg-1);

sol_lst := [];

for i1 from 1 to n_var do
 sol_lst := [ op(sol_lst), -op(i1,temp_var_lst)*gdm1/gd + op(i1,point_lst) ];
end do;

return sol_lst;

end proc:

# ---------------------------------------------------------------------------------------------
#
# Input: poly is a polynomial in the variables given by the list var_lst.
#
# Output: A list of points of multiplicity (deg-1)
#
# Remark: We may assume that deg(poly)>=2.
#
FindPoints := proc(poly,var_lst)

local n_var,deg,eq_lst,my_degree,composition_lst,i1,temp_sublst,j1,diff_lst,my_sol_temp,my_sol;

n_var := nops(var_lst);
deg := degree( expand(poly), {op(var_lst)});

if (deg<2) then
 error "polynomial has degree %1", deg;
end if;

# set up equations up to degree (deg-2)
eq_lst := [];
for my_degree from 0 to (deg-2) do
 composition_lst := composition(my_degree+n_var, n_var);
 for i1 from 1 to nops(composition_lst) do
  temp_sublst := op(i1,composition_lst);
  diff_lst := [];
  for j1 from 1 to n_var do
   diff_lst := [ op(diff_lst), op(j1,var_lst)$(op(j1,temp_sublst)-1) ];
  end do;
  eq_lst := [ op(eq_lst), diff(poly,diff_lst) ];
 end do;
end do;

my_sol_temp := solve( eq_lst, var_lst);

# remove points of higher multiplicity
eq_lst := [];
my_degree := deg-1;
composition_lst := composition(my_degree+n_var, n_var);
for i1 from 1 to nops(composition_lst) do
 temp_sublst := op(i1,composition_lst);
 diff_lst := [];
 for j1 from 1 to n_var do
  diff_lst := [ op(diff_lst), op(j1,var_lst)$(op(j1,temp_sublst)-1) ];
 end do;
 eq_lst := [ op(eq_lst), diff(poly,diff_lst) ];
end do;

my_sol := [];
for i1 from 1 to nops(my_sol_temp) do
 if ( map(normal,subs(op(i1,my_sol_temp),eq_lst)) <> [seq(0, j1=1..n_var)] ) then
  my_sol := [ op(my_sol), op(i1,my_sol_temp) ];
 end if;
end do:

return my_sol;

end proc:

# ---------------------------------------------------------------------------------------------
#
# Rationalise u^2 - F2^2 + 4*F1*F3
# A rational parametrisation of F1+F2+F3 can be used to obtain a rational parametrisation of u^2-F2^2+4*F1*F3
#
# Options: 
#
#  MultipleSolutions  = true           multiple solutions returned
#                     = false          the first solution found is returned
#  GeneralC           = true           solutions may depend on free parameters C_1, C_2, ...
#                     = false          the value zero is substituted for all occurring free parameters
#  GeneralT           = true           rational parametrisation is given in terms of homogeneous projective coordinates
#                     = false          rational parametrisation is given in terms of affine coordinates
#  AllPoints          = true           multiple solutions returned, if multiple points of multiplicity (d-1) were found
#                     = false          the first solution found is returned
#  AllCharts          = true           all charts in projective space are searched for solutions
#                     = false          only the default chart is searched for solutions
#  AllFDecomposition  = true           multiple solutions from FDecomposition returned
#                     = false          the first solution found is returned
#  FPolynomials       = []             use a heuristic algorithm to find possible F-decompositions
#                     = [F1,F2,F3]     use this F-decomposition, where poly=u^2-F2^2+4*F1*F3
#  OutputRoot         = true           output contains parametrisations for all variables
#                     = false          first variable is the named root, suppress output for this variable
#
# Output: A list of rational parametrisations as substitution lists [ [x_1 = f_1(y_1,...,y_m), ..., x_n = f_n(y_1,...,y_m)], ... ].
#         If MultipleSolutions=false the list contains one substitution list.
#         If no rational parametrisation is found, the empty list is returned.
#
ParametrizeFDecomposition := proc(poly,var_lst,new_var_lst, 
                                  { MultipleSolutions::boolean := false, 
                                    GeneralC::boolean          := false, 
                                    GeneralT::boolean          := false, 
                                    AllPoints::boolean         := true, 
                                    AllCharts::boolean         := true, 
                                    AllFDecomposition::boolean := true, 
                                    FPolynomials::list         := [],
                                    OutputRoot::boolean        := true })

local n_var,kk,my_vars,my_vars1,my_vars2,my_poly,F_lst, i1,j1, my_sol, my_F1,my_F2,my_F3, my_temp_sol, deg_max;

n_var := nops(var_lst);

my_sol := [];

for kk from 1 to (n_var-2) do
 # check for variable u^2
 if ( (degree(expand(poly),{op(kk,var_lst)})=2) and (diff(poly,(op(kk,var_lst))$2)=2) and (coeff(collect(expand(poly),op(kk,var_lst)),op(kk,var_lst),1)=0) ) then

  # original polynomial is u^2 - F2^2 + 4*F1*F3, convert to F2^2 - 4*F1*F3
  my_poly  := subs( op(kk,var_lst)=0, -poly );
  my_vars  := [ op(1..(kk-1),var_lst), op((kk+1)..n_var,var_lst) ];
  my_vars1 := [op(1..(nops(my_vars)-1),my_vars)];
  my_vars2 := [op(1..(nops(my_vars)-2),my_vars)];

  if ( FPolynomials = [] ) then
   F_lst := FindFDecomposition(my_poly, my_vars1 );
  else
   my_F1 := op(1,FPolynomials);
   my_F2 := op(2,FPolynomials);
   my_F3 := op(3,FPolynomials);

   if ( expand( my_poly - my_F2^2 + 4*my_F1*my_F3 ) = 0 ) then
    deg_max := GetDegMax(my_F1,my_F2,my_F3,my_vars2);

    F_lst := [[ConvertAffinePolynomialToProjectivePolynomialOfDegree(my_F1,my_vars1,deg_max-2),
               ConvertAffinePolynomialToProjectivePolynomialOfDegree(my_F2,my_vars1,deg_max-1),
               ConvertAffinePolynomialToProjectivePolynomialOfDegree(my_F3,my_vars1,deg_max)]];
   else # wrong kk, sanity check failed
    F_lst := [];
   end if;

  end if; # not( FPolynomials = [] )

  for i1 from 1 to nops(F_lst) do
   if ( AllFDecomposition or (i1=1) ) then

    my_F1 := op(1,op(i1,F_lst));
    my_F2 := op(2,op(i1,F_lst));
    my_F3 := op(3,op(i1,F_lst));

    my_temp_sol := RationalizeFDecomposition(my_F1,my_F2,my_F3,my_vars,new_var_lst,op(kk,var_lst),
                                             ':-MultipleSolutions' = MultipleSolutions, 
                                             ':-GeneralC'          = GeneralC, 
                                             ':-GeneralT'          = GeneralT, 
                                             ':-AllPoints'         = AllPoints, 
                                             ':-AllCharts'         = AllCharts, 
                                             ':-OutputRoot'        = OutputRoot );

    if ( (MultipleSolutions=false) and (my_temp_sol <> [] ) ) then
     return my_temp_sol;
    end if;

    for j1 from 1 to nops(my_temp_sol) do
     my_sol := [ op(my_sol), op(j1,my_temp_sol) ]; 
    end do; # j1 

   end if; # ( AllFDecomposition or (i1=1) )
  end do; # i1

 end if; # check for u^2
end do; # kk

# remove identical entries
my_sol := [op({op(my_sol)})];

return my_sol;

end proc:

# ---------------------------------------------------------------------------------------------
#
# Input: poly is a polynomial in the first (n-1) variables given by the list var_lst = [x1,x2,...,xn].
#
# Output: A list of decompositions [ [F1_1,F2_1,F3_1], [F1_2,F2_2,F3_2], ...]
#         such that poly = F2^2 - 4*F1*F3 for xn=1
#         F1 is homogenous in [x1,x2,...,xn] of degree deg_max-2
#         F2 is homogenous in [x1,x2,...,xn] of degree deg_max-1
#         F3 is homogenous in [x1,x2,...,xn] of degree deg_max
#         In addition deg_max < deg(poly)
#
FindFDecomposition := proc(poly, var_lst)

local i1,deg,F_lst,my_res,my_F1,my_F2,my_F3,deg_max;

deg := degree( expand(poly), {op(var_lst)});

F_lst := GetF1F2F3(poly,var_lst);

my_res := [];

for i1 from 1 to nops(F_lst) do

 my_F1 := op(1,op(i1,F_lst));
 my_F2 := op(2,op(i1,F_lst));
 my_F3 := op(3,op(i1,F_lst));

 deg_max := GetDegMax(my_F1,my_F2,my_F3,var_lst);

 # require lower degree
 if ( deg_max < deg ) then
  my_res := [ op(my_res), [ ConvertAffinePolynomialToProjectivePolynomialOfDegree(my_F1,var_lst,deg_max-2),
                            ConvertAffinePolynomialToProjectivePolynomialOfDegree(my_F2,var_lst,deg_max-1),
                            ConvertAffinePolynomialToProjectivePolynomialOfDegree(my_F3,var_lst,deg_max)] ]; 
 end if;

end do;

return my_res;

end proc:

# ---------------------------------------------------------------------------------------------
#
# Input: F1, F2, F3 are homogeneous polynomial in the first (n-1) variables given by the list var_lst = [x1,x2,...,xn].
#        deg(F1) = deg_max-2
#        deg(F2) = deg_max-1
#        deg(F3) = deg_max
#        new_var_lst = [t1,...,t_{n-2}] is a list of (n-2) variables parametrising the rationalisation.
#
# Options: 
#
#  MultipleSolutions = true           multiple solutions returned
#                    = false          the first solution found is returned
#  GeneralC          = true           solutions may depend on free parameters C_1, C_2, ...
#                    = false          the value zero is substituted for all occurring free parameters
#  GeneralT          = true           rational parametrisation is given in terms of homogeneous projective coordinates
#                    = false          rational parametrisation is given in terms of affine coordinates
#  AllPoints         = true           multiple solutions returned, if multiple points of multiplicity (d-1) were found
#                    = false          the first solution found is returned
#  AllCharts         = true           all charts in projective space are searched for solutions
#                    = false          only the default chart is searched for solutions
#  OutputRoot        = true           output contains parametrisations for all variables
#                    = false          first variable is the named root, suppress output for this variable
#
RationalizeFDecomposition := proc(F1,F2,F3,var_lst,new_var_lst,u0,
                                  { MultipleSolutions::boolean := false, 
                                    GeneralC::boolean          := false, 
                                    GeneralT::boolean          := false, 
                                    AllPoints::boolean         := true, 
                                    AllCharts::boolean         := true,
                                    OutputRoot::boolean        := true })

local n_var,my_F,i1,j1,my_sol,my_W_sol,my_W_expr,my_V_expr,my_subs_lst;

n_var := nops(var_lst);

my_F := ConvertAffinePolynomialToProjectivePolynomial( F1+F2+F3, var_lst );

my_W_sol := RationalizeRootProjective(my_F,var_lst,new_var_lst,
                                      ':-MultipleSolutions' = MultipleSolutions, 
                                      ':-GeneralC'          = GeneralC, 
                                      ':-GeneralT'          = GeneralT, 
                                      ':-AllPoints'         = AllPoints, 
                                      ':-AllCharts'         = AllCharts, 
                                      ':-IsProjective'      = false, 
                                      ':-OutputRoot'        = true );

if ( my_W_sol = [] ) then
 return [];
end if;

my_sol := [];

for i1 from 1 to nops(my_W_sol) do

 my_W_expr := op(i1,my_W_sol);

 my_subs_lst := [];
 for j1 from 1 to (n_var-2) do
  my_subs_lst := [ op(my_subs_lst), op(j1,var_lst) = normal( subs(my_W_expr, op(j1,var_lst)/op(n_var-1,var_lst)) ) ];
 end do;
 my_subs_lst := [ op(my_subs_lst), op(n_var-1,var_lst) = 1 ];

 if ( OutputRoot ) then
  my_V_expr := [ u0 = normal( 2 * subs(my_subs_lst,F3) * subs(my_W_expr,op(n_var-1,var_lst)) + subs(my_subs_lst,F2) ) ];
 else
  my_V_expr := [];
 end if;

 my_V_expr := [ op(my_V_expr), op(1..(n_var-2),my_subs_lst) ]; 

 my_sol := [ op(my_sol), my_V_expr ];

 if ( MultipleSolutions=false ) then
  return [my_V_expr];
 end if;

end do; # i1

return my_sol;

end proc:

# ---------------------------------------------------------------------------------------------
#
# Input: poly is a polynomial in the variables given by var_lst
#
# Output: A list [ [F1_1,F2_1,F3_1], [F1_2,F2_2,F3_2], ...]
# such that
#  P = F2^2 - 4*F1*F3
#
GetF1F2F3 := proc(poly,var_lst)

local i1,my_poly,n_terms,my_seq,my_power_seq,my_res,my_subset_seq,part1,part2,sqrt_part1,my_F1F3_fact;

my_poly := expand(poly);

n_terms := nops(my_poly);

my_seq := [seq(i1,i1=1..n_terms)];

my_res := [];

my_power_seq := subsets(my_seq);
while ( not(my_power_seq[finished]) ) do

 my_subset_seq := my_power_seq[nextvalue]();
 part1 := 0;
 for i1 from 1 to nops(my_subset_seq) do
  part1 := part1 + op(op(i1,my_subset_seq),my_poly);
 end do;
 part2 := my_poly - part1;

 sqrt_part1 := CheckSquare(part1,var_lst);

 if ( (sqrt_part1 <> 0) or (part1 = 0) ) then

  my_F1F3_fact := GetF1F3(part2,var_lst);
  for i1 from 1 to nops(my_F1F3_fact) do
   my_res := [ op(my_res), [op(1,op(i1,my_F1F3_fact)), sqrt_part1, op(2,op(i1,my_F1F3_fact))] ]; 
  end do;

 end if; # ( (sqrt_part1 <> 0) or (part1 = 0) )

end do; # while ( not(my_power_seq[finished]) )


return my_res;

end proc:

# ---------------------------------------------------------------------------------------------
#
# Input: poly is a polynomial in the variables given by var_lst
#
# The routine returns a list of possible factorisations [ [F1_1,F3_1], [F1_1,F3_1], ...] (including trivial factorisations like F1=1)
# such that P = - 4*F1*F3
#
GetF1F3 := proc(poly,var_lst)

local my_content,my_poly,my_factor_lst,my_exp_lst,i1,j1,k,my_index,flag_counter,temp,res_lst,my_F1,my_F3;


my_content := content(poly,var_lst);

my_poly := factor(poly/my_content);

if ( type(my_poly,`*`) ) then
 my_factor_lst := [];
 my_exp_lst := [];
 for i1 from 1 to nops(my_poly) do 
  if ( type(op(i1,my_poly),`^`) ) then 
   my_factor_lst := [ op(my_factor_lst), op(1,op(i1,my_poly)) ];
   my_exp_lst := [ op(my_exp_lst), op(2,op(i1,my_poly)) ];
  else
   my_factor_lst := [ op(my_factor_lst), op(i1,my_poly) ];
   my_exp_lst := [ op(my_exp_lst), 1 ];
  end if
 end do;
elif ( type(my_poly,`*`) ) then
 my_factor_lst := [op(1,my_poly)];
 my_exp_lst := [op(2,my_poly)];
else
 my_factor_lst := [my_poly];
 my_exp_lst := [1];
end if;

k := nops(my_factor_lst);
my_index := [ seq(0,i1=1..k) ];
flag_counter := true;

res_lst := [];

while (flag_counter) do

 my_F1 := 1;
 my_F3 := 1;

 for i1 from 1 to k do
  my_F1 := my_F1 * ( op(i1,my_factor_lst) )^(op(i1,my_index));
  my_F3 := my_F3 * ( op(i1,my_factor_lst) )^(op(i1,my_exp_lst)-op(i1,my_index));
 end do;

 res_lst := [ op(res_lst), [my_content/(-4)*my_F1, my_F3] ];

 temp := NextMultiIndex( my_index, my_exp_lst );
 my_index := op(2,temp); 
 flag_counter := not(op(1,temp));
end do;

return res_lst;

end proc:

# ----------------------------------------------------------------------------
#
# F1,F2,F3 are polynomials in var_lst.
#
# The routine returns
#  max( deg(F3), deg(F2)+1, deg(F1)+2 )
#
GetDegMax := proc(F1,F2,F3,var_lst)

local deg_F1,deg_F2,deg_F3,deg_max;

deg_F1 := degree( expand(F1), {op(var_lst)});
deg_F2 := degree( expand(F2), {op(var_lst)});
deg_F3 := degree( expand(F3), {op(var_lst)});

deg_max := deg_F3;
if ( deg_F2 > deg_max-1 ) then
 deg_max := deg_F2 + 1;
end if;
if ( deg_F1 > deg_max-2 ) then
 deg_max := deg_F1 + 2;
end if;

return deg_max;

end proc:

# ---------------------------------------------------------------------------------------------
#
# Input: poly is a polynomial in the variables given by var_lst
#
# If poly is the square of another polynomial, this polynomial is returned.
# Otherwise, zero is returned.
#
CheckSquare := proc(poly,var_lst)

local res;

res := sqrt(poly,symbolic);

if ( type(res,polynom(anything,var_lst)) = false ) then
 return 0;
end if;

return res;

end proc:


# ---------------------------------------------------------------------------------------------
#
# Input: var_lst ist a list of variables, e.g. var_lst = [x1,x2,x3]
#        point_subs_lst ist a substitution list, e.g. point_subs_lst = [ x1=3, x2=7, x3=x3 ]
#
# Options: 
#
#  GeneralC          = true           solutions may depend on free parameters C_1, C_2, ...
#                    = false          the value zero is substituted for all occurring free parameters
#
# Ouput: A list of specific values for the variables, obtained by applying point_subs_lst to var_lst,
#        e.g. ConvertSubsLstToPoint([x1,x2,x3], [ x1=3, x2=7, x3=x3 ], GeneralC=true)  gives [3,7,C_3]
#             ConvertSubsLstToPoint([x1,x2,x3], [ x1=3, x2=7, x3=x3 ], GeneralC=false) gives [3,7,0]
#
ConvertSubsLstToPoint := proc(var_lst,point_subs_lst, {GeneralC::boolean := true})

local n_var,i1,point_lst,para_lst,flag_division_by_zero;

n_var := nops(var_lst);

point_lst := [];

for i1 from 1 to n_var do
 point_lst := [ op(point_lst), normal(subs(point_subs_lst,op(i1,var_lst))) ];
end do;

para_lst := [];
for i1 from 1 to n_var do
 if ( GeneralC ) then
  para_lst := [ op(para_lst), op(i1,var_lst) = SymbolFactory('C',i1) ];
 else
  para_lst := [ op(para_lst), op(i1,var_lst) = 0 ];
 end if;
end do;

# check for division by zero
flag_division_by_zero := false;
for i1 from 1 to n_var do
 if ( subs(para_lst,denom(op(i1,point_lst))) = 0 ) then
  flag_division_by_zero := true;
 end if;
end do;

# repeat with default value 1
if ( not GeneralC ) then
 para_lst := [];
 for i1 from 1 to n_var do
  para_lst := [ op(para_lst), op(i1,var_lst) = 1 ];
 end do;
end if;

# check again for division by zero
flag_division_by_zero := false;
for i1 from 1 to n_var do
 if ( subs(para_lst,denom(op(i1,point_lst))) = 0 ) then
  flag_division_by_zero := true;
 end if;
end do;

# give up
if ( flag_division_by_zero ) then
 return [];
end if;

return subs(para_lst,point_lst);

end proc:

# ---------------------------------------------------------------------------------------------
#
# Input: A point in projective space projective_lst = [ x1, x2, ..., xn ]
#
# Output: A point in the affine chart xn=1, i.e. [ x1/xn, x2/xn, ..., x_{n-1}/xn ]
#
ConvertProjectiveCoordinatesToAffineCoordinates := proc(projective_lst)

local n_var,i1,xn,affine_lst:

n_var := nops(projective_lst);

xn := op(n_var,projective_lst);

if ( xn=0 ) then
 error "cannot convert to standard affine coordinates";
end if;

affine_lst := [];

for i1 from 1 to (n_var-1) do
 affine_lst := [ op(affine_lst), op(i1,projective_lst)/xn ];
end do:

return affine_lst;

end proc:

# ----------------------------------------------------------------------------
#
# Input: poly is a polynomial in the first (n-1) variables of the list var_lst = [x1,x2,...,xn].
#
# Output: A homogeneous polynomial in the n variables [x1,x2,...,xn].
#
ConvertAffinePolynomialToProjectivePolynomial := proc(poly,var_lst)

local n_var,var_set,deg,my_poly_in,my_poly_out,i1,temp;

n_var := nops(var_lst);
var_set := {op(var_lst)};

my_poly_in := expand(poly);
my_poly_out := 0;

deg := degree( my_poly_in, var_set);

if ( type(my_poly_in,`+`) = false ) then
 return my_poly_in;
end if;

for i1 from 1 to nops(my_poly_in) do
 temp := op(i1,my_poly_in);
 my_poly_out := my_poly_out + temp * (op(n_var,var_lst))^(deg-degree(temp,var_set));
end do;

return my_poly_out;

end proc:

# ----------------------------------------------------------------------------
#
# Input: poly is a polynomial in the first (n-1) variables of the list var_lst = [x1,x2,...,xn].
#
# Output: A homogeneous polynomial of degree deg in the n variables [x1,x2,...,xn].
#
ConvertAffinePolynomialToProjectivePolynomialOfDegree := proc(poly,var_lst,deg)

local n_var,var_set,my_poly_in,my_poly_out,i1,temp;

n_var := nops(var_lst);
var_set := {op(var_lst)};

my_poly_in := expand(poly);
my_poly_out := 0;

if ( type(my_poly_in,`+`) = false ) then
 return my_poly_in * (op(n_var,var_lst))^(deg-degree(my_poly_in,var_set));
end if;

for i1 from 1 to nops(my_poly_in) do
 temp := op(i1,my_poly_in);
 my_poly_out := my_poly_out + temp * (op(n_var,var_lst))^(deg-degree(temp,var_set));
end do;

return my_poly_out;

end proc:

# ----------------------------------------------------------------------------
#
# This routine returns a symbol name.
# For example, SymbolFactory(c,[1,2,3])
# returns c_1_2_3.
#
# This works also with an atom:
# For example SymbolFactory(c,1)
# returns c_1.

SymbolFactory := proc(base_name,lst_indices)

local i1, temp;

temp := base_name;
for i1 from 1 to nops(lst_indices) do
 temp := cat(temp,"_",op(i1,lst_indices));
end do;

return convert(temp,symbol);
end proc:

# ----------------------------------------------------------------------------
#
# current_index = [ i1, i2, ..., ik ] is a multi-index.
# upper_limit = [ i1_max, i2_max, ..., ik_max ] defines the upper limit.
#
# The entries of current_index satisfy 0 <= ij <= ij_max
#
# The routine increases the counter by one (starting with i1) and returns 
# the carry_flag and the new value of the counter.
#
NextMultiIndex := proc( current_index, upper_limit )

local k,new_index,carry_flag,i1;

k := nops(current_index);

carry_flag := true;

new_index := [];

for i1 from 1 to k do
 if ( carry_flag = false ) then
   new_index := [ op(new_index), op(i1,current_index) ];
 else
  if ( op(i1,current_index) = op(i1,upper_limit) ) then
   new_index := [ op(new_index), 0 ]; # carry_flag remains to be set
  else
   new_index := [ op(new_index), op(i1,current_index)+1 ];
   carry_flag := false;
  end if;
 end if;
end do;

return [carry_flag,new_index];

end proc:

