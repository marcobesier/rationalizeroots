(* ::Package:: *)

BeginPackage["RationalizeRoots`"];
RationalizeRoot::usage="\"RationalizeRoot\" takes a the n-th root of a rational expression as input, possibly muliplied
by another rational expression and returns rational parametrizations if found.";

FindPoints::usage="Find all singular points on the hypersurface \"poly\[Equal]0\" with multiplicity \"deg-1\". 
The points may depend on free parameters C[1],...,C[n].";

ParametrizePolynomial::usage="\"ParametrizePolynomial\" takes a polynomial in n variables and returns a rational parametrization in terms of n-1 variables.
It automatically searches for a singular point of mulitplicity d-1 (where d is the degree of the polynomial) in the procective closure of the space.
If no such point is found it searches for all possible decomposisions in F1, F2, F3, such that u=Sqrt[F2^2-4*F1*F3] and uses the theorem in the appendix 
of arxiv:809.10983 to rationalize F1+F2+F3. 
Setting Variables\[Rule]{...} the list of actual variables can be specified. All other variables will the be treated as parameters that do not need to be rationized.
As an output a list of solutions is provided. By default only one solution is provided. Setting the options AllCharts, AllPoints and/or AllFDecomposition to 
	true, several solutions will be provides. The solution may contains free parameters (C[1], C[2],...). By default they are set to simple values, 
but setting GeneralC to True the solution will leave them open.
Setting GeneralT to true, the solution will depend on n instead of n-1 variables where one variable must be set to some value (usually one)";

ProjectiveClosure::usage="Turns poly into a polynomial of homogenous degree with the variables of poly and pvar";

RationalParametrization::usage;
GeneralC::usage;
GeneralT::usage;
AllPoints::usage;
AllCharts::usage;
MultipleSolutions::usage;
RationalizeFDecomposition::usage;
FindFDecomposition::usage;
TryFDecomposition::usage;
AllFDecompositions::usage;
ForceFDecomposition::usage;
FPolynomials::usage;
OutputVariables::usage;
ParamT::usage;
RootOutput::usage;


Begin["`Private`"]


DebugMode=False;
DebugPrint[expr_]:=If[DebugMode,Print[expr],Null];


(*Turns the affine polynomial "poly" into a homogenous polynomial with the additional variable "pvar"*)
Options[ProjectiveClosure]={Variables->Null};
ProjectiveClosure::argx="ProjectiveClosure called with 1 argument. 2 or 3 arguments expected.";
ProjectiveClosure[a_]:=Message[ProjectiveClosure::argx];
ProjectiveClosure[poly_, pvar_,OptionsPattern[]]:=Module[{vs,deg,w},
	vs=If[OptionValue[Variables]===Null,Variables[poly],OptionValue[Variables]];
	deg=Exponent[poly/.MapThread[Rule,{vs,w vs}],w];
	pvar^deg (poly/.MapThread[Rule,{vs,vs/pvar}])//Expand
]


(*Find all singular points on the hypersurface "poly\[Equal]0" with multiplicity "deg-1". If
poly is a projective polynomial the trivial solution x1=0,...,xn=0 must be excluded. The points may
depend on free parameters C[1],...,C[n].*)
Options[FindPoints]={Variables->Null};
FindPoints[poly_,OptionsPattern[]]:=Module[{deg,deriv,w,vs,sol,points,free},
	vs=If[OptionValue[Variables]===Null,Variables[poly],OptionValue[Variables]];
	deg=Exponent[poly/.MapThread[Rule,{vs,w vs}],w];
	deriv=Prepend[Table[D[poly,{vs,i}],{i,deg-1}],poly];
	sol=Quiet[Solve[deriv[[1;;-2]]==0,vs]];
	(*Delete solutions, where next derivative vanishes*)
	sol=DeleteCases[sol,a_/;Factor[deriv[[-1]]/.a]==(0*deriv[[-1]])];
	If[Length[sol]==0,Return[{}]];
	(*sol=DeleteCases[sol,a_/;a[[All,2]]==ConstantArray[0,Length[vs]]];*)
	points=vs/.sol;
	Table[free=Intersection[vs,Variables[points[[i]]]];MapThread[Rule,{vs,points[[i]]/.MapThread[Rule,{free,C/@Range[Length[free]]}]}],{i,Length[points]}]
]


FindPoints[poly_,z_,OptionsPattern[]]:=Module[{vs},
	vs=If[OptionValue[Variables]===Null,Union[Variables[poly],{z}],Union[OptionValue[Variables],{z}]];
	FindPoints[ProjectiveClosure[poly,z],Variables->vs]
]


RationalParametrization::OutputVariables="OutputVariables must either be a list of the same length as the number of input variables or a symbol. Now using t or ParamT as output variable."
RationalParametrization::VariableConflict="OutputVariables must not have same values as the input variables, which are specified by 'point'"
Options[RationalParametrization]={OutputVariables->Null};
RationalParametrization[ppoly_,point_,cvar_,OptionsPattern[]]:=Module[{p,wpoly,affineshifted,vs,w,we,cl,fnum,fden,out,tfree,outv},
	(*Assign output variables. OutputVariables can be a list of variables or a symbol (say "s"). In the latter case the output variables are {s[1], s[2], ...}.
	If the specifaction is wrong, a warning is displayed and the default choice is used. The default choice is {t[1],t[2],...} if 
	none of them are already defined by the use and otherwise {ParamT[1],ParamT[2],...}*)
	(*The output variables must not be the same as the input variables for RationalParametrization[]. Note that for ParametrizePolynomial[] and RationalizeRoot[] they may be the same*)	
	out=Null;
	tfree=And@@Table[ToString[Global`t[i]]===("t["<>ToString[i]<>"]"),{i,Length[point]}];
	If[OptionValue[OutputVariables]===Null,out=ParamT/@Range[Length[point]-1]/.If[tfree,ParamT->Global`t,{}]
		,
		If[Head[OptionValue[OutputVariables]]===Symbol,
			out=OptionValue[OutputVariables]/@Range[Length[point]-1],
			If[Head[OptionValue[OutputVariables]]===List&&Length[OptionValue[OutputVariables]]==Length[point]-1,out=OptionValue[OutputVariables]]
		]
	];
	If[out===Null,Message[RationalParametrization::OutputVariables]; out=ParamT/@Range[Length[point]-1]/.If[tfree,ParamT->Global`t,{}]];	

	If[Length[Intersection[out,point[[All,1]]]]>0,Message[RationalParametrization::VariableConflict];Return[False]];
	If[(cvar/.point)===0,Print[cvar," must not be zero"];Abort[],p=MapThread[Rule,{point[[All,1]],point[[All,2]]/(cvar/.point)}]];
	If[wpoly=(ppoly/.MapThread[Rule,{point[[All,1]],w point[[All,1]]}]);(Factor[wpoly]/.w^(Exponent[wpoly,w]):>0)=!=0,
		Print["Internal error. Polynomial must be of homogenous degree in all variables"];Abort[];];

	affineshifted=Expand[ppoly/.cvar->1/.MapThread[Rule,{p[[All,1]],p[[All,1]]+p[[All,2]]}]];
	vs=DeleteCases[point[[All,1]],cvar];
	we=affineshifted/.MapThread[Rule,{vs,w vs}];
	cl=CoefficientList[we,w]//Factor;
	If[cl[[1;;-3]]=!=ConstantArray[0,Length[cl]-2],Print["Degree of singularity too small"];Abort[]];
	{fnum,fden}=cl[[{-2,-1}]]//Expand;
	(*t=Global`t;*)
	(vs/cvar)==Factor[((out)*(-fnum/fden/.MapThread[Rule,{vs,out}])+(vs/.p))]
]


MyCollect[expr_, vars_, func_] := Module[{Ruler},
	FromCoefficientRules[CoefficientRules[expr, vars]/.Rule->Ruler/.Ruler[a_,b_]:>Rule[a,func[b]], vars]
]


SumToList[term_] := If[Head[term] === Plus, List @@ term, {term}] 
ProductToList[term_] := If[Head[term] === Times, List @@ term, {term}]


(*
**This function takes a polynomial in n variables and returns a rational parametrization in terms of n-1 variables
**It automatically searches for a singular point of mulitplicity d-2 (where d is the degree of the polynomial) in the procective closure of the space
**If no such point is found it searches for all possible decomposisions in F1, F2, F3, such that u=Sqrt[F2^2-4*F1*F3] and uses the theorem in the appendix to
	rationalize F1+F2+F3. 
**Setting Variables\[Rule]{...} the list of actual variables can be specified. All other variables will the be treated as parameters that do not need to be rationized
**As an output a list of solutions is provided. By default only one solution is provided. Setting the options AllCharts, AllPoints and/or AllFDecomposition to 
	true, several solutions will be provided
**Sometimes the solution contains free parameters (C[1], C[2],...). By default they are set to simple values, but setting GeneralC to True the solution will leave them open
**Setting GeneralT to true, the solution will depend on n instead of n-1 variables where one variable must be set to some value (usually one)
*)
Options[ParametrizePolynomial]={Variables->Null,GeneralC->False,GeneralT->False,AllCharts->False,AllPoints->False,TryFDecomposition->True,
	AllFDecompositions->False, OutputVariables->Null, MultipleSolutions->False,ForceFDecomposition->False,FPolynomials->Null};
ParametrizePolynomial::variables="Polynomial must depend on all specified variables."
ParametrizePolynomial::polynomial="Input is not a polynomial in all (or specified) variables."
ParametrizePolynomial::Fpolynomials="FPolynomials must be a list of 3 polynomials {F1,F2,F3}, "<>
	"such that poly=u^2-(F2^2-4F1*F3) or the square root has the form Sqrt[F2^2-4F1*F3], and {F1,F2,F3} must not depend on 'u'."
ParametrizePolynomial::OutputVariables="OutputVariables must either be a list of the same length as the number of input variables or a symbol. Now using t or ParamT as output variable."
ParametrizePolynomial[poly_,OptionsPattern[]]:=Module[{hpoly,vs,points,z,pp,vv,sol,fd,ret,polyexpanded,rvs,pointsC, soll, out, tfree, outT, outt,
		variables, interntov, vtointern, ipoly, vvv, allcharts, allpoints, fpols, uv},
	(*Switch to intern variables*)
	(*v=Global`v;*)
	(*OptionValue[AllCharts]=True;*)
	allcharts=OptionValue[AllCharts]||OptionValue[MultipleSolutions];
	allpoints=OptionValue[AllPoints]||OptionValue[MultipleSolutions];
	variables=Variables[poly];

	vtointern=Table[variables[[i]]->vvv[i],{i,Length[variables]}];
	interntov=Table[vvv[i]->variables[[i]],{i,Length[variables]}];
	ipoly=poly/.vtointern;

	(*If no variables are specified by default all variables in the polynomials are taken as variables*)	
	vs=If[OptionValue[Variables]===Null,Variables[ipoly],OptionValue[Variables]/.vtointern];

	(*Assign output variables. OutputVariables can be a list or a symbol (say "s"). In the latter case the output variables are {s[1], s[2], ...}.
	If the specifaction is wrong, a warning is displayed and the default choice is used. The default choice is {t[1],t[2],...} if 
	none of them are defined differently and otherwise {ParamT[1],ParamT[2],...}. The length of "out" is one more element than input variables in case the GeneralT option is used.
	If for n input variables only n-1 output variables are specified, the n-th variable is by default ParamT[n]*)
	(*Internally this function uses the private variables 'outT' as the output variables to avoid conflicts with the input variables. Only in the last step they are replaced by 
	the output variables 'out'. This means the output variables may also be the same as the input variables*)
	outT=Array[outt,Length[vs]];
	out=Null;
	tfree=And@@Table[ToString[Global`t[i]]===("t["<>ToString[i]<>"]"),{i,Length[vs]}];
	If[OptionValue[OutputVariables]===Null,out=ParamT/@Append[Range[Length[vs]-1],0]/.If[tfree,ParamT->Global`t,{}]
		,
		If[Head[OptionValue[OutputVariables]]===Symbol,
			out=OptionValue[OutputVariables]/@Range[Length[vs]],
			If[Head[OptionValue[OutputVariables]]===List&&MemberQ[{Length[vs],Length[vs]-1},Length[OptionValue[OutputVariables]]],out=OptionValue[OutputVariables]]
		]
	];
	If[out===Null,Message[ParametrizePolynomial::OutputVariables]; out=ParamT/@Append[Range[Length[vs]-1],0]/.If[tfree,ParamT->Global`t,{}]];
	If[Length[out]==Length[vs]-1,out=Append[out,ParamT[0]/.If[tfree,ParamT->Global`t,{}]]];
	(*Expand polynomial in all variables and factorize the coefficients*)
	polyexpanded=MyCollect[ipoly,vs,Factor];
	(*Check if poly is really a polynomial in all variables*)
	If[!PolynomialQ[polyexpanded,vs],Message[ParametrizePolynomial::polynomial,vs];Return[{}]];
	(*Check if poly depends on all specified variables*)
	If[!SubsetQ[Variables[ipoly],vs],Message[ParametrizePolynomial::variables,vs];Return[{}]];
	hpoly=ProjectiveClosure[ipoly,z,Variables->vs];
	(*If polynomial is accidentally homogenous, z will still be included as a variable, but we have to make sure not 
	to choose a point, where all variables (except z, which can be arbitraty) are zero*)
	vs=Union[vs,{z}];

	(*points=DeleteCases[FindPoints[hpoly,Variables->vs],a_/;Union[a[[All,2]]]==={0}];*)
	(*Find singular points of multiplicity d-2*)
	pointsC=If[OptionValue[ForceFDecomposition]||OptionValue[FPolynomials]=!=Null,
		{},
		FindPoints[hpoly,Variables->vs]
	];
	If[Length[pointsC]==0,
		ret={},
		(*Set free parameters C[1],C[2],... if GeneralC option is not used*)
		If[OptionValue[GeneralC],
			points=pointsC;,
			points=Join@@Table[Quiet[Factor[pointsC[[i]]/.C[j]->1/.C[_]:>0],{Power::infy,Infinity::indet}],{i,Length[pointsC]},{j,Length[Cases[Variables[pointsC[[i,All,2]]],_C]]}];
		];
		(*Delete cases with infinity or inderminate values*)
		points=DeleteCases[points,a_/;!FreeQ[a[[All,2]],Infinity|ComplexInfinity|Indeterminate]];
		(*There might be cases where no solution is found unless GeneralC is set to True*)
		(*Delete duplicate solutions and solutions the non projective solutions at all variables zero.*)
		rvs=If[FreeQ[hpoly,z],DeleteCases[vs,z],vs];
		points=DeleteDuplicates[DeleteCases[points,a_/;Union[rvs/.a]==={0}]];
		If[Length[points]==0,ret={},		
			ret=Join@@Table[
				(*pp=If[OptionValue[GeneralC],point,point/.C[1]->1/.C[_]:>0];*)
				pp=point;
				(*Each non-zero variable can be used as a chart. If option AllCharts is True, all charts are used. 
				Otherwise only the last non-zero variable of the point*)
				(*vv: non-zero variables*)
				vv=pp[[DeleteCases[Position[pp[[All,2]],a_/;a=!=0,{1}][[All,1]],0],1]];
				(*In cases where g_{d-1}=0 by accident, there is no solution to the equations given by RationalParametrization*)
				sol=Table[soll=Solve[RationalParametrization[hpoly,pp,v,OutputVariables->outT]/.z->1,DeleteCases[vs,z]];
					If[Length[soll]!=1,Nothing,soll[[1]]],{v,If[allcharts,vv,{vv[[-1]]}]}]/.
						If[OptionValue[GeneralT],{},outT[[Length[vs]-1]]:>1]
				(*If[OptionValue[GeneralSolution],sol,sol/.t[Length[Variables]]\[RuleDelayed]1]*)
					,{point,If[allpoints,points,{points[[1]]}]}	
			];
		];
	];
	(*If no solution was found to this stage we try to apply the Theorem in appendix C of arxiv:1809.10983*)
	If[Length[ret]==0 && OptionValue[TryFDecomposition],
		DebugPrint["Try FDecomposition"];
		fd=If[OptionValue[FPolynomials]===Null,
			FindFDecomposition[ipoly,Variables->(OptionValue[Variables]/.vtointern)],
			fpols=OptionValue[FPolynomials]/.vtointern;
			If[Head[fpols]=!=List||Length[fpols]!=3||!(And@@(PolynomialQ[#,vs]&/@fpols)),Message[ParametrizePolynomial::Fpolynomials];Return[{}]];
			uv=Factor[ipoly+fpols[[2]]^2-4fpols[[1]]fpols[[3]]];
			If[MatchQ[uv,a_. b_^2/;(MemberQ[vs,b]&&FreeQ[a,Alternatives@@vs]&&FreeQ[fpols,b])], 
				{{PowerExpand[Sqrt[uv]],fpols[[2]],fpols[[1]],fpols[[3]]}}
				,
				Message[ParametrizePolynomial::Fpolynomials];
				{}
			]			
		];
		ret={};
		If[Length[fd]>0,
			Do[
				ret=Join[ret,RationalizeFDecomposition[f,Variables->(OptionValue[Variables]/.vtointern),
					GeneralC->OptionValue[GeneralC],GeneralT->OptionValue[GeneralT],
					AllCharts->allcharts,AllPoints->allpoints,MultipleSolutions->OptionValue[MultipleSolutions],OutputVariables->outT]];
				If[Length[ret]>0&&(!OptionValue[AllFDecompositions]),Break[]];
				,{f,fd}
			];
		];
	];

	(*Return 'ret' and replace internal output variables 'outT' by 'out'*)
	Table[MapThread[Rule,{ret[[i,All,1]],ret[[i,All,2]]/.MapThread[Rule,{outT,out}]}],{i,Length[ret]}]/.interntov
	
]


Options[RationalizeRoot]={Variables->Null,GeneralC->False,GeneralT->False,AllCharts->False,AllPoints->False,TryFDecomposition->True,
	AllFDecompositions->False, OutputVariables->Null,RootOutput->False,MultipleSolutions->False,FPolynomials->Null,ForceFDecomposition->False};
RationalizeRoot[pref_. Power[root_,exp_Rational],OptionsPattern[]]:=Module[{r,f,u,vs,res},
	r=Factor[root pref^(1/exp)];
	f=u^(1/exp)Denominator[r]-Numerator[r];
	vs=If[OptionValue[Variables]===Null,DeleteCases[Variables[f],u],OptionValue[Variables]];
	res=ParametrizePolynomial[f,Variables->Append[vs,u],GeneralC->OptionValue[GeneralC],GeneralT->OptionValue[GeneralT],
		AllCharts->OptionValue[AllCharts],AllPoints->OptionValue[AllPoints],TryFDecomposition->OptionValue[TryFDecomposition],
		AllFDecompositions->OptionValue[AllFDecompositions],OutputVariables->OptionValue[OutputVariables],MultipleSolutions->OptionValue[MultipleSolutions],
		FPolynomials->OptionValue[FPolynomials],ForceFDecomposition->OptionValue[ForceFDecomposition]
	];
	If[OptionValue[RootOutput],res/.u->pref Power[root,exp],DeleteCases[#,a_/;a[[1]]===u]&/@res]
]


(*This function expands "poly" and search for all possible decompositions of the form u^2-(F2^2-4F1 F3).
Returns a list of lists {{u,F2,F1,F3},...}, where the dots indicate further solutions.*)
Options[FindFDecomposition]={Variables->Null};
FindFDecomposition[poly_,OptionsPattern[]]:=Module[{u,vs,polyexp,termlist,ucands,polynew,listpairs,Fdecs},

	vs=If[OptionValue[Variables]===Null,Variables[poly],OptionValue[Variables]];
	(*Expands the polynomial as a sum of monomials of the variables "vs". The coefficients of the monomials
	may depend on further constant parameters and are factorized*)
	polyexp=MyCollect[poly,vs,Factor];
	termlist=If[Head[polyexp]===Plus,List@@polyexp,{polyexp}];
	DebugPrint["Length termlist ",Length[termlist]];
	(*If the polynomial has too many terms this function is very inefficient and is aborted*)
	If[Length[termlist]>20,DebugPrint["Polynomial too large. FDecomposition will not be searched"];Return[]];
	(*Search squared terms as candidates for "u^2"*)
	ucands=Cases[termlist,a_. b_^2/;(MemberQ[vs,b]&&FreeQ[a,Alternatives@@vs]&&Length[Position[termlist,b]]==1)];
	If[Length[ucands]==0,DebugPrint["No FDecomposition found"];Return[]];
	(*For every candidate "u" search for decompositions of the rest of the polynomial in terms of F2^2-4F1 F3 by considering
	all possible decomposition into two expressions*)
	DeleteDuplicates@(Join@@Table[
		polynew=-polyexp+uc;
		termlist=If[Head[polynew]===Plus,List@@polynew,{polynew}];
		(*List of pairs (a,b), where a is perfect square, and b is a list of factors of the rest of the polynomial*)
		listpairs=DeleteCases[{PowerExpand[Sqrt[Factor[Total[#]]]],(#[[1]]^#[[2]])&/@FactorList[Total[Complement[termlist,#]]]}&
			/@DeleteCases[Subsets[termlist],{}|termlist],a_/;!FreeQ[a,Power[b_,n_]/;!FreeQ[b,Alternatives@@vs]&&!IntegerQ[n]]];
		(*expand factors like {y^3} out to {y,y,y}*)
		listpairs=Replace[listpairs, {a_, {b1___, b2_^n_, b3___}} /; (IntegerQ[n] && n > 1) :> {a, {b1, Sequence @@ ConstantArray[b2, n],  b3}}, {1}];
		(*DebugPrint[{"Length listpairs ",Length[listpairs]}];*)
		Fdecs=Join@@Table[{2PowerExpand[Sqrt[Factor[uc]]],2listpairs[[i,1]],-Times@@#,Times@@ListComplement[listpairs[[i,2]],#]}&
			/@DeleteCases[Subsequences[listpairs[[i,2]]],{}|listpairs[[i,2]]],{i,Length[listpairs]}]
		,{uc,ucands}
	])
]

ListComplement[a_List, b_List] := Module[{aa, fp},
  aa = a;
  Do[
   fp = FirstPosition[aa, u, Null, {1}][[1]];
   aa = Delete[aa, fp]
   , {u, b}];
  aa
]

(*Takes as input a list of list of polynomials {{u,F2,F1,F3},...}, which is the output format of "FindFDecomposition" function.
Returns *)
Options[RationalizeFDecomposition]={Variables->Null,GeneralC->False,GeneralT->False,AllCharts->False,AllPoints->False,
	OutputVariables->Null,MultipleSolutions->False};
RationalizeFDecomposition[fdec_,OptionsPattern[]]:=Module[{vs,u,F1,F2,F3,uv,degs,rh,vdz,z},
	vs=If[OptionValue[Variables]===Null,Variables[fdec],OptionValue[Variables]];
	{u,F2,F1,F3}=Expand[{fdec[[1]],ProjectiveClosure[fdec[[2]],z,Variables->vs],
		ProjectiveClosure[fdec[[3]],z,Variables->vs],ProjectiveClosure[fdec[[4]],z,Variables->vs]}];
	uv=Intersection[Variables[u],vs];
	degs=Exponent[((SumToList/@{F1,F2,F3})[[All,1]])/.MapThread[Rule,{vs,ConstantArray[z,Length[vs]]}],z]+{2,1,0};
	degs=-degs+ConstantArray[Max[degs],Length[degs]];
	{F1,F2,F3}={F1,F2,F3}*(z^degs);
	(*If F1, F2, F3 accidentally do not depend on z, return no solution.*)
	If[FreeQ[{F1,F2,F3},z],DebugPrint["No z dependence in F1, F2, F3"]; Return[{}]];
	rh=ParametrizePolynomial[F1+F2+F3,Variables->Append[DeleteCases[vs,Alternatives@@Variables[fdec[[1]]]],z],
		GeneralC->OptionValue[GeneralC],GeneralT->OptionValue[GeneralT],
		AllCharts->OptionValue[AllCharts],AllPoints->OptionValue[AllPoints],TryFDecomposition->True,
		MultipleSolutions->OptionValue[MultipleSolutions],OutputVariables->OptionValue[OutputVariables]
	];
	Table[
		vdz=MapThread[Rule,{rs[[All,1]],rs[[All,1]]/z}];
		Join[{Factor@Solve[fdec[[1]]-(2(F3/.vdz) z+(F2/.vdz)/.rs)==0,uv][[1,1]]},
			MapThread[Rule,{DeleteCases[rs[[All,1]],z],DeleteCases[rs[[All,1]],z]/z/.rs}]]
		,{rs,rh}
	]
]


End[];
EndPackage[];
Print["RationalizeRoots: Marco Besier, Pascal Wasser, and Stefan Weinzierl (2019)"];
