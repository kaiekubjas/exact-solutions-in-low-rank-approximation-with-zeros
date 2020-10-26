restart;

# Load the linear algebra package the FGb library for Groebner basis computations

with(LinearAlgebra):
with(FGb):

# A procedure to compute the dimension and the degree of an algebraic variety
# defined by the polynomials in the list L

fgb_degree := proc(L)
local p, vars, gb, h, DD, dim;
    p:=65521:

    vars := [ op(indets(L)) ]:
    gb:=fgb_gbasis_lm([op(L)], p, [], vars, {"index"=10^7, "verb"=0});
    # gb:=fgb_gbasis_lm([op(L)], p, [], drl( op(vars))):

    gb := gb[2];
    # vars := [ op(indets( gb))];
    h := fgb_hilbert( gb, p, [], vars, T, {"index"=10^7, "verb"=0}); # {"index"=10^7, "verb"=1});
    DD := subs(T=1, h[1]):

    dim := h[2];
    # printf( "h: %a \n", h);
    # print(DD, dim);
    # printf( "D = %2d \t dim: %2d \n", DD, dim);

    return DD;
end proc:




# m x n matrix of rank r
# L is the linear constraints

GetDegree := proc(m, n, r, L)

local X, vars, Y1, Z1, Imr, Inr, Y, Z, s, N, NN, U, tmp1, f_obj, M1, M2, M, W, eq1, eq2, eq3, eqns;
  X := Matrix(m, n, symbol=x);
  vars := [ op(indets(X))];

  Y1 := Matrix(r, m-r, symbol=y);
  Z1 := Matrix(r, n-r, symbol=z);
  Imr := IdentityMatrix(m-r):
  Inr := IdentityMatrix(n-r):
  Y := Matrix([[Imr], [Y1]]); # print("Y", Y);
  Z := Matrix([[Inr], [Z1]]); # print("Z", Z);

  # Linear constraints
  # L := [x[1,1]]:
  s := nops(L);

  N := [ seq( seq( Column(Y, i) . Row(Transpose(Z), j), i=1..m-r), j=1..n-r) ];
  NN := seq( convert(Transpose(N[i]), list), i=1..nops(N));
 # print("NN", N, NN);

  U := RandomMatrix(m, n);
  tmp1 := ListTools[Flatten]( convert(X - U, list));
  f_obj := add( tmp1[i]^2, i=1..nops(tmp1));

  M1 := Matrix([NN]);
  M2 := VectorCalculus[Jacobian] ([op(L), f_obj], vars);
# print(Dimensions(M1), Dimensions(M2));
  M := Matrix([[M1], [M2]]);

  W := < seq(w[i], i=1..(m-r)*(n-r)+s), 1 >;

  eq1 := convert( Transpose(W) . M, list);
  eq2 := convert( Transpose(Y) . X, list);
  eq3 := convert( X . Z, list);;

  eqns := [ op(eq1), op(eq2), op(L), op(eq3) ];

  fgb_degree( convert(eqns, list));

end proc:


Diagonal_constraints := [x[1,1], x[1,2]];

N:= 5:

for m from 3 to N do
  for n from 3 to N do
    DD[m,n] := GetDegree(m, n, 2, Diagonal_constraints);
od; od;



# print the results a table format
printf("    ");
for n from 3 to N do
  printf("%4d ", n);
od;
printf("\n");

for m from 3 to N do

  printf("%4d", m);
  for n from 3 to N do
    printf("%4d ", DD[m,n]);
  od;
  printf("\n");

od;
