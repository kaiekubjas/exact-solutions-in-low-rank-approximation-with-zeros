-- m x n case, U picked randomly
-- we consider the constrained low-rank approximation problem,
-- when the constraints are given by l random linear relations in the entries x_(i,j)
-- Goal: compute all linear (and affine) relations among critical points
-- Conjecture: if l is nonzero, the unique linear relations are the first l given.
-- In other words, the linear span of the critical points coincides with the subspace L
-- of codimension l defined by the linear relations.r
-- we assume m <= n for simplicity
-- l is the number of random linear relations considered

restart

linearRelations = (m,n,r,l) -> (
    KK := ZZ/101;
    R := KK[x_(1,1)..x_(m,n),
	    y_(1,1)..y_(r,m-r),z_(1,1)..z_(r,n-r),
	    w_1..w_((m-r)*(n-r)+l)];
    X := matrix for i in 1..m list for j in 1..n list x_(i,j);
    U := matrix for i in 1..m list for j in 1..n list random(KK);
    Y := id_(KK^(m-r))||(matrix for i in 1..r list for j in 1..m-r list y_(i,j));
    Z := id_(KK^(n-r))||(matrix for i in 1..r list for j in 1..n-r list z_(i,j));
    for i in 1..(m-r) do for j in 1..(n-r) do N_((m-r)*(j-1)+i) = matrix{flatten entries((matrix Y_(i-1))*(transpose matrix Z_(j-1)))};
    matN := N_1; for i in 2..(m-r)*(n-r) do matN = matN||N_i;
    L := flatten entries(random(KK^l,KK^(m*n))*transpose(matrix{{x_(1,1)..x_(m,n)}}));
    D := matrix{flatten entries(X-U)};
    W := matrix{toList(w_1..w_((m-r)*(n-r)+l))|{1}};
    if #L!=0 then (
	dL := diff(transpose matrix{L},matrix{{x_(1,1)..x_(m,n)}});
	M = W*(matN||dL||D);
	) else (
	M = W*(matN||D)
	);
    J := ideal((flatten entries((transpose Y)*X))|(flatten entries(X*Z))|L|(flatten entries M));
    el := eliminate(toList(w_1..w_((m-r)*(n-r)+l))|toList(y_(1,1)..y_(r,m-r))|toList(z_(1,1)..z_(r,n-r)),J);
    G := first entries gens el;
    result := {};
    for f in G do if (degree(f))#0==1 then result=append(result,f);
    return (result, L, U)
    )

-- examples
dim1 = 3
dim2 = 4
r = 2
l = 1
S = ZZ/101[x_(1,1)..x_(dim1,dim2)]

time (Res, LL, dataU) = linearRelations(dim1, dim2, r, l);
#Res
#LL
rank contract(matrix{{x_(1,1)..x_(dim1,dim2)}}, transpose matrix{Res|LL})
