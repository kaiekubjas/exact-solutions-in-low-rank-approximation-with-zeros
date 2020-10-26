# First of all we need a start system, namely one solution of our problem.
# To this aim, we pick a random $n\times n$ matrix $X\in\mathcal{L}_{n-1}^S$
# and we choose a random data point $U$ in the normal space $N_X\mathcal{L}_{n-1}^S$.
# In particular, $X$ is a critical point of $\delta_U$:

using LinearAlgebra;
n = 4; a = 10*randn(Float64,n-1); b = 10*randn(Float64,n-2);
push!(b, -sum([a[i]*b[i] for i=1:n-2])/a[n-1]);
A = [pushfirst!(10*randn(Float64,n-1), a[j]) for j=1:n-1];
B = [pushfirst!(10*randn(Float64,n-1), b[j]) for j=1:n-1];
X = sum([A[i]*(B[i])' for i=1:n-1]);
L, S, R = svd(X, full=true);
Lk = L[:,n]/(L[n,n]); Rk = R[:,n]/(R[n,n]);
N = vec(Lk*(Rk'));
dL = pushfirst!([0.0 for i=1:n^2-1],1);
x0 = vec(X); y0 = Lk[1:n-1]; z0 = Rk[1:n-1]; w0 = 10*randn(Float64,2);
u0 = vec(X)+hcat(N,dL)*w0;

# Hence the start system is defined by the start solution (x0,y0,z0,w0) and by the start data point u0.
# Secondly, we define symbolically the polynomial system whose solutions
# are the critical points of $\delta_U$ on $\mathcal{L}_{n-1}^S$:

using HomotopyContinuation;
@polyvar x[1:n,1:n] u[1:n,1:n] y[1:n-1] z[1:n-1] w[1:2];
D = vec(x-u); Y = [y;1]; Z = [z;1];
N = vec(Y*(Z'));
l = x[1,1];
dL = differentiate([l], x)';
W = [w; 1];
M = hcat(N,dL,D)*W;
eq = [vec(Y'*x); vec(x*Z); [l]; vec(M)];

# Now we are ready to solve the polynomial system \verb|eq| numerically.
# First, we pick a random data matrix \verb|eU| and we find a new critical point
# using the start solution computed before:

eU = 10*randn(ComplexF64, n^2);
res_nice = solve(eq, [x0; y0; z0; w0]; parameters = vec(u),
start_parameters = u0, target_parameters = eU);

# Finally we recover all critical points from the first critical point obtained:

sol = monodromy_solve(eq, solutions(res_nice), eU, parameters = vec(u),
target_solutions_count = 5*(n-1)-2);
listsol = solutions(sol);
