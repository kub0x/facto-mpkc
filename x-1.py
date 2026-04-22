#!/usr/bin/env python3
from sage.all import *
from itertools import combinations

p = 11
n = 7
F = GF(p)

def random_invertible_matrix(F, r):
    while True:
        M = random_matrix(F, r, r)
        if M.is_invertible():
            return M

def multiplication_matrix(F, n):
    """
    Build M of size (2n-1) x n^2 such that for u,v in F^n,

        M * (u ⊗ v)

    gives the coefficient vector of u(t)*v(t).

    Tensor order:
        u ⊗ v = [u0*v0, u0*v1, ..., u0*v_{n-1},
                 u1*v0, u1*v1, ..., u1*v_{n-1},
                 ...
                 u_{n-1}*v_{n-1}]
    """
    rows = 2 * n - 1
    cols = n * n
    M = matrix(F, rows, cols)

    for i in range(n):
        for j in range(n):
            k = i + j
            col = i * n + j
            M[k, col] = F(1)

    return M

def tensor_vec_symbolic(a, b):
    """
    Return a ⊗ b as a column vector over the same parent ring.
    """
    R = a[0].parent()
    return vector(R, [ai * bj for ai in a for bj in b])

def poly_from_coeffs(coeffs, t):
    s = t.parent()(0)
    for i, c in enumerate(coeffs):
        s += c * (t ** i)
    return s

def find_subset_sum(terms, target):
    for r in range(1, len(terms) + 1):
        for subset in combinations(range(len(terms)), r):
            if sum(terms[i] for i in subset) == target:
                yield subset

def extract_coeffs_in_var_as_vector(P, var, n):
    """
    Given a multivariate polynomial P and a variable 'var',
    extract the coefficients as a vector of length n:
    [coeff of t^0, coeff of t^1, ..., coeff of t^{n-1}].

    If a degree is missing, the coefficient is 0.
    """
    R = P.parent()
    vars = R.gens()
    idx = vars.index(var)
    
    coeffs = {}

    for mon, c in P.dict().items():
        d = mon[idx]
        mon_wo_var = list(mon)
        mon_wo_var[idx] = 0
        coeff_poly = R({tuple(mon_wo_var): 1}) * c

        if d not in coeffs:
            coeffs[d] = coeff_poly
        else:
            coeffs[d] += coeff_poly
            
    coeff_list = []
    for i in range(n):
        coeff_list.append(coeffs.get(i, R(0)))

    return vector(R, coeff_list)

R = PolynomialRing(F, [var('x' + str(i)) for i in range(1, 2*n+1)] + ['t'])
vars = R.gens()
t = vars[-1] 

X = vector(R, vars[:2 * n])

S1 = random_invertible_matrix(F, n)

S =  block_matrix([[S1, zero_matrix(n,n)],[zero_matrix(n,n), S1]])
invS = S.inverse()

L1 = S[0:n, :]
L2 = S[n:2*n, :]

T = random_invertible_matrix(F, 2 * n - 1)

M = multiplication_matrix(F, n)

K = L1.tensor_product(L2)

Pmat = T * M * K

XX = tensor_vec_symbolic(X, X)
P_of_X = Pmat * XX

eqs = []
for i, eq in enumerate(P_of_X):
    eqs.append(eq)
print()

in_x = vector(F, [F.random_element() for _ in range(2 * n)])
print("Random input X0 =")
print(in_x)
print()

eval_x = {f'{vars[i]}' : in_x[i]  for i in range(2*n)}

output = vector(P_of_X.subs(**eval_x))


print("output = P(X0) =")
print(output)
print()

out2 = T.inverse()*output
outpoly = sum(out2[i] * t**i for i in range(2*n-1))
factorization = outpoly.factor()
terms = [f.degree() * e for (f, e) in factorization]

seen = set()

found = False

if outpoly.degree() != 2*n-2:
    print("re-run, cannot decompose into two n-1 polynomials")
    sys.exit(-1)

for solution_indices in find_subset_sum(terms, n-1):
    if solution_indices not in seen:
        seen.add(solution_indices)

        solution1 = 1
        solution2 = 1

        for i in range(len(solution_indices)):
            solution1 *= (factorization[solution_indices[i]][0])**(factorization[solution_indices[i]][1])

        for i in range(len(factorization)):
            if (i not in solution_indices):
                solution2 *= (factorization[i][0])**(factorization[i][1])

        assert(factorization.unit()*solution1*solution2 == outpoly)

        y1 = extract_coeffs_in_var_as_vector(solution1, t, n)
        y2 = extract_coeffs_in_var_as_vector(solution2, t, n)
        

        u = factorization.unit()

        uy1 = u*y1
        uy2 = u*y2

        sol1 = invS * vector(uy1.list() + y2.list())
        sol2 = invS * vector(y1.list() + uy2.list())

        sol3 = invS * vector(y2.list() + uy1.list())
        sol4 = invS * vector(uy2.list() + y1.list())

        candidates = [sol1, sol2, sol3, sol4]

        def eval_candidate(sol):
            eval_x = {str(vars[i]): sol[i] for i in range(2*n)}
            return vector(P_of_X.subs(**eval_x))

        matching = [(idx, sol, eval_candidate(sol)) for idx, sol in enumerate(candidates, 1)]
        matching = [(idx, sol, val) for idx, sol, val in matching if val == output]

        if len(matching) > 0:
            found = True
            solX = matching[0][1]

            print("signature point:", solX)
            break
    
if not found:
    print("Y has no signature point X")
    sys.exit(-1)

eqs_eval = tuple(eqs[i] - R(output[i]) for i in range(2 * n - 1))

field_eqs = vector(R,[X[i]**p-X[i] for i in range(2*n)])

magma_code = f'''
SetNthreads(15);
q:={p};
K:=FiniteField(q);
verbose:=1;
n:={2*n};
R<{str(X).strip('(').strip(')')}> := PolynomialRing(K, {2*n}, "grevlex");
eqs := {list(eqs_eval)};
I := ideal<R | eqs>;
SetVerbose("Groebner",verbose);
time GB:=GroebnerBasis(I);
'''
with open(f'systems/facto-X1-{p}-{2*n-1}-{2*n}.txt','w') as f:
    f.write(magma_code)

eqs_eval += tuple(field_eqs)

magma_code = f'''
SetNthreads(15);
q:={p};
K:=FiniteField(q);
verbose:=1;
n:={2*n};
R<{str(X).strip('(').strip(')')}> := PolynomialRing(K, {2*n}, "grevlex");
eqs := {list(eqs_eval)};
I := ideal<R | eqs>;
SetVerbose("Groebner",verbose);
time GB:=GroebnerBasis(I);
'''
with open(f'systems/facto-X1-feqs-{p}-{2*n-1}-{2*n}.txt','w') as f:
    f.write(magma_code)