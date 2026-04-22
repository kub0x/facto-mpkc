#!/usr/bin/env python3
from sage.all import *
from itertools import combinations

p = 11
n = 7
k=n-1
F = GF(p)

def random_invertible_matrix(F, r):
    while True:
        M = random_matrix(F, r, r)
        if M.is_invertible():
            return M

def multiplication_matrix(F, n):
    """
    Build M of size (2n) x n^2 such that for u,v in F^n,

        M * (u ⊗ v)

    gives the coefficient vector of u(t)*v(t).

    Tensor order:
        u ⊗ v = [u0*v0, u0*v1, ..., u0*v_{n-1},
                 u1*v0, u1*v1, ..., u1*v_{n-1},
                 ...
                 u_{n-1}*v_{n-1}]
    """
    rows = 2 * n
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

S = random_invertible_matrix(F, 2 * n)

invS = S.inverse()

S1 = S[0:n, :]
S2 = S[n:2*n, :]

L = block_matrix([ [zero_matrix(n,n),zero_matrix(n,n)],[identity_matrix(n),identity_matrix(n)] ])

T = random_invertible_matrix(F, 2 * n)

M = multiplication_matrix(F, n)

K = S1.tensor_product(S2)

Q1 = M * K

Pmat = T * block_matrix([[Q1, L*S]])

XX = vector(tensor_vec_symbolic(X, X).list() + X.list())

P_of_X = Pmat * XX 

eqs = []
for i, eq in enumerate(P_of_X):
    eqs.append(eq)
    
in_x = random_vector(GF(p), 2*n)
print("Random input X0 =")
print(in_x)
print()

eval_x = {f'{vars[i]}' : in_x[i]  for i in range(2*n)}

eval_P = P_of_X.subs(**eval_x)

output = vector(eval_P)
print("Output Y0 =")
print(output)
print()

out2 = T.inverse()*output
outpoly = sum(out2[i] * t**i for i in range(2*n)) + t**(2*n)
factorization = outpoly.factor()
terms = [f.degree() * e for (f, e) in factorization]


partition =[]
irred_factors = []
for i in range(len(terms)): 
    for j in range(factorization[i][1]):
        partition.append(factorization[i][0].degree())
        irred_factors.append(factorization[i][0])

print("original partition", terms)
print("extended partition", partition)
seen = set()

found = False

for solution_indices in find_subset_sum(partition, n):
    if solution_indices not in seen:
        seen.add(solution_indices) 

        solution1 = 1
        solution2 = 1

        for i in range(len(solution_indices)):
            solution1 *= irred_factors[solution_indices[i]]

        for i in range(len(partition)):
            if (i not in solution_indices):
                solution2 *= irred_factors[i]

        assert(solution1*solution2 == outpoly)
        
        sol1 = extract_coeffs_in_var_as_vector(solution1, t, n)
        sol2 = extract_coeffs_in_var_as_vector(solution2, t, n)
        xx = invS*vector(sol1.list() + sol2.list())

        if (in_x == xx or in_x == -xx):
            print("solution found 1")
            found = True
            break
        
        xx = invS*vector(sol2.list() + sol1.list())

        if (in_x == xx or in_x == -xx):
            print("solution found 2")
            found = True
            break

if not found:
    print("Y has no plaintext point X")
    sys.exit(-1)
    
eqs_eval = tuple(eqs[i] - R(output[i]) for i in range(2 * n))

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
with open(f'systems/facto-EDF-{p}-{2*n}-{2*n}.txt','w') as f:
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
with open(f'systems/facto-EDF-feqs-{p}-{2*n}-{2*n}.txt','w') as f:
    f.write(magma_code)