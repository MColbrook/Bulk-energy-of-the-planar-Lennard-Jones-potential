"""Regenerate the signed reference in the original 20-coordinate system.

The basis, normalizations, fixed columns, and finite reference equations
come from OpenAI, An atomic certificate for triangular-lattice universal
optimality (2026), Sections 3--5.  The new input is the signed Lennard--Jones
Laplace weight.  Its moments and the two 20-by-20 inverse matrices are
enclosed with Arb at 160-bit precision.

Input: the rational lattice-sum enclosures from lattice_constants.py.
The analytic interpolation results are specified in ../MATHEMATICAL_SCOPE.txt.
"""
if not __debug__:
    raise RuntimeError("Certificate checks require assertions; do not use -O or -OO.")
import argparse
import hashlib
import json
from fractions import Fraction
from pathlib import Path
from flint import arb, acb, arb_mat, fmpq, ctx

ctx.prec = 160
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--lattice-constants", required=True, type=Path)
parser.add_argument("--output", required=True, type=Path)
args = parser.parse_args()


def rat(value):
    value = Fraction(value)
    return arb(fmpq(value.numerator, value.denominator))


def interval(lo, hi):
    lo, hi = Fraction(lo), Fraction(hi)
    return rat((lo + hi) / 2) + arb(0, rat((hi - lo) / 2).abs_upper())


def convolution(first, second):
    out = [acb(0) for _ in range(len(first) + len(second) - 1)]
    for i, x in enumerate(first):
        for j, y in enumerate(second):
            out[i + j] += x * y
    return out


pi, imaginary = arb.pi(), acb(0, 1)
b, h, H = arb(3).sqrt() / 2, rat("17/50"), rat("27/50")
kappa, lam = pi / 36, pi / b
residues = [0, 1, 3, 4, 7, 9, 12, 13, 16, 19, 21, 25, 27, 28, 31]
nodes = residues[1:11]
frequencies = [rat(Fraction(j - 15, 18)) for j in range(31)]
factors = [[-(imaginary * pi * a / 18).exp(), acb(2),
            -(-imaginary * pi * a / 18).exp()] for a in residues]
product = [acb(1)]
for factor in factors:
    product = convolution(product, factor)
columns = []
for a in range(11):
    base = [acb(1)]
    for r in range(15):
        if r != a:
            base = convolution(base, factors[r])
    columns.append([acb(0)] + [4 * kappa**2 * z for z in base] + [acb(0)])
    derivative = [imaginary * kappa * (imaginary * pi * residues[a] / 18).exp(),
                  acb(0),
                  -imaginary * kappa * (-imaginary * pi * residues[a] / 18).exp()]
    columns.append(convolution(derivative, base))
columns.append(product)

Q, D = {}, {}
for n in nodes:
    q, d = 4 * kappa**2, arb(0)
    for a in residues:
        if a != n:
            angle = kappa * (n - a)
            q *= 4 * angle.sin()**2
            d += 2 * kappa * angle.cos() / angle.sin()
    Q[n], D[n] = q, d

fourier = arb_mat(20, 23)
for ni, n in enumerate(nodes):
    for ci, column in enumerate(columns):
        damping = h if ci == 22 else H
        value, derivative = acb(0), acb(0)
        for j, frequency in enumerate(frequencies):
            denominator = acb(damping) - imaginary * frequency
            wave = (-pi * n / (b**2 * denominator) + pi * H * n).exp()
            wave /= b * denominator * Q[n]
            value -= column[j] * wave
            derivative += column[j] * wave * (
                D[n] + pi / (b**2 * denominator) - pi * H)
        fourier[2 * ni, ci] = value.real
        fourier[2 * ni + 1, ci] = derivative.real

R = arb_mat([[fourier[i, j + 2] for j in range(20)] for i in range(20)])
identity = arb_mat(20, 20)
for i in range(20):
    identity[i, i] = 1
inverse_plus, inverse_minus = (identity - R).inv(), (identity + R).inv()


def matrix_norm(matrix):
    return max(sum(abs(matrix[i, j]) for i in range(matrix.nrows())).upper()
               for j in range(matrix.ncols()))


norm_plus, norm_minus = matrix_norm(inverse_plus), matrix_norm(inverse_minus)
assert norm_plus < rat("16.16") and norm_minus < rat("3.89")
fixed_sum = [rat("22/25") * fourier[i, 0] for i in range(20)]
fixed_difference = [rat("6/500") * fourier[i, 22] for i in range(20)]

constants = json.loads(args.lattice_constants.read_text())
Z6 = interval(constants["Z6"]["lower_exact"], constants["Z6"]["upper_exact"]) * b**3
Z12 = interval(constants["Z12"]["lower_exact"], constants["Z12"]["upper_exact"]) * b**6
A, B = pi**6 / (1440 * Z12), pi**3 / (12 * Z6)


def moment(m, rate):
    return rate.gamma_upper(m + 1) / rate**(m + 1)


def integrated_reference(sign):
    positive_power, negative_power = (5, 2) if sign == "plus" else (-6, -3)

    def weighted_moment(shift, rate):
        return (A * moment(positive_power + shift, rate)
                - B * moment(negative_power + shift, rate))

    mass = (pi * h).exp() / Q[1] * (
        lam * weighted_moment(1, lam) - pi * h * weighted_moment(0, lam))
    source = [arb(0) for _ in range(20)]
    for j, n in enumerate(nodes[:4]):
        fac = (pi * H * n).exp() / Q[n]
        zeroth, first = weighted_moment(0, lam * n), weighted_moment(1, lam * n)
        source[2 * j] = fac * zeroth
        source[2 * j + 1] = fac * (-lam * first + (pi * H - D[n]) * zeroth)
    sum_solution = inverse_plus * arb_mat([
        [source[i] + mass * fixed_sum[i]] for i in range(20)])
    difference_solution = inverse_minus * arb_mat([
        [source[i] + mass * fixed_difference[i]] for i in range(20)])
    result = []
    for side in range(2):
        coefficients = [arb(0) for _ in range(23)]
        coefficients[0] = rat("11/25") * mass
        coefficients[22] = rat("-3/500" if side == 0 else "3/500") * mass
        for i in range(20):
            coefficients[i + 2] = (sum_solution[i, 0]
                + (1 if side == 0 else -1) * difference_solution[i, 0]) / 2
        result.append(coefficients)
    return result, mass


positive, Iplus = integrated_reference("plus")
negative, Iminus = integrated_reference("minus")
uref = [a - b for a, b in zip(positive[0], negative[1])]
vref = [a - b for a, b in zip(positive[1], negative[0])]
u_norm, v_norm = sum(abs(x) for x in uref), sum(abs(x) for x in vref)
assert u_norm < rat(".509") and v_norm < rat(".619")
alpha0 = (B / A)**rat("1/3")
assert B > A and A > 0 and alpha0 > 1


def truncated_moment(power):
    return (lam * alpha0).gamma_upper(power + 1) / lam**(power + 1)


positive_mass = (pi * h).exp() / Q[1] * (
    A * (lam * truncated_moment(6) - pi * h * truncated_moment(5))
    - B * (lam * truncated_moment(3) - pi * h * truncated_moment(2)))
absolute_mass = 2 * positive_mass - Iplus - Iminus
assert absolute_mass < rat(".425")
assert absolute_mass / 10**5 < rat("4.25e-6")
assert 3 * absolute_mass / 10**9 < rat("1.275e-9")
output = {
    "scope": "Finite reference and signed moments; imported infinite interpolation is an analytic premise",
    "finite_checks_passed": True,
    "precision_bits": ctx.prec,
    "source_input": args.lattice_constants.name,
    "source_input_sha256": hashlib.sha256(args.lattice_constants.read_bytes()).hexdigest(),
    "nodes": nodes,
    "A": str(A), "B": str(B),
    "Iplus": str(Iplus), "Iminus": str(Iminus),
    "slow_coefficient": str(vref[22]),
    "inverse_norm_plus": str(norm_plus), "inverse_norm_minus": str(norm_minus),
    "u_reference_coefficients": [str(x) for x in uref],
    "v_reference_coefficients": [str(x) for x in vref],
    "reference_u_l1": str(u_norm), "reference_v_l1": str(v_norm),
    "positive_weight_zero": str(alpha0),
    "absolute_mixture_mass": str(absolute_mass),
    "inherited_finite_error": str(absolute_mass / 10**5),
    "inherited_tail_error": str(3 * absolute_mass / 10**9),
}
destination = args.output
destination.parent.mkdir(parents=True, exist_ok=True)
destination.write_text(json.dumps(output, indent=2) + "\n")
print("20-coordinate inverse norms:", norm_plus, norm_minus)
print("Reference coefficient norms:", u_norm, v_norm)
print("Slow-column coefficient:", vref[22])
print("Absolute signed-mixture mass:", absolute_mass)
print("Generated", destination.name)
