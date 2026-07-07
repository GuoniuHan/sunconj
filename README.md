# A constructive approach to unit fractions with shifted prime denominators

**Guo-Niu Han** — 2026/07/08

## Abstract

A conjecture of Zhi-Wei Sun states that every positive rational number
can be written as a finite sum of unit fractions whose denominators are
of the form *p*+1 or *p*−1, where *p* runs over distinct primes.
Both cases have recently been proved in full generality by Thomas Bloom.
However, even for small integers such as *r*=2 or *r*=3,
finding an explicit decomposition is highly non-trivial: the required
primes can have hundreds of decimal digits, and the construction
involves hundreds of elimination steps.
We exhibit explicit decompositions for these two integer cases:
the case *r*=2 with denominators *p*+1, using a set
of 1963 distinct primes, and the case *r*=3 with denominators
*p*−1, using a set of 698 distinct primes.
The algorithm behind these constructions — based on beam search
and splitting tables — may be of independent interest for related
problems.

## The full paper

[sunconj.pdf](sunconj.pdf)

## The program

[sunconj.py](sunconj.py)

## Proofs of Theorems

We collect the outputs of the program `sunconj.py` for the Sun's conjecture
studied in the paper.

> **Note.** The extracted prime-set files `data_frac_*_1000.txt` are provided
> here gzip-compressed (`.gz`); run `gunzip <file>.gz` to obtain the plain
> text. The full raw outputs `output_frac_*_1000.txt` are very large and are
> not included in the repository, but can be regenerated with the commands
> below.

### Proof of Theorem 3

To reproduce the result of Theorem 3, run the following command:

```
python sunconj.py +1 2 false true > data_plus2.txt
```

This generates the output file [`data_plus2.txt`](data_plus2.txt),
from which the prime set stated in Theorem 3 is extracted and provided in
[`prime_set_P_plus2.txt`](prime_set_P_plus2.txt).

### Proof of Theorem 4

To reproduce the result of Theorem 4, run the following command:

```
python sunconj.py -1 3 false true > data_minus3.txt
```

This generates the output file [`data_minus3.txt`](data_minus3.txt),
from which the prime set stated in Theorem 4 is extracted and provided in
[`prime_set_P_minus3.txt`](prime_set_P_minus3.txt).

### Proof of Theorem 5.1 for `s = +1`

To reproduce the result of Theorem 5.1 for `s = +1`, run the following command
using the program [`sunconj_frac.py`](sunconj_frac.py):

```
python sunconj_frac.py 1000 +1 true false > output_frac_plus_1000.txt
```

This command generates the large output file `output_frac_plus_1000.txt`
(not included here); the prime sets extracted from it are provided in
[`data_frac_plus_1000.txt.gz`](data_frac_plus_1000.txt.gz).

> **Note.** The computation above (major bound `M = 1000`) takes a very long
> time. For a first try, run it with a small parameter, for example
> ```
> python sunconj_frac.py 100 +1 true false
> ```

### Proof of Theorem 5.1 for `s = -1`

To reproduce the result of Theorem 5.1 for `s = -1`, run the following command
using the program [`sunconj_frac.py`](sunconj_frac.py):

```
python sunconj_frac.py 1000 -1 true false > output_frac_minus_1000.txt
```

This command generates the large output file `output_frac_minus_1000.txt`
(not included here); the prime sets extracted from it are provided in
[`data_frac_minus_1000.txt.gz`](data_frac_minus_1000.txt.gz).

> **Note.** The computation above (major bound `M = 1000`) takes a very long
> time. For a first try, run it with a small parameter, for example
> ```
> python sunconj_frac.py 100 +1 true false
> ```

### Proof of Theorem 5.2 for `s = +1`

To reproduce the result of Theorem 5.2 for `s = +1`, run the following command
using the program [`sunconj_prime.py`](sunconj_prime.py):

```
python sunconj_prime.py 100000 +1 true true > output_prime_plus_100000.txt
```

This generates the output file [`output_prime_plus_100000.txt`](output_prime_plus_100000.txt),
from which the prime sets are extracted and provided in
[`data_prime_plus_100000.txt`](data_prime_plus_100000.txt).

> **Note.** The computation above (major bound `M = 100000`) takes a very long
> time. For a first try, run it with a small parameter, for example
> ```
> python sunconj_prime.py 1000 +1 true true
> ```

### Proof of Theorem 5.2 for `s = -1`

To reproduce the result of Theorem 5.2 for `s = -1`, run the following command
using the program [`sunconj_prime.py`](sunconj_prime.py):

```
python sunconj_prime.py 100000 -1 true true > output_prime_minus_100000.txt
```

This generates the output file [`output_prime_minus_100000.txt`](output_prime_minus_100000.txt),
from which the prime sets are extracted and provided in
[`data_prime_minus_100000.txt`](data_prime_minus_100000.txt).

> **Note.** The computation above (major bound `M = 100000`) takes a very long
> time. For a first try, run it with a small parameter, for example
> ```
> python sunconj_prime.py 1000 +1 true true
> ```

### Proof of Theorem 5.3 for `s = +1`

To reproduce the result of Theorem 5.3 for `s = +1`, run the following command
using the program [`sunconj_power.py`](sunconj_power.py):

```
python sunconj_power.py 40 +1 true true > output_power_plus_40.txt
```

This generates the output file [`output_power_plus_40.txt`](output_power_plus_40.txt),
from which the prime sets are extracted and provided in
[`data_power_plus_40.txt`](data_power_plus_40.txt).

### Proof of Theorem 5.3 for `s = -1`

To reproduce the result of Theorem 5.3 for `s = -1`, run the following command
using the program [`sunconj_power.py`](sunconj_power.py):

```
python sunconj_power.py 40 -1 true true > output_power_minus_40.txt
```

This generates the output file [`output_power_minus_40.txt`](output_power_minus_40.txt),
from which the prime sets are extracted and provided in
[`data_power_minus_40.txt`](data_power_minus_40.txt).

## Conjecture 7 (length of reciprocal decompositions)

The paper conjectures that every reciprocal `1/d` (`d ≥ 2`) admits a
shifted-prime decomposition `1/d = 1/(p_1 + s) + ... + 1/(p_k + s)` into
distinct primes whose length satisfies `|P| = k ≤ log2(d) + 1`, for both
shifts `s ∈ {+1, -1}`.

The self-contained program [`sunconj_reciprocal.py`](sunconj_reciprocal.py)
checks this conjecture. For each `d` from 2 to `M` it constructs a
decomposition of `1/d`, verifies the identity in exact rational arithmetic,
and checks the length bound `|P| ≤ log2(d) + 1`, printing any violation and a
final summary. It embeds the decomposition algorithm, so it does **not**
depend on `sunconj.py` at run time (only `gmpy2` is required).

```
python sunconj_reciprocal.py 3000 +1
python sunconj_reciprocal.py 3000 -1
```

The first argument is the range bound `M`; the second is the shift
`s ∈ {+1, -1}` (default `+1`).

> **Note.** Checking up to `M = 3000` takes a very long time. For a first try,
> run it with a small parameter, for example
> ```
> python sunconj_reciprocal.py 100 +1
> ```
