"""
sunconj_reciprocal.py — self-contained checker for Conjecture 7
================================================================

Conjecture 7 (length bound for reciprocals).  For every integer d >= 2 and
shift s in {+1, -1}, the reciprocal 1/d admits a decomposition

        1/d = 1/(p_1 + s) + ... + 1/(p_k + s),     p_1 < ... < p_k  prime,

into DISTINCT shifted-prime unit fractions whose length satisfies

        |P| = k  <=  log2(d) + 1.

This program is self-contained: it embeds the decomposition algorithm of
sunconj.py (greedy + beam search + splitting tables + deduplication), so it
does NOT call sunconj.py as a subprocess and prints no RESGOOD lines.  For
each d it finds a decomposition of 1/d, verifies the identity in exact
rational arithmetic, and checks the length bound, reporting any violation.

Usage
-----
    python sunconj_reciprocal.py M [SHIFT]
    python sunconj_reciprocal.py 3000 +1

Requires: gmpy2.
"""

import sys
import math
import gmpy2
from gmpy2 import mpz, mpq, is_prime

sys.setrecursionlimit(10_000)

# ---------------------------------------------------------------------------
# Global configuration (SHIFT is set by the checker before each run)
# ---------------------------------------------------------------------------

SHIFT: int = +1          # +1 -> 1/(p+1),  -1 -> 1/(p-1)

# Splitting mode used by the decomposer.  prefer_large=False keeps the largest
# prime small (it avoids the rare pathological blow-ups seen with True) while
# giving the same short lengths, so it is the value used throughout.
PREFER_LARGE: bool = False

_split_cache: dict = {}


# ---------------------------------------------------------------------------
# Utility helpers
# ---------------------------------------------------------------------------

def unit_fraction(p) -> mpq:
    """Return the exact rational 1 / (p + SHIFT)."""
    return mpq(1, int(p) + SHIFT)


def prime_sum(primes: list) -> mpq:
    """Exact sum  Sum 1/(p + SHIFT)  over p in *primes*."""
    total = mpq(0)
    for p in primes:
        total += unit_fraction(p)
    return total


def multiplicity_dict(primes: list) -> dict:
    counts: dict = {}
    for p in primes:
        counts[p] = counts.get(p, 0) + 1
    return counts


def count_duplicates(primes: list) -> int:
    """Total number of excess occurrences; 0 means every prime is distinct."""
    excess = 0
    for cnt in multiplicity_dict(primes).values():
        if cnt >= 2:
            excess += cnt - 1
    return excess


# ---------------------------------------------------------------------------
# Split tables  (ways to write 1 as a sum of 2 or 3 unit fractions)
# ---------------------------------------------------------------------------

def _reduced_fractions_up_to(bound: int) -> set:
    fracs = set()
    for k in range(2, bound):
        for j in range(1, k):
            if gmpy2.gcd(j, k) == 1:
                fracs.add(mpq(j, k))
    return fracs


def _build_one_as_two(bound: int) -> list:
    fracs = _reduced_fractions_up_to(bound)
    pairs = []
    for a in fracs:
        if a < mpq(1, 2) and (1 - a) in fracs:
            pairs.append([a, 1 - a])
    return pairs


def _build_one_as_three(bound: int) -> list:
    fracs = _reduced_fractions_up_to(bound)
    frac_list = sorted(f for f in fracs if f < mpq(1, 2))
    triples = []
    for a in frac_list:
        if a >= mpq(1, 3):
            break
        for b in frac_list:
            if b <= a:
                continue
            c = 1 - a - b
            if c in fracs and c > b:
                triples.append([a, b, c])
    return triples


_SPLIT_TABLES: dict = {}

def _init_split_tables():
    global _SPLIT_TABLES
    _SPLIT_TABLES["Fast"]   = [[mpq(1, 1)]] + _build_one_as_two(10) + _build_one_as_three(7)
    _SPLIT_TABLES["Middle"] = [[mpq(1, 1)]] + _build_one_as_two(20) + _build_one_as_three(13)
    _SPLIT_TABLES["Deep"]   = [[mpq(1, 1)]] + _build_one_as_two(30) + _build_one_as_three(17)

_init_split_tables()


# ---------------------------------------------------------------------------
# Core greedy prime search
# ---------------------------------------------------------------------------

def smallest_feasible_prime(remainder: mpq, min_prime):
    """Smallest prime p > min_prime with 1/(p+SHIFT) <= remainder, or None."""
    if remainder <= 0:
        return None
    p_lb = max(int(1 / remainder) - 2, 1)
    p_lb = max(p_lb, int(min_prime) + 1)
    p = mpz(p_lb)
    while True:
        if is_prime(p) and remainder >= unit_fraction(p):
            return p
        p += 1


def greedy_decompose(remainder: mpq, min_prime=1, cutoff=10**8, max_steps=99):
    """Greedy: subtract the largest feasible 1/(p+SHIFT) until 0. list | False."""
    current = remainder
    last_prime = min_prime
    primes_found = []
    for _ in range(max_steps):
        p = smallest_feasible_prime(current, last_prime)
        if p is None:
            return False
        primes_found.append(int(p))
        current -= unit_fraction(p)
        last_prime = p
        if current == 0:
            return primes_found
        if p > cutoff:
            if cutoff < 10**30:
                return False
            tail = greedy_decompose(current, min_prime=int(last_prime),
                                    cutoff=cutoff**2, max_steps=max_steps)
            if tail is False:
                return False
            return primes_found[:-1] + [primes_found[-1]] + tail
    return False


# ---------------------------------------------------------------------------
# Beam / candidate expansion
# ---------------------------------------------------------------------------

def _num_candidates(mode: str) -> int:
    return {"Fast": 3, "Middle": 4, "Deep": 6}[mode]


def generate_prime_candidates(remainder: mpq, min_prime, mode: str = "Fast") -> list:
    n = _num_candidates(mode)
    candidates = []
    prev_p = min_prime
    for _ in range(n):
        p = smallest_feasible_prime(remainder, prev_p)
        if p is None:
            break
        candidates.append([[p], remainder - unit_fraction(p)])
        prev_p = p
    return candidates


def beam_search_decompose(remainder: mpq, min_prime=1, max_steps=99,
                          mode: str = "Fast"):
    """Structured beam search for a decomposition of *remainder*. list | False."""
    expansion_rounds = 1 if mode == "Fast" else 3
    frontier = generate_prime_candidates(remainder, min_prime, mode)
    solutions = []
    for _ in range(expansion_rounds):
        next_frontier = []
        for node_primes, node_rem in frontier:
            if node_rem == 0:
                solutions.append(node_primes)
                continue
            last_p = node_primes[-1]
            children = generate_prime_candidates(node_rem, last_p, mode)
            for child_primes, child_rem in children:
                merged = sorted(child_primes + node_primes)
                next_frontier.append([merged, child_rem])
        frontier = next_frontier
    for node_primes, node_rem in frontier:
        if node_rem == 0:
            solutions.append(node_primes)
            continue
        tail = greedy_decompose(node_rem, min_prime=node_primes[-1])
        if tail is not False:
            solutions.append(node_primes + tail)
    if not solutions:
        return False
    return min(solutions, key=sum)


# ---------------------------------------------------------------------------
# Split-table based decomposition
# ---------------------------------------------------------------------------

def _decompose_via_split_table(target: mpq, min_prime: int, mode: str):
    table = _SPLIT_TABLES[mode]
    best_partial = None
    for weights in table:
        sub_targets = [w * target for w in weights]
        results = [beam_search_decompose(st, min_prime, mode=mode)
                   for st in sub_targets]
        failed = results.count(False)
        if failed == 0:
            combined = []
            for res in results:
                combined.extend(res)
            return combined
        if failed == 1 and best_partial is None:
            todo = [sub_targets[i] for i, r in enumerate(results) if r is False]
            solved = [p for res in results if res for p in res]
            best_partial = ["PARTIAL", todo, solved]
    if best_partial and best_partial[2]:
        return best_partial
    return False


def _decompose_list(targets: list, min_prime: int, mode: str):
    todo: list = []
    solved: list = []
    for t in targets:
        result = _decompose_via_split_table(t, min_prime, mode)
        if result is False:
            return False
        if isinstance(result, list) and result[0] == "PARTIAL":
            todo.extend(result[1])
            solved.extend(result[2])
        else:
            solved.extend(result)
    return ["PARTIAL", todo, solved] if todo else solved


def decompose_rational(target: mpq, min_prime: int, mode: str = "Fast"):
    pending = [target]
    all_solved: list = []
    for _ in range(5):
        result = _decompose_list(pending, min_prime, mode)
        if result is False:
            return False
        if not (isinstance(result, list) and result[0] == "PARTIAL"):
            all_solved.extend(result)
            return sorted(all_solved)
        pending = result[1]
        all_solved.extend(result[2])
    return False


# ---------------------------------------------------------------------------
# Prime splitting  (replace m copies of 1/(p+SHIFT) by distinct primes)
# ---------------------------------------------------------------------------

def _split_prime_multiplicity(prime: int, multiplicity: int, mode: str = "Fast"):
    global _split_cache
    for kept in range(multiplicity, 0, -1):
        residual = mpq(kept, int(prime) + SHIFT)
        result = decompose_rational(residual, prime, mode=mode)
        if result is not False:
            _split_cache[(prime, kept)] = result
            extra_copies = [prime] * (multiplicity - kept)
            return sorted(result + extra_copies)
    return False


def _lookup_split_cache(prime: int, multiplicity: int):
    if (prime, multiplicity) in _split_cache:
        return sorted(_split_cache[(prime, multiplicity)])
    return False


def expand_repeated_prime(prime: int, multiplicity: int, mode: str = "Fast"):
    cached = _lookup_split_cache(prime, multiplicity)
    if cached is not False:
        return cached
    return _split_prime_multiplicity(prime, multiplicity, mode=mode)


# ---------------------------------------------------------------------------
# Deduplication loop  (make all primes distinct)
# ---------------------------------------------------------------------------

def reduce_duplicates(primes: list, mode: str = "Fast"):
    L = sorted(primes)
    repeated = None
    for i in range(len(L) - 1):
        if L[i] == L[i + 1]:
            repeated = L[i]
            break
    if repeated is None:
        return L
    count = L.count(repeated)
    for _ in range(count - 1):
        L.remove(repeated)
    replacement = expand_repeated_prime(repeated, count - 1, mode=mode)
    if replacement is False:
        return False
    L.extend(replacement)
    return sorted(L)


# ---------------------------------------------------------------------------
# Top-level decomposition (silent: returns a sorted list of distinct primes,
# or False).  No printing, no RESGOOD.
# ---------------------------------------------------------------------------

def decompose(target: mpq, prefer_large: bool = False):
    primes = decompose_rational(target, 1, mode="Fast")
    if primes is False:
        primes = decompose_rational(target, 1, mode="Middle")
    if primes is False:
        primes = decompose_rational(target, 1, mode="Deep")
    if primes is False:
        return False

    split_mode = "Middle" if prefer_large else "Fast"
    for _ in range(10_000):
        if count_duplicates(primes) == 0:
            break
        updated = reduce_duplicates(primes, mode=split_mode)
        if updated is False:
            updated = reduce_duplicates(primes, mode="Middle")
        if updated is False:
            updated = reduce_duplicates(primes, mode="Deep")
        if updated is False:
            break
        primes = updated

    if count_duplicates(primes) > 0:
        return False
    return sorted(set(int(p) for p in primes))


# ---------------------------------------------------------------------------
# Conjecture 7 checker
# ---------------------------------------------------------------------------

def length_bound(d: int) -> float:
    """Conjectured length bound: |P| <= log2(d) + 1."""
    return math.log2(d) + 1.0


def check_reciprocals(M, shift=1):
    global SHIFT
    SHIFT = shift

    failures = []      # d with no valid decomposition found
    violations = []    # (d, |P|, bound) with |P| > log2(d)+1
    max_len = 0
    max_len_d = None
    max_slack = -1e18
    max_slack_d = None

    for d in range(2, M + 1):
        sol = decompose(mpq(1, d), prefer_large=PREFER_LARGE)

        # --- validity checks (exact arithmetic) ---
        if sol is False:
            failures.append(d)
            print(f"d={d}: FAILED (no decomposition found)", flush=True)
            continue
        if prime_sum(sol) != mpq(1, d):
            failures.append(d)
            print(f"d={d}: FAILED (sum != 1/{d})", flush=True)
            continue
        if len(set(sol)) != len(sol) or not all(is_prime(mpz(p)) for p in sol):
            failures.append(d)
            print(f"d={d}: FAILED (non-distinct or non-prime entry)", flush=True)
            continue

        # --- length bound check ---
        l = len(sol)
        lb = length_bound(d)
        ok = (l <= lb)
        if l > max_len:
            max_len, max_len_d = l, d
        slack = l - lb
        if slack > max_slack:
            max_slack, max_slack_d = slack, d
        if not ok:
            violations.append((d, l, round(lb, 3)))

        if not ok:
            print(f"d={d}: |P|={l}  (bound log2 d + 1 = {lb:.2f})   <-- VIOLATION",
                  flush=True)

    # ------------------------------- summary ------------------------------
    print("\n" + "=" * 64)
    print(f"Conjecture 7 check:  1/d = sum 1/(p+s),  d = 2..{M},  s = {shift:+d}")
    print(f"  bound tested:            |P| <= log2(d) + 1")
    print(f"  values tested:           {M - 1}")
    print(f"  decomposition failures:  {len(failures)}"
          + (f"  -> {failures}" if failures else ""))
    print(f"  length-bound violations: {len(violations)}"
          + (f"  -> {violations}" if violations else ""))
    print(f"  longest decomposition:   |P| = {max_len} at d = {max_len_d}"
          f"  (bound {length_bound(max_len_d):.2f})" if max_len_d else "")
    print(f"  tightest case:           slack {max_slack:+.3f} at d = {max_slack_d}")
    if not failures and not violations:
        print("  RESULT: Conjecture 7 holds for every tested d.")
    else:
        print("  RESULT: see failures/violations listed above.")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python sunconj_reciprocal.py M [SHIFT]")
        print("Example: python sunconj_reciprocal.py 3000 +1")
        sys.exit(1)

    M = int(sys.argv[1])
    shift = int(sys.argv[2]) if len(sys.argv) > 2 else 1

    if shift not in (+1, -1):
        print("Error: SHIFT must be +1 or -1")
        sys.exit(1)

    check_reciprocals(M, shift)
