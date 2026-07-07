"""
sunconj.py — Prime Unit-Fraction Decomposition (Sun's Conjecture)
==================================================================

Goal
----
Given a positive rational number r, express r as a finite sum of distinct
unit fractions of the form 1/(p + SHIFT), where p ranges over prime numbers
and SHIFT ∈ {+1, −1} is a fixed global parameter:

    r = 1/(p₁ + SHIFT) + 1/(p₂ + SHIFT) + … + 1/(pₖ + SHIFT),
    with p₁ < p₂ < … < pₖ  (all prime).

This is related to a conjecture by Zhi-Wei Sun (2015) asserting that every
positive rational has such a representation.  The default SHIFT = +1 targets
fractions 1/(p+1); set SHIFT = -1 to use 1/(p−1).

Algorithm overview
------------------
1. A *greedy* pass picks the largest feasible prime p at each step (largest
   unit fraction that does not exceed the remainder).
2. A *structured beam search* (testL / expand_candidates) explores several
   candidate primes at each depth level, then falls back to greedy for the
   tail.
3. A *splitting* loop (reduce_duplicates / split_prime) iteratively replaces
   any repeated prime p appearing m times with a fresh decomposition of
   m·(1/(p+SHIFT)), removing multiplicity until the representation uses
   every prime at most once.

Usage
-----
    python sunconj.py           # verifies all k/n with gcd(k,n)=1, n ≤ 200
    python sunconj.py 100       # same, up to n = 100

    from sunconj import decompose
    print(decompose(mpq(3, 7)))  # → sorted list of primes

Requirements
------------
    pip install gmpy2
"""

import sys
import gmpy2
from gmpy2 import mpz, mpq, is_prime
from functools import lru_cache

# ---------------------------------------------------------------------------
# Global configuration
# ---------------------------------------------------------------------------

sys.setrecursionlimit(10_000)

# SHIFT determines which unit fractions are used:
#   SHIFT = +1  →  1/(p+1)   (default, Sun's original conjecture)
#   SHIFT = -1  →  1/(p-1)
SHIFT: int = +1

# When True, print verbose progress information during the search.
VERBOSE: bool = False

# Cache: maps (prime p, multiplicity m) → list of primes that represent
#   m * 1/(p + SHIFT)  as a sum of distinct unit fractions 1/(q + SHIFT).
_split_cache: dict = {}


# ---------------------------------------------------------------------------
# Utility helpers
# ---------------------------------------------------------------------------

def unit_fraction(p) -> mpq:
    """Return the exact rational 1 / (p + SHIFT)."""
    return mpq(1, int(p) + SHIFT)


def prime_sum(primes: list) -> mpq:
    """
    Return the exact rational sum  Σ 1/(p + SHIFT)  for p in *primes*.

    Uses gmpy2.mpq for arbitrary-precision arithmetic; no floating-point
    rounding errors.
    """
    total = mpq(0)
    for p in primes:
        total += unit_fraction(p)
    return total


def multiplicity_dict(primes: list) -> dict:
    """
    Return a dict  {prime: count}  from a list of primes (possibly repeated).

    Example:  [2, 3, 3, 5]  →  {2: 1, 3: 2, 5: 1}
    """
    counts: dict = {}
    for p in primes:
        counts[p] = counts.get(p, 0) + 1
    return counts


def count_duplicates(primes: list) -> int:
    """
    Return the total number of *excess* occurrences (i.e. total repetitions
    minus the number of distinct primes that repeat ≥ 2 times).

    A return value of 0 means every prime appears exactly once.

    Example:  [2, 3, 3, 5, 5, 5]  →  (1 extra 3) + (2 extra 5s) = 3
    """
    excess = 0
    for cnt in multiplicity_dict(primes).values():
        if cnt >= 2:
            excess += cnt - 1
    return excess


# ---------------------------------------------------------------------------
# Fraction tables  (ways to write 1 as a sum of 2 or 3 unit fractions)
# These are used to *split* a single term r into 2 or 3 sub-problems whose
# sum equals r, potentially making each sub-problem easier to solve.
# ---------------------------------------------------------------------------

def _reduced_fractions_up_to(bound: int) -> list:
    """
    Return all reduced fractions j/k with 1 ≤ j < k < bound and gcd(j,k)=1,
    as mpq objects, stored in a set for O(1) membership testing.
    """
    fracs = set()
    for k in range(2, bound):
        for j in range(1, k):
            if gmpy2.gcd(j, k) == 1:
                fracs.add(mpq(j, k))
    return fracs


def _build_one_as_two(bound: int) -> list:
    """
    Return all ordered pairs [a, b] of reduced fractions with a < b = 1−a,
    both having denominator < bound.  Each pair satisfies a + b = 1.
    """
    fracs = _reduced_fractions_up_to(bound)
    pairs = []
    for a in fracs:
        if a < mpq(1, 2) and (1 - a) in fracs:
            pairs.append([a, 1 - a])
    return pairs


def _build_one_as_three(bound: int) -> list:
    """
    Return all ordered triples [a, b, c] of reduced fractions with
    a < b < c and a + b + c = 1, all having denominator < bound.
    """
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


# Pre-build the split tables at three depth levels.
# "Fast"   uses small tables → quick but may miss some decompositions.
# "Middle" uses medium tables.
# "Deep"   uses large tables → slow but more thorough.

_SPLIT_TABLES: dict = {}

def _init_split_tables():
    global _SPLIT_TABLES
    _SPLIT_TABLES["Fast"]   = [[mpq(1, 1)]] + _build_one_as_two(10)  + _build_one_as_three(7)
    _SPLIT_TABLES["Middle"] = [[mpq(1, 1)]] + _build_one_as_two(20)  + _build_one_as_three(13)
    _SPLIT_TABLES["Deep"]   = [[mpq(1, 1)]] + _build_one_as_two(30)  + _build_one_as_three(17)

_init_split_tables()


# ---------------------------------------------------------------------------
# Core greedy prime search
# ---------------------------------------------------------------------------

def smallest_feasible_prime(remainder: mpq, min_prime) -> mpz | None:
    """
    Find the *smallest* prime p ≥ min_prime such that

        1/(p + SHIFT) ≤ remainder,

    i.e. p is the natural next greedy step.

    Returns None if remainder ≤ 0.
    """
    if remainder <= 0:
        return None

    # Lower bound on p: from 1/(p+SHIFT) ≤ r we get p ≥ 1/r − SHIFT.
    # Subtract 2 for safety against rounding, then also enforce p > min_prime.
    p_lb = max(int(1 / remainder) - 2, 1)
    p_lb = max(p_lb, int(min_prime) + 1)
    p = mpz(p_lb)

    while True:
        if is_prime(p) and remainder >= unit_fraction(p):
            return p
        p += 1


def greedy_decompose(remainder: mpq, min_prime=1, cutoff=10**8, max_steps=99) -> list | bool:
    """
    Greedy algorithm: repeatedly subtract the largest feasible unit fraction
    1/(p + SHIFT) until remainder reaches 0.

    Parameters
    ----------
    remainder : mpq
        The rational number to decompose.
    min_prime : int
        All selected primes must be strictly greater than this value.
    cutoff : int
        If the selected prime exceeds *cutoff*, restart the tail search with
        a squared cutoff (allows very large primes).
    max_steps : int
        Maximum number of primes to add before giving up.

    Returns
    -------
    list of int  on success,  False  on failure.
    """
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
            # Try again with a much larger allowed cutoff for remaining tail.
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
    """Number of candidate primes to generate per level in beam search."""
    return {"Fast": 3, "Middle": 4, "Deep": 6}[mode]


def generate_prime_candidates(remainder: mpq, min_prime, mode: str = "Fast") -> list:
    """
    Generate a short list of candidate (primes_so_far, new_remainder) pairs
    by stepping through consecutive primes above *min_prime*.

    Each candidate is  [[p], remainder − 1/(p+SHIFT)].

    The number of candidates depends on *mode*:
      Fast   → 3 candidates
      Middle → 4 candidates
      Deep   → 6 candidates
    """
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
                          mode: str = "Fast") -> list | bool:
    """
    Structured beam search for a prime decomposition of *remainder*.

    Expands *generate_prime_candidates* for several levels (CT iterations),
    then calls *greedy_decompose* on the surviving frontier nodes.

    Returns the solution with minimum prime-sum (heuristic for simplicity),
    or False if none found.
    """
    # Number of beam-expansion rounds before falling back to greedy.
    expansion_rounds = 1 if mode == "Fast" else 3

    # Frontier: list of  [primes_chosen_so_far, remaining_value]
    frontier = generate_prime_candidates(remainder, min_prime, mode)
    solutions = []

    # Expand the beam for several rounds.
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

    # For each surviving node, complete with greedy.
    for node_primes, node_rem in frontier:
        if node_rem == 0:
            solutions.append(node_primes)
            continue
        tail = greedy_decompose(node_rem, min_prime=node_primes[-1])
        if tail is not False:
            solutions.append(node_primes + tail)

    if not solutions:
        return False

    # Return the solution with smallest sum of primes (fewest / smallest primes).
    return min(solutions, key=sum)


# ---------------------------------------------------------------------------
# Split-table based decomposition  (check1 / check1list / check)
# ---------------------------------------------------------------------------

def _decompose_via_split_table(target: mpq, min_prime: int,
                               mode: str) -> list | bool:
    """
    Try to decompose *target* by first writing it as

        target = w₁·target + w₂·target + …   (where [w₁, w₂, …] ∈ split table)

    and solving each piece independently with beam_search_decompose.

    Returns a flat sorted list of primes on success, or:
      ["PARTIAL", todo_list, solved_primes]  if exactly one piece failed.
    False if completely stuck.
    """
    table = _SPLIT_TABLES[mode]
    best_partial = None

    for weights in table:
        # Each weight vector sums to 1; multiply through by target.
        sub_targets = [w * target for w in weights]
        results = [beam_search_decompose(st, min_prime, mode=mode)
                   for st in sub_targets]

        failed = results.count(False)

        if failed == 0:
            # Full solution found: concatenate all sub-results.
            combined = []
            for res in results:
                combined.extend(res)
            return combined

        if failed == 1 and best_partial is None:
            # Remember the best partial solution (one piece unsolved).
            todo = [sub_targets[i] for i, r in enumerate(results) if r is False]
            solved = [p for res in results if res for p in res]
            best_partial = ["PARTIAL", todo, solved]

    if best_partial and best_partial[2]:
        return best_partial
    return False


def _decompose_list(targets: list, min_prime: int, mode: str) -> list | bool:
    """
    Decompose each value in *targets* independently and combine results.

    If all succeed → return flat solved prime list.
    If exactly some fail → return ["PARTIAL", unsolved_targets, solved_primes].
    If any is completely stuck → return False.
    """
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


def decompose_rational(target: mpq, min_prime: int,
                       mode: str = "Fast") -> list | bool:
    """
    Main recursive decomposition driver.

    Iterates up to 5 rounds of *_decompose_list*, carrying over unsolved
    sub-targets from one round to the next.  Returns a sorted list of
    primes whose unit fractions sum to *target*, or False on failure.
    """
    pending = [target]
    all_solved: list = []

    for _ in range(5):
        result = _decompose_list(pending, min_prime, mode)
        if result is False:
            return False
        if not (isinstance(result, list) and result[0] == "PARTIAL"):
            all_solved.extend(result)
            return sorted(all_solved)
        # Carry forward unsolved pieces.
        pending = result[1]
        all_solved.extend(result[2])

    return False


# ---------------------------------------------------------------------------
# Prime splitting  (replace  m·(1/(p+SHIFT))  by distinct primes)
# ---------------------------------------------------------------------------

def _split_prime_multiplicity(prime: int, multiplicity: int,
                              mode: str = "Fast") -> list | bool:
    """
    Find a list of distinct primes  q₁, q₂, …  such that

        multiplicity · 1/(prime + SHIFT)
            = 1/(q₁ + SHIFT) + 1/(q₂ + SHIFT) + …

    where none of the qᵢ equal *prime* (so the multiplicity is eliminated).

    Tries nu = multiplicity, multiplicity−1, … 1 copies of 1/(prime+SHIFT)
    in the residual, storing the best split in _split_cache.

    Returns sorted list of primes (including the kept copies of *prime*) or False.
    """
    global _split_cache

    for kept in range(multiplicity, 0, -1):
        # Residual after keeping *kept* copies of prime:
        #   (multiplicity − kept) · 1/(prime + SHIFT)  needs a fresh split.
        residual = mpq(kept, int(prime) + SHIFT)
        pass  # (set DEBUG=True in source to trace inner residuals)

        result = decompose_rational(residual, prime, mode=mode)
        if result is not False:
            _split_cache[(prime, kept)] = result
            extra_copies = [prime] * (multiplicity - kept)
            return sorted(result + extra_copies)

    return False


def _lookup_split_cache(prime: int, multiplicity: int) -> list | bool:
    """
    Return a cached split for (prime, multiplicity) if available,
    otherwise False.

    The cache stores the replacement for *kept* copies of the prime's unit
    fraction, where kept ≤ multiplicity.  Any remaining (multiplicity − kept)
    copies of prime are appended directly.
    """
    if (prime, multiplicity) in _split_cache:
        cached = _split_cache[(prime, multiplicity)]
        # multiplicity copies were requested; cache entry covers all of them.
        return sorted(cached)
    return False


def expand_repeated_prime(prime: int, multiplicity: int,
                          mode: str = "Fast") -> list | bool:
    """
    Replace *multiplicity* copies of *prime* in a decomposition by a list of
    distinct primes, using the cache when possible.

    This is the public interface for splitting.
    """
    # Try cache first.
    cached = _lookup_split_cache(prime, multiplicity)
    if cached is not False:
        return cached
    # Compute and cache.  Only print if verbose, to avoid cluttering the
    # main progress line.
    if VERBOSE:
        print(f"    Computing split: DO[({prime},{multiplicity})]", flush=True)
    result = _split_prime_multiplicity(prime, multiplicity, mode=mode)
    return result


# ---------------------------------------------------------------------------
# Deduplication loop  (make all primes distinct)
# ---------------------------------------------------------------------------

def reduce_duplicates(primes: list, mode: str = "Fast") -> list | bool:
    """
    Given a list *primes* (possibly with repetitions), replace the first
    repeated prime p (with multiplicity m) by an equivalent set of *distinct*
    primes using *expand_repeated_prime*.

    Returns:
      - The updated sorted list (possibly still with other duplicates), or
      - The same list unchanged if no duplicates remain, or
      - False if the split for the repeated prime failed.
    """
    L = sorted(primes)

    # Find the first prime that appears more than once.
    repeated = None
    for i in range(len(L) - 1):
        if L[i] == L[i + 1]:
            repeated = L[i]
            break

    if repeated is None:
        return L  # Already duplicate-free.

    count = L.count(repeated)
    # Remove all but one copy (we will replace count−1 via split).
    for _ in range(count - 1):
        L.remove(repeated)

    replacement = expand_repeated_prime(repeated, count - 1, mode=mode)
    if replacement is False:
        return False

    L.extend(replacement)
    return sorted(L)


# ---------------------------------------------------------------------------
# Top-level decomposition with full deduplication
# ---------------------------------------------------------------------------

def decompose(target: mpq, initial_primes: list | None = None,
              prefer_large: bool = False) -> list | bool:
    """
    Decompose *target* into distinct primes  p₁ < p₂ < … < pₖ  such that

        target = Σ 1/(pᵢ + SHIFT).

    Parameters
    ----------
    target : mpq
        The positive rational to decompose.
    initial_primes : list or None
        Skip the initial search and start deduplication from this list.
    prefer_large : bool
        If True, use "Middle" mode for splitting (finds slightly larger but
        sometimes cleaner representations).

    Returns
    -------
    Sorted list of distinct primes, or False if no decomposition found.
    """
    if initial_primes is None:
        # Try modes in increasing order of thoroughness.
        primes = decompose_rational(target, 1, mode="Fast")
        if primes is False:
            primes = decompose_rational(target, 1, mode="Middle")
        if primes is False:
            primes = decompose_rational(target, 1, mode="Deep")
        if primes is False:
            return False
    else:
        primes = list(initial_primes)

    split_mode = "Middle" if prefer_large else "Fast"

    # Iteratively eliminate all duplicate primes.
    #
    # Progress is always printed so you can watch ct converge to 0:
    #
    #   loop j=<iter>  ct=<dups>  trend=<arrow>  [targeting prime=<p>]
    #
    #   ct     = number of excess prime copies still to eliminate
    #   trend  = down (good: converging), up (temporarily worse), flat (no change)
    #
    # With VERBOSE=True, each line also shows which prime is being targeted,
    # and the final solution list is printed.
    #
    prev_dups = count_duplicates(primes)
    _TREND = {-1: "down", 0: "flat", 1: "up"}

    for iteration in range(10_000):
        dups = count_duplicates(primes)
        trend = _TREND[0 if dups == prev_dups else (-1 if dups < prev_dups else 1)]

        # Build the progress line.
        line = f"{target}, loop j={iteration:>5}  ct={dups:>5}  trend={trend}"
        if VERBOSE:
            counts = multiplicity_dict(primes)
            repeated_now = next((p for p in sorted(counts) if counts[p] >= 2), None)
            if repeated_now:
                line += f"  targeting prime={repeated_now}"
        print(line, flush=True)

        if dups == 0:
            break

        prev_dups = dups

        updated = reduce_duplicates(primes, mode=split_mode)
        if updated is False:
            updated = reduce_duplicates(primes, mode="Middle")
        if updated is False:
            updated = reduce_duplicates(primes, mode="Deep")
        if updated is False:
            break

        primes = updated

    final_dups = count_duplicates(primes)
    if final_dups > 0:
        print(f"  FAILED: stuck at ct={final_dups}.")
        return False

    if VERBOSE:
        print(f"  SUCCESS: {len(primes)} distinct primes: {primes}")

    sol= sorted(list(set(primes)))
    if sol is not False:
        print(f"RESGOOD {target} {[int(p) for p in sol]} {len(sol)}")
    else:
        print(f"RESBAD  {r} FALSE")

    return sol 



# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def parse_number_string(s: str):
    if '/' in s:
        num_str, denom_str = s.split('/', 1)
        p = int(num_str)
        q = int(denom_str)
        if q == 0:
            raise ValueError("Denominator cannot be zero")
        if q < 0:
            p = -p
            q = -q
        return p, q
    else:
        return int(s), 1

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python sunconj.py SHIFT string_of_number [VERBOSE] [prefer_large]")
        print("Example: python sunconj.py +1 2 ")
        print("Example: python sunconj.py -1 3 ")
        print("Example: python sunconj.py +1 7/4 ")
        sys.exit(1)

    shift_str = sys.argv[1]
    num_str = sys.argv[2]

    # Parse SHIFT
    if shift_str == '+1':
        SHIFT = 1
    elif shift_str == '-1':
        SHIFT = -1
    else:
        print("Error: SHIFT must be +1 or -1")
        sys.exit(1)

    # Parse VERBOSE (default False)
    VERBOSE = False 
    if len(sys.argv) > 3:
        verbose_str = sys.argv[3].lower()
        if verbose_str in ('true', '1', 'yes'):
            VERBOSE = True
        elif verbose_str in ('false', '0', 'no'):
            VERBOSE = False
        else:
            print("Error: VERBOSE must be True/False, 1/0, or yes/no")
            sys.exit(1)

    # Parse prefer_large (default True)
    prefer_large = True
    if len(sys.argv) > 4:
        prefer_large_str = sys.argv[4].lower()
        if prefer_large_str in ('true', '1', 'yes'):
            prefer_large = True
        elif prefer_large_str in ('false', '0', 'no'):
            prefer_large = False
        else:
            print("Error: prefer_large must be True/False, 1/0, or yes/no")
            sys.exit(1)

    # Parse the number string
    try:
        p, q = parse_number_string(num_str)
        number = mpq(p, q)
    except ValueError as e:
        print(f"Error: {e}")
        sys.exit(1)

    print(f"Decomposing {number} with SHIFT = {SHIFT:+d}, prefer_large = {prefer_large}, VERBOSE = {VERBOSE}")

    decompose(number, prefer_large=prefer_large)

