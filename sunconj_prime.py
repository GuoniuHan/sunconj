import os
import sys
import gmpy2  # For prime checking

def check_prime_reciprocals(M, shift=1, verbose=True, prefer_large=True):
    """
    Check 1/p for all prime numbers p ≤ M.
    Calls sunconj.py for each prime and prints results in real-time.
    """
    failed = []
    for p in range(2, M + 1):
        if gmpy2.is_prime(p):  # Check if p is prime
            number_str = f"1/{p}"
            shift_str = f"{shift:+d}"
            verbose_str = str(verbose).lower()
            prefer_large_str = str(prefer_large).lower()

            # Call sunconj.py using os.system (real-time output)
            cmd = (
                f"python sunconj.py {shift_str} {number_str} {verbose_str} {prefer_large_str}"
            )
            print(f"\n--- Checking 1/{p} ---")  # Optional: Add a separator for clarity
            os.system(cmd)  # Output from sunconj.py appears immediately

    print(f"\nChecked all primes up to {M}.")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python check_prime_reciprocals.py M [SHIFT] [VERBOSE] [prefer_large]")
        print("Example: python check_prime_reciprocals.py 100 +1 True True")
        sys.exit(1)

    M = int(sys.argv[1])
    shift = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    verbose = sys.argv[3].lower() in ('true', '1', 'yes') if len(sys.argv) > 3 else True
    prefer_large = sys.argv[4].lower() in ('true', '1', 'yes') if len(sys.argv) > 4 else True

    check_prime_reciprocals(M, shift, verbose, prefer_large)

