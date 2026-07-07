import os
import sys
from gmpy2 import gcd

def verify_all_fractions(bound: int, shift: int = 1, prefer_large: bool = False, verbose: bool = False):

    for n in range(2, bound + 1):
        for k in range(1, n):
            if gcd(n, k) == 1:
                #p = n + k  # TODO
                p = k  # TODO
                q = n
                number_str = f"{p}/{q}"
                shift_str = f"{shift:+d}"
                prefer_large_str = str(prefer_large).lower()
                verbose_str = str(verbose).lower()

                # Call sunconj.py using os.system (real-time output)
                cmd = (
                    f"python sunconj.py {shift_str} {number_str} "
                    f"{prefer_large_str} {verbose_str}"
                )
                os.system(cmd)  # Output appears in real-time

if __name__ == "__main__":
    bound = int(sys.argv[1]) if len(sys.argv) > 1 else 200
    shift = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    prefer_large = sys.argv[3].lower() in ('true', '1', 'yes') if len(sys.argv) > 3 else False
    verbose = sys.argv[4].lower() in ('true', '1', 'yes') if len(sys.argv) > 4 else False

    verify_all_fractions(bound, shift, prefer_large, verbose)
