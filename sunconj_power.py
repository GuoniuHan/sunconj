import os
import sys

def check_powers_of_two(MaxK, shift=1, verbose=True, prefer_large=True):
    """
    Check 1/2^k for k from 1 to MaxK.
    Calls sunconj.py for each power and prints results in real-time.
    """
    for k in range(1, MaxK + 1):
        print("RES k=", k)
        denominator = 2 ** k
        number_str = f"1/{denominator}"
        shift_str = f"{shift:+d}"
        verbose_str = str(verbose).lower()
        prefer_large_str = str(prefer_large).lower()

        # Call sunconj.py using os.system (real-time output)
        cmd = (
            f"python sunconj.py {shift_str} {number_str} {verbose_str} {prefer_large_str}"
        )
        print(f"\n--- Checking 1/2^{k} = 1/{denominator} ---")
        os.system(cmd)  # Output from sunconj.py appears immediately

    print(f"\nChecked all powers of two up to k = {MaxK}.")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python check_powers_of_two.py MaxK [SHIFT] [VERBOSE] [prefer_large]")
        print("Example: python check_powers_of_two.py 10 +1 True True")
        sys.exit(1)

    MaxK = int(sys.argv[1])
    shift = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    verbose = sys.argv[3].lower() in ('true', '1', 'yes') if len(sys.argv) > 3 else True
    prefer_large = sys.argv[4].lower() in ('true', '1', 'yes') if len(sys.argv) > 4 else True

    check_powers_of_two(MaxK, shift, verbose, prefer_large)
