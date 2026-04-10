import re
from fractions import Fraction
from pathlib import Path

# Matches things like '<1, 1>' (allowing spaces and optional +/-)
TUPLE_RE = re.compile(r"<\s*([-+]?\d+)\s*,\s*([-+]?\d+)\s*>")

def main():
    folder = Path("data_generic")
    total = Fraction(0, 1)

    if not folder.is_dir():
        raise FileNotFoundError(f"Folder not found: {folder}")

    for infile in sorted(folder.glob("*.txt")):
        with open(infile, "r", encoding="utf-8") as f:
            for lineno, line in enumerate(f, start=1):
                m = TUPLE_RE.search(line)
                if not m:
                    continue  # skip lines without a <x, y> tuple

                x = int(m.group(1))
                if x == 0:
                    raise ZeroDivisionError(
                        f"File {infile}, line {lineno}: first coordinate is 0, cannot add 1/x."
                    )

                total += Fraction(1, x)

    print(f"Exact sum: {total}")
    print(f"Float sum: {float(total)}")

if __name__ == "__main__":
    main()

