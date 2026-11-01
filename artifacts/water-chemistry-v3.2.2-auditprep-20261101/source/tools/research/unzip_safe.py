"""Extract members of a zip archive into a target directory without trusting member paths.

Usage: python3 -I unzip_safe.py <zip_path> <out_dir>
Only the basename of each member is used; directories and names containing
path separators or leading dots are rejected. Existing files are not overwritten.
"""
import os
import sys
import zipfile


def main():
    zpath, out = sys.argv[1], sys.argv[2]
    with zipfile.ZipFile(zpath) as z:
        for info in z.infolist():
            if info.is_dir():
                continue
            name = os.path.basename(info.filename.replace("\\", "/"))
            if not name or name.startswith(".") or "/" in name:
                raise SystemExit(f"rejected member name: {info.filename!r}")
            dest = os.path.join(out, name)
            if os.path.exists(dest):
                raise SystemExit(f"refusing to overwrite {dest}")
            if info.file_size > 50_000_000:
                raise SystemExit(f"member too large: {info.filename!r}")
            with z.open(info) as src, open(dest, "wb") as dst:
                dst.write(src.read())
            print(name, info.file_size)


if __name__ == "__main__":
    main()
