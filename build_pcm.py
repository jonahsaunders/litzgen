"""Build the KiCad Plugin and Content Manager (PCM) package.

    python build_pcm.py            -> dist/litzgen-pcm.zip  (+ dist/litzgen-pcm-<version>.zip)

The version comes from litzgen/__init__.py and is written into the packaged
metadata.json, so a release only needs a version bump and a tag.
The unversioned file name is what the README's "latest release" link points at.
"""
import hashlib
import json
import os
import re
import shutil
import zipfile

ROOT = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.join(ROOT, "litzgen")
DIST = os.path.join(ROOT, "dist")


def version() -> str:
    txt = open(os.path.join(PKG, "__init__.py"), encoding="utf-8").read()
    return re.search(r'__version__\s*=\s*"([^"]+)"', txt).group(1)


def build() -> str:
    ver = version()
    meta = json.load(open(os.path.join(ROOT, "metadata.json"), encoding="utf-8"))
    meta["versions"][0]["version"] = ver
    os.makedirs(DIST, exist_ok=True)
    out = os.path.join(DIST, "litzgen-pcm.zip")
    install_size = 0
    with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("metadata.json", json.dumps(meta, indent=2) + "\n")
        z.write(os.path.join(ROOT, "resources", "icon.png"), "resources/icon.png")
        for d, dirs, files in os.walk(PKG):
            dirs[:] = [x for x in dirs if x != "__pycache__"]
            for f in sorted(files):
                if f.endswith((".py", ".png")):
                    full = os.path.join(d, f)
                    install_size += os.path.getsize(full)
                    z.write(full, os.path.join("plugins", os.path.relpath(full, PKG)).replace(os.sep, "/"))
    shutil.copyfile(out, os.path.join(DIST, "litzgen-pcm-%s.zip" % ver))
    sha = hashlib.sha256(open(out, "rb").read()).hexdigest()
    print("built %s (version %s)" % (out, ver))
    print("download_sha256: %s" % sha)
    print("download_size:   %d" % os.path.getsize(out))
    print("install_size:    %d" % install_size)
    return out


if __name__ == "__main__":
    build()
