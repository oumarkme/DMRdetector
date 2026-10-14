"""The `jit` decorator: numba nopython, no fastmath, no parallel reductions.

The compiled code has the same float64 arithmetic as the Python source, so
results do not depend on the thread count. DMRDETECTOR_NO_JIT=1 runs the plain
Python source. Compiled code is cached in NUMBA_CACHE_DIR; when that is unset,
`set_default_cache_dir` points it to a private folder of the user in the
temporary directory. The command line calls it at start, before the temporary
directory moves into the work folder of a run; `jit` calls it for every other
caller. The image holds a cache of the kernels, compiled at build time;
`set_default_cache_dir` copies it into the folder in use on first use, so a
fresh container does not compile.

Numba unpickles the files of its cache, so whoever can write into the cache
folder can run code as the user. The per-user folder is therefore used only
when it is a real folder of the user with no group or world write bit, in a
temporary directory that no other user can rename it out of; otherwise the
program makes a new private folder for the process and removes it at exit.
"""

from __future__ import annotations

import atexit
from collections.abc import Callable
import functools
import os
import shutil
import stat
import tempfile


# The cache baked into the image (see env/warm_cache.sh). Numba cannot use it
# in place: it tests that a cache folder is writable, and the image is not.
SEED_DIR = "/opt/dmrdetector/numba-cache"
SEED_ENV = "DMRDETECTOR_NUMBA_SEED"  # another seed folder; empty: no seeding
_SEED_ID = "SEED"  # file of the seed folder that names its build
_SEEDED = ".seeded"  # file of the cache folder that names the seed copied
_NO_SEED = "noseed"  # identity of a cache folder without a seed


def seed_identity(seed: str | None = None) -> str:
    """Returns the short identity of a seed: a hash of its SEED file.

    The SEED file is written when the image is built and names the program
    version, the CPU the kernels were compiled for and a hash of the cache
    files (env/warm_cache.sh), so two builds never share an identity and no
    cache file is read to find it.

    Args:
        seed: Seed folder; default `DMRDETECTOR_NUMBA_SEED`, else the folder of
            the image.

    Returns:
        12 hexadecimal digits, or "noseed" when there is no seed.
    """
    import hashlib  # not needed by --help and --version
    if seed is None:
        seed = os.environ.get(SEED_ENV, SEED_DIR)
    if not seed:
        return _NO_SEED
    try:
        with open(os.path.join(seed, _SEED_ID), "rb") as fh:
            data = fh.read(4096)
    except OSError:
        return _NO_SEED
    return hashlib.sha256(data).hexdigest()[:12]


def cache_folder_name(uid: int) -> str:
    """Returns the name of the per-user cache folder for the current seed.

    Versions and builds of the image name different folders, so a folder never
    holds the entries of two seeds.
    """
    return f"dmrdetector-numba-{uid}-{seed_identity()}"


def _copy_private(src: str, dst: str) -> None:
    """Copies `src` to `dst` through a unique temporary file, mode 0600.

    The temporary name comes from `tempfile.mkstemp` (exclusive creation), so
    processes, threads and containers that share the folder never write the
    same file. `os.replace` makes the file appear complete or not at all.
    """
    fd, tmp = tempfile.mkstemp(dir=os.path.dirname(dst),
                               prefix=os.path.basename(dst) + ".tmp")
    try:
        with os.fdopen(fd, "wb") as out, open(src, "rb") as inp:
            shutil.copyfileobj(inp, out)
        os.replace(tmp, dst)
    except BaseException:
        try:
            os.unlink(tmp)
        except OSError:
            pass
        raise


def seed_cache(path: str) -> bool:
    """Copies the compiled kernels of the image into the cache folder `path`.

    The copy happens once per folder and seed: `path/.seeded` holds the content
    of the seed's SEED file, and nothing is copied when it is the same. Files
    are copied under a unique temporary name and moved into place, so
    processes that start together never read a half-written file; the data
    files go first and the index files last, and `.seeded` after all of them.
    Links in the seed are ignored. Any error leaves the folder as it was
    (numba then compiles), and the function returns False.

    Args:
        path: Cache folder, writable by the current user.

    Returns:
        True when files were copied, False when there was nothing to do or
        the copy failed.
    """
    seed = os.environ.get(SEED_ENV, SEED_DIR)
    if not seed:
        return False
    try:
        with open(os.path.join(seed, _SEED_ID), encoding="utf-8") as fh:
            name = fh.read()
        try:
            with open(os.path.join(path, _SEEDED), encoding="utf-8") as fh:
                if fh.read() == name:
                    return False
        except FileNotFoundError:
            pass
        files = []
        for top, _, names in os.walk(seed, followlinks=False):
            for n in names:
                src = os.path.join(top, n)
                if n.endswith((".nbc", ".nbi")) and not os.path.islink(src):
                    files.append(src)
        files.sort(key=lambda f: (f.endswith(".nbi"), f))
        for src in files:
            dst = os.path.join(path, os.path.relpath(src, seed))
            os.makedirs(os.path.dirname(dst), mode=0o700, exist_ok=True)
            _copy_private(src, dst)
        fd, tmp = tempfile.mkstemp(dir=path, prefix=f"{_SEEDED}.tmp")
        try:
            with os.fdopen(fd, "w", encoding="utf-8") as fh:
                fh.write(name)
            os.replace(tmp, os.path.join(path, _SEEDED))
        except BaseException:
            try:
                os.unlink(tmp)
            except OSError:
                pass
            raise
        return True
    except OSError:
        return False


@functools.lru_cache(maxsize=None)
def _seed_once(path: str) -> None:
    """Seeds `path` on the first call with it in the process."""
    seed_cache(path)


def _parent_is_safe(parent: str) -> bool:
    """Returns True unless others could rename a folder out of `parent`.

    That is the case for a folder that group or world can write to and that has
    no sticky bit (like a shared scratch folder with mode 0777). A folder made
    there can be replaced by one with other cache files, and numba unpickles
    them.
    """
    try:
        mode = os.stat(parent).st_mode
    except OSError:
        return False
    return not (mode & 0o022 and not mode & stat.S_ISVTX)


def _fallback_folder() -> str | None:
    """Makes a private folder for this process with mkdtemp, or returns None.

    It is made in the temporary directory, or, when others could rename a
    folder out of that, in the user cache directory ($XDG_CACHE_HOME or
    ~/.cache). Without a safe parent there is no folder, and numba uses its
    own default.
    """
    cache = os.environ.get("XDG_CACHE_HOME") or os.path.join(
        os.path.expanduser("~"), ".cache")
    for parent in (tempfile.gettempdir(), cache):
        if not _parent_is_safe(parent):
            continue
        try:
            return tempfile.mkdtemp(prefix="dmrdetector-numba-", dir=parent)
        except OSError:
            continue
    return None


def _is_private_folder(path: str, uid: int) -> bool:
    """Returns True when `path` is a folder that only `uid` can write to.

    The folder is made (mode 0700) when it is missing. It is refused, and
    never repaired, when it is not a folder, is a link, belongs to another
    user, or has a group or world write bit: someone else may have put cache
    files in it, and numba unpickles them. The same holds when the parent can
    be written by others and is not sticky, because anyone could then rename
    the folder away and put another in its place. A folder that others may
    only read or enter is closed to 0700.
    """
    try:
        if not _parent_is_safe(os.path.dirname(path)):
            return False
        os.makedirs(path, mode=0o700, exist_ok=True)
        st = os.lstat(path)
        if not stat.S_ISDIR(st.st_mode) or st.st_uid != uid:
            return False
        if st.st_mode & 0o022:
            return False
        if st.st_mode & 0o077:
            os.chmod(path, 0o700)
        return True
    except OSError:
        return False


def _remove_folder(path: str, pid: int) -> None:
    """Removes `path` at exit, only in the process that made it."""
    if os.getpid() == pid:
        shutil.rmtree(path, ignore_errors=True)


def set_default_cache_dir() -> str | None:
    """Sets NUMBA_CACHE_DIR to a private per-user folder when it is unset.

    The folder is <temporary directory>/dmrdetector-numba-<euid>-<seed
    identity> with mode 0700 (`cache_folder_name`). It is used only when it is
    a real folder of the effective user without group or world write bit, in a
    temporary directory that others cannot rename it out of
    (`_is_private_folder`). Otherwise the process makes a new folder with
    mkdtemp (mode 0700) and removes it, with the kernels copied into it, when
    it exits; the per-user folder is never removed. The new folder is made in
    the temporary directory, or in the user cache directory when others could
    rename a folder out of the temporary directory; when neither is safe,
    NUMBA_CACHE_DIR stays unset. A variable that the caller
    set is never changed. The folder in use, the caller's or ours, gets the
    kernels of the image once (`seed_cache`).

    Returns:
        The folder set, or None when nothing was changed.
    """
    if os.environ.get("DMRDETECTOR_NO_JIT") == "1":
        return None
    if os.environ.get("NUMBA_CACHE_DIR"):
        _seed_once(os.environ["NUMBA_CACHE_DIR"])
        return None
    uid = os.geteuid()
    path = os.path.join(tempfile.gettempdir(), cache_folder_name(uid))
    if not _is_private_folder(path, uid):
        folder = _fallback_folder()
        if folder is None:
            return None  # numba falls back to its own default
        path = folder
        atexit.register(_remove_folder, path, os.getpid())
    os.environ["NUMBA_CACHE_DIR"] = path
    _seed_once(path)
    return path


def jit(*args: object, **kwargs: object) -> Callable:
    """Compiles a function with numba.njit(cache=True, nogil=True).

    Used as `@jit` or `@jit(...)`. With DMRDETECTOR_NO_JIT=1 the function is
    returned unchanged.

    Args:
        *args: The function when used as `@jit`.
        **kwargs: Options for numba.njit; fastmath is always off.

    Returns:
        The compiled function, or a decorator when called with options.
    """
    if os.environ.get("DMRDETECTOR_NO_JIT") == "1":

        def passthrough(fn):
            """Returns `fn` unchanged."""
            return fn

        if args and callable(args[0]):
            return passthrough(args[0])
        return passthrough
    # a caller that does not go through the command line (the image smoke test,
    # pytest, an import) gets the private cache folder too
    set_default_cache_dir()
    import numba

    # numba read the variable when it was imported; a caller that imported it
    # before this point (the image smoke test) needs the new value
    if numba.config.CACHE_DIR != os.environ.get("NUMBA_CACHE_DIR", ""):
        numba.config.reload_config()
    kwargs.setdefault("cache", True)
    kwargs.setdefault("nogil", True)
    kwargs["fastmath"] = False
    if args and callable(args[0]):
        return numba.njit(**kwargs)(args[0])
    return numba.njit(**kwargs)
