# DMRdetector image. Build from the repository root:
#   docker build -f env/Dockerfile -t dmrdetector:3.2.0 .
# Release builds require the compiled kernels (the build stops when the kernel
# cache cannot be made, and the smoke test asserts that it loads):
#   docker build -f env/Dockerfile --build-arg DMRDETECTOR_REQUIRE_CACHE=1 -t dmrdetector:3.2.0 .
# env/dmrdetector.def is the Apptainer recipe of the same image. The two files are written by hand; keep them in step
# (tests/unit/test_container_recipes.py fails when the base image, the install command, the ENV lines, the libdeflate
# pin or the version label differ).

# Stage 1: install the hash-locked packages into a venv, put the package source (with its data files) at
# /opt/dmrdetector/src, make it importable through a .pth file, add libdeflate, prune, strip, byte-compile, then
# compile the numerical kernels once (env/warm_cache.sh) into /opt/dmrdetector/numba-cache for a generic x86-64 CPU.
FROM python:3.12.15-slim-trixie@sha256:02108f5d322dd89f1c9e552442c25acb0543dfdbc455693a5599624f20d9155d AS build
ENV PIP_NO_CACHE_DIR=1 PIP_DISABLE_PIP_VERSION_CHECK=1 PYTHONDONTWRITEBYTECODE=1
ENV LIBDEFLATE_DEB_VERSION=1.23-2 LIBDEFLATE_DEB_SHA256=5fe37fb79d53f5c91da6ba3ac4c1a8453e1f2b1ff0e0ea45189c5338bbfdbf82
COPY env/requirements.lock /build/requirements.lock
COPY env/launcher.py /opt/dmrdetector/launcher.py
COPY src /opt/dmrdetector/src
COPY env/warm_cache.sh /build/warm_cache.sh
COPY tests/smoke/generate_data.py tests/smoke/run_smoke.py env/warm_big.py /build/smoke/
RUN apt-get update \
 && apt-get install -y --no-install-recommends binutils \
 && rm -rf /var/lib/apt/lists/*
RUN python -m venv --without-pip /opt/venv \
 && pip --python /opt/venv/bin/python install --require-hashes --only-binary=:all: --no-compile -r /build/requirements.lock \
 && find /opt/dmrdetector/src -name __pycache__ -prune -exec rm -rf {} + \
 && chmod -R a+rX /opt/dmrdetector \
 && echo /opt/dmrdetector/src > /opt/venv/lib/python3.12/site-packages/dmrdetector.pth \
 && install -m 755 /opt/dmrdetector/launcher.py /opt/venv/bin/dmrdetector
RUN mkdir /build/deb && cd /build/deb \
 && apt-get update \
 && apt-get download libdeflate0=$LIBDEFLATE_DEB_VERSION \
 && echo "$LIBDEFLATE_DEB_SHA256  libdeflate0_${LIBDEFLATE_DEB_VERSION}_amd64.deb" | sha256sum -c - \
 && dpkg -x libdeflate0_${LIBDEFLATE_DEB_VERSION}_amd64.deb /build/deb/root \
 && cp -L /build/deb/root/usr/lib/x86_64-linux-gnu/libdeflate.so.0 /opt/venv/lib/libdeflate.so.0 \
 && echo "$LIBDEFLATE_DEB_VERSION" > /opt/venv/lib/libdeflate.version \
 && rm -rf /build/deb /var/lib/apt/lists/*
RUN find /opt/venv -depth -type d \( -name tests -o -name test -o -name __pycache__ -o -name include \) -exec rm -rf {} + \
 && find /opt/venv -type f \( -name '*.pyc' -o -name '*.pyi' -o -name '*.pxd' -o -name '*.pyx' -o -name '*.a' \) -delete \
 && rm -f /opt/venv/bin/f2py /opt/venv/bin/numpy-config \
 && find /opt/venv -type f \( -name '*.so' -o -name '*.so.*' \) -exec strip --strip-unneeded {} + \
 && python -m compileall -q --invalidation-mode unchecked-hash /opt/venv /opt/dmrdetector/src
ENV NUMBA_CPU_NAME=generic
ARG DMRDETECTOR_REQUIRE_CACHE=0
RUN bash /build/warm_cache.sh /opt/dmrdetector/numba-cache /build/smoke

# Stage 2: the same base image (the venv interpreter is a symlink into /usr/local), no apt, no compilers, no pip.
FROM python:3.12.15-slim-trixie@sha256:02108f5d322dd89f1c9e552442c25acb0543dfdbc455693a5599624f20d9155d AS final
COPY --from=build /opt/venv /opt/venv
COPY --from=build /opt/dmrdetector/src /opt/dmrdetector/src
COPY --from=build /opt/dmrdetector/numba-cache /opt/dmrdetector/numba-cache
COPY env/requirements.lock /opt/dmrdetector/requirements.lock
COPY env/smoke_test.py /opt/dmrdetector/smoke_test.py
RUN rm -rf /usr/local/lib/python3.12/site-packages/* /usr/local/lib/python3.12/ensurepip \
    /usr/local/lib/python3.12/idlelib /usr/local/lib/python3.12/tkinter /usr/local/lib/python3.12/turtledemo \
    /usr/local/lib/python3.12/turtle.py /usr/local/lib/python3.12/lib2to3 /usr/local/include \
    /usr/local/bin/pip* /usr/local/bin/idle* /usr/local/bin/pydoc* /usr/local/share/man \
    /var/lib/apt/lists/* /var/cache/apt /usr/share/doc/* /usr/share/man/* /usr/share/info/*
RUN python -m compileall -q --invalidation-mode unchecked-hash /usr/local/lib/python3.12
ENV PATH=/opt/venv/bin:$PATH
ENV PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1 PYTHONPATH=
ENV OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 NUMBA_CPU_NAME=generic
# glibc picks the FMA or the plain variant of exp, log, pow, sin and cos from the CPU; scipy.fft and
# scipy.stats call them, so the image fixes the plain variant (doc/dependencies.md)
ENV GLIBC_TUNABLES=glibc.cpu.hwcaps=-FMA,-FMA4
LABEL org.opencontainers.image.title="DMRdetector" org.opencontainers.image.version="3.2.0" org.opencontainers.image.authors="Jen-Hsiang Ou" org.opencontainers.image.licenses="MIT"
# Smoke test (the Apptainer recipe runs the same command as %test)
ARG DMRDETECTOR_REQUIRE_CACHE=0
RUN dmrdetector --version && python /opt/dmrdetector/smoke_test.py
ENTRYPOINT ["dmrdetector"]
