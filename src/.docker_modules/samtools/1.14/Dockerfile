FROM quay.io/biocontainers/samtools:1.14--hb421002_0 AS build

# ldd /usr/local/bin/samtools
#         linux-vdso.so.1 (0x00007ffec91ca000)
#         libhts.so.3 => /usr/local/bin/../lib/libhts.so.3 (0x00007fa2fc31f000)
#         libtinfow.so.6 => /usr/local/bin/../lib/libtinfow.so.6 (0x00007fa2fc2e0000)
#         libncursesw.so.6 => /usr/local/bin/../lib/libncursesw.so.6 (0x00007fa2fc2a3000)
#         libm.so.6 => /lib/x86_64-linux-gnu/libm.so.6 (0x00007fa2fc11f000)
#         libz.so.1 => /usr/local/bin/../lib/libz.so.1 (0x00007fa2fc105000)
#         libpthread.so.0 => /lib/x86_64-linux-gnu/libpthread.so.0 (0x00007fa2fc0e4000)
#         libc.so.6 => /lib/x86_64-linux-gnu/libc.so.6 (0x00007fa2fbf21000)
#         libdeflate.so.0 => /usr/local/bin/../lib/./libdeflate.so.0 (0x00007fa2fbf11000)
#         liblzma.so.5 => /usr/local/bin/../lib/./liblzma.so.5 (0x00007fa2fbee8000)
#         libbz2.so.1.0 => /usr/local/bin/../lib/./libbz2.so.1.0 (0x00007fa2fbed4000)
#         libdl.so.2 => /lib/x86_64-linux-gnu/libdl.so.2 (0x00007fa2fbecf000)
#         /lib64/ld-linux-x86-64.so.2 (0x00007fa2fc4af000)
#         librt.so.1 => /lib/x86_64-linux-gnu/librt.so.1 (0x00007fa2fbec3000)

ENV SAMTOOLS_VERSION=1.14
ENV PACKAGES bash
FROM alpine:3.8 AS env
RUN apk update && \
    apk add ${PACKAGES}
COPY --from=build \
    /usr/local/bin/samtools /usr/local/bin/
COPY --from=build \
    /usr/local/bin/../lib/libhts.so.3 \
    /usr/local/lib/libtinfow.so.6 \
    /usr/local/lib/libncursesw.so.6 \
    /usr/local/lib/libz.so.1 \
    /usr/local/lib/libdeflate.so.0 \
    /usr/local/lib/liblzma.so.5 \
    /usr/local/lib/libbz2.so.1.0 /usr/local/lib/
COPY --from=build \
    /lib/x86_64-linux-gnu/libm.so.6 \
    /lib/x86_64-linux-gnu/libpthread.so.0 \
    /lib/x86_64-linux-gnu/libc.so.6 \
    /lib/x86_64-linux-gnu/libdl.so.2 \
    /lib/x86_64-linux-gnu/librt.so.1 /lib/x86_64-linux-gnu/
COPY --from=build \
    /lib64/ld-linux-x86-64.so.2 /lib64/