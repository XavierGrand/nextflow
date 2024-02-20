FROM lbmc/samtools:1.14 as build

ENV HISAT2_VERSION=2.2.1

FROM quay.io/biocontainers/hisat2:2.2.1--h1b792b2_3
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