FROM quay.io/biocontainers/samtools:0.1.18--hfb9b9cc_10 AS samtools 
# /usr/local/bin/samtools
# / # ldd /usr/local/bin/samtools
#         /lib64/ld-linux-x86-64.so.2 (0x7efddcdcc000)
#         libncurses.so.6 => /usr/local/bin/../lib/libncurses.so.6 (0x7efddcfad000)
#         libtinfo.so.6 => /usr/local/bin/../lib/libtinfo.so.6 (0x7efddcf6f000)
#         libm.so.6 => /lib64/ld-linux-x86-64.so.2 (0x7efddcdcc000)
#         libz.so.1 => /usr/local/bin/../lib/libz.so.1 (0x7efddcf55000)
#         libc.so.6 => /lib64/ld-linux-x86-64.so.2 (0x7efddcdcc000)

FROM quay.io/biocontainers/minimap2:2.24--h7132678_1
MAINTAINER adapted by Xavier Grand from Laurent Modolo template
ENV MINIMAP2_VERSION=2.24

COPY --from=samtools /usr/local/bin/samtools /usr/local/bin/
COPY --from=samtools /usr/local//lib/libncurses.so.6 /usr/local/lib/
COPY --from=samtools /usr/local//lib/libtinfo.so.6 /usr/local/lib/

# /usr/local/bin/minimap2
# / # ldd /usr/local/bin/minimap2
#         /lib64/ld-linux-x86-64.so.2 (0x7fe14f5a8000)
#         libm.so.6 => /lib64/ld-linux-x86-64.so.2 (0x7fe14f5a8000)
#         libz.so.1 => /usr/local/bin/../lib/libz.so.1 (0x7fe14f7c4000)
#         libpthread.so.0 => /lib64/ld-linux-x86-64.so.2 (0x7fe14f5a8000)
#         libc.so.6 => /lib64/ld-linux-x86-64.so.2 (0x7fe14f5a8000)