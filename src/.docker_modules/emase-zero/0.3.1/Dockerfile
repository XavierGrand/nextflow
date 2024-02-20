FROM debian:buster as build
MAINTAINER Laurent Modolo

ENV EMASEZERO_VERSION=0.3.1
ENV PACKAGES build-essential \
    ca-certificates \
    zlib1g-dev \
    git

RUN apt-get update && \
    apt-get install -y --no-install-recommends ${PACKAGES} && \
    apt-get clean

RUN git clone https://github.com/churchill-lab/emase-zero.git \
    && cd emase-zero/src/ \
    && make \
    && mv emase-zero /usr/local/bin/ \
    && ldd /usr/local/bin/emase-zero
# linux-vdso.so.1 (0x00007ffe8e35f000)
# libz.so.1 => /lib/x86_64-linux-gnu/libz.so.1 (0x00007fbd358b5000)
# libstdc++.so.6 => /usr/lib/x86_64-linux-gnu/libstdc++.so.6 (0x00007fbd35731000)
# libm.so.6 => /lib/x86_64-linux-gnu/libm.so.6 (0x00007fbd355ae000)
# libgcc_s.so.1 => /lib/x86_64-linux-gnu/libgcc_s.so.1 (0x00007fbd35594000)
# libc.so.6 => /lib/x86_64-linux-gnu/libc.so.6 (0x00007fbd353d3000)
# /lib64/ld-linux-x86-64.so.2 (0x00007fbd35af7000)

FROM debian:buster-slim

COPY --from=build /usr/local/bin/emase-zero /usr/local/bin/
COPY --from=build /lib/x86_64-linux-gnu/libz.so.1 /lib/x86_64-linux-gnu/
COPY --from=build /lib/x86_64-linux-gnu/libm.so.6 /lib/x86_64-linux-gnu/
COPY --from=build /lib/x86_64-linux-gnu/libgcc_s.so.1 /lib/x86_64-linux-gnu/

RUN apt-get update && \
    apt-get install -y procps bash
