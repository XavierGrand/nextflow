FROM quay.io/mbhall88/rasusa:0.6.0 AS quay_source

FROM alpine:3.13
COPY --from=quay_source /bin/rasusa /bin/
RUN apk add --update --no-cache bash procps

CMD [ "bash" ]
