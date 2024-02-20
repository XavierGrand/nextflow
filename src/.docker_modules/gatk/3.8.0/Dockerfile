FROM broadinstitute/gatk3:3.8-0
MAINTAINER Laurent Modolo

ENV GATK_VERSION=3.8.0

RUN echo "#\!/bin/sh\njava -jar /usr/GenomeAnalysisTK.jar \$@" > /bin/gatk3 
RUN chmod +x /bin/gatk3
RUN sed -i 's/java/java -Xmx128g/g' /bin/gatk3
