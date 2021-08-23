FROM kbase/sdkbase2:python
MAINTAINER KBase Developer
# -----------------------------------------
# In this section, you can install any system dependencies required
# to run your App.  For instance, you could place an apt-get update or
# install line here, a git checkout to download code, or run any other
# installation scripts.

WORKDIR /kb/module

# vcfstats needs python 3.6+, so need to move from stretch to buster
RUN echo "deb http://deb.debian.org/debian stable main contrib" > /etc/apt/sources.list
RUN apt-get update && apt-get install -y bowtie samtools bcftools bedtools emboss python3-pip wget zlib1g-dev emacs bc gawk vcftools r-base r-cran-ggplot2 curl libcurl4-openssl-dev libssl-dev gnuplot python-pip python-numpy r-cran-pheatmap poppler-utils python3-cyvcf2
RUN wget https://github.com/bedops/bedops/releases/download/v2.4.35/bedops_linux_x86_64-v2.4.35.tar.bz2 && tar jxf bedops_linux_x86_64-v2.4.35.tar.bz2 && cp bin/* /usr/local/bin/
RUN pip3 install openopt numpy scipy FuncDesigner DerApproximator Cython pandas biopython vcfstats
# still need some python2 code for strainfinder v1
RUN pip2 install openopt numpy scipy FuncDesigner DerApproximator Cython
RUN git clone https://github.com/lh3/bwa.git && cd bwa && git reset --hard 3ddd7b87d41f89a404d95f083fb37c369f70d783 && make && cp bwa /usr/local/bin/
RUN pip3 install --upgrade numpy
RUN git clone -b py3 https://github.com/caozhichongchong/meta_decoder.git && cd meta_decoder && git reset --hard 192bb32b13b0d0b73cdea872602171d8b5032649 && git clone https://bitbucket.org/yonatanf/strainfinder && cd strainfinder && git reset --hard 720087d84b4998a0ab7308880e48249c184dd0c2 && python2 setup_cython.py build_ext --inplace
ENV PATH $PATH:/kb/module/meta_decoder:/kb/module/meta_decoder/strainfinder

# -----------------------------------------

COPY ./ /kb/module
RUN mkdir -p /kb/module/work
RUN chmod -R a+rw /kb/module

WORKDIR /kb/module

RUN make all

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
