FROM kbase/sdkbase2:python
MAINTAINER KBase Developer
# -----------------------------------------
# In this section, you can install any system dependencies required
# to run your App.  For instance, you could place an apt-get update or
# install line here, a git checkout to download code, or run any other
# installation scripts.

RUN apt-get update && apt-get install -y bowtie samtools bcftools bedtools emboss python3-pip python3-numpy wget zlib1g-dev emacs bc gawk
# need newer bcftools than in repo for this old image:
RUN wget http://ftp.us.debian.org/debian/pool/main/h/htslib/libhts2_1.9-11_amd64.deb && wget http://ftp.us.debian.org/debian/pool/main/b/bcftools/bcftools_1.9-1_amd64.deb && wget http://ftp.us.debian.org/debian/pool/main/libd/libdeflate/libdeflate0_1.3-3_amd64.deb && dpkg --install *.deb
RUN wget https://github.com/bedops/bedops/releases/download/v2.4.35/bedops_linux_x86_64-v2.4.35.tar.bz2 && tar jxf bedops_linux_x86_64-v2.4.35.tar.bz2 && cp bin/* /usr/local/bin/
RUN pip3 install --upgrade pip
RUN pip3 install openopt
RUN pip3 install biopython
RUN pip3 install pandas
RUN pip3 install numpy
RUN pip3 install scipy
RUN pip3 install FuncDesigner
RUN pip3 install DerApproximator
RUN pip3 install openopt
RUN git clone https://github.com/lh3/bwa.git && cd bwa && make && cp bwa /usr/local/bin/
RUN git clone -b py3 https://github.com/caozhichongchong/meta_decoder.git

# -----------------------------------------

COPY ./ /kb/module
RUN mkdir -p /kb/module/work
RUN chmod -R a+rw /kb/module

WORKDIR /kb/module

RUN make all

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]