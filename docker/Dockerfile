# Ubuntu base image 
FROM ubuntu:14.04
MAINTAINER Benjamin Schubert <schubert@infomratik.uni-tuebingen.de>

#installation of software
RUN apt-get update && apt-get install -y software-properties-common \
&& add-apt-repository ppa:george-edison55/cmake-3.x \
&& add-apt-repository ppa:ubuntu-toolchain-r/test \
&& apt-get update && apt-get install -y \
    gcc-4.9 \
    g++-4.9 \
    build-essential \
    coinor-cbc \
    git \
    mercurial \
    wget \
    pkg-config \
    python \
    python-pip \
    python-dev \
    python-matplotlib \
    libmysqlclient-dev \
    zlib1g-dev \
    tcsh \
    gawk \
    cmake \
    r-base \
    bowtie \
    libbz2-dev \
    libboost-dev \
&& update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-4.9 60 --slave /usr/bin/g++ g++ /usr/bin/g++-4.9 \
&& apt-get clean \
&& apt-get purge \
&& rm -rf /var/lib/apt/lists/*

#python stuff
RUN pip install --upgrade pip

#LKH solver here
RUN mkdir -p /usr/src/LKH \
    && wget --secure-protocol=TLSv1 http://www-bs3.informatik.uni-tuebingen.de/extern/schubert/LKH-2.0.7.tgz -O LKH.tgz \
    && tar -xzf LKH.tgz -C /usr/src/LKH \
    && make -C /usr/src/LKH/LKH-2.0.7 \
    && mv /usr/src/LKH/LKH-2.0.7/LKH /usr/local/bin/ \
    && rm LKH.tgz


#HLA Typing
#OptiType dependecies
RUN wget --secure-protocol=TLSv1 https://www.hdfgroup.org/ftp/HDF5/current/bin/linux-centos7-x86_64-gcc485/hdf5-1.8.17-linux-centos7-x86_64-gcc485-shared.tar.gz \
    && tar -xvf hdf5-1.8.17-linux-centos7-x86_64-gcc485-shared.tar.gz \
    && mv hdf5-1.8.17-linux-centos7-x86_64-gcc485-shared/bin/* /usr/local/bin/ \
    && mv hdf5-1.8.17-linux-centos7-x86_64-gcc485-shared/lib/* /usr/local/lib/ \
    && mv hdf5-1.8.17-linux-centos7-x86_64-gcc485-shared/include/* /usr/local/include/ \
    && mv hdf5-1.8.17-linux-centos7-x86_64-gcc485-shared/share/* /usr/local/share/ \
    && rm -rf hdf5-1.8.17-linux-centos7-x86_64-gcc485-shared/ \
    && rm -f hdf5-1.8.17-linux-centos7-x86_64-gcc485-shared.tar.gz

ENV LD_LIBRARY_PATH /usr/local/lib:$LD_LIBRARY_PATH
ENV HDF5_DIR /usr/local/

RUN pip install \
    numpy \
    pyomo \
    pysam \
    matplotlib \
    tables \
    pandas

#installing optitype form git repository (version Dec 09 2015) and wirtig config.ini
RUN git clone https://github.com/FRED-2/OptiType.git \
    && sed -i -e '1i#!/usr/bin/env python\' OptiType/OptiTypePipeline.py \
    && mv OptiType/ /usr/local/bin/ \
    && chmod 777 /usr/local/bin/OptiType/OptiTypePipeline.py \
    && echo "[mapping]\n\
razers3=/usr/local/bin/razers3 \n\
threads=1 \n\
\n\
[ilp]\n\
solver=cbc \n\
threads=1 \n\
\n\
[behavior]\n\
deletebam=true \n\
unpaired_weight=0 \n\
use_discordant=false\n" >> /usr/local/bin/OptiType/config.ini


#installing razers3
RUN git clone https://github.com/seqan/seqan.git seqan-src \
    && cd seqan-src \
    && cmake -DCMAKE_BUILD_TYPE=Release \
    && make razers3 \
    && cp bin/razers3 /usr/local/bin \
    && cd .. \
    && rm -rf seqan-src


#Polysolver


#Seq2HLA
RUN hg clone https://bitbucket.org/sebastian_boegel/seq2hla \
    && sed -i -e '1i#!/usr/bin/env python\' seq2hla/seq2HLA.py \
    && mv seq2hla/ /usr/local/bin/ \
    && chmod 777 /usr/local/bin/seq2hla/seq2HLA.py 
    

ENV PATH /usr/local/bin/OptiType:/usr/local/bin/seq2hla:$PATH

#Prediction Tools (generate a zip with all prediction tool binaries and wget it form our servers)
#UN wget http://abi.informatik.uni-tuebingen.de/pkg_predictors.tar.gz \

RUN wget --secure-protocol=TLSv1 http://www-bs3.informatik.uni-tuebingen.de/extern/schubert/pkg_predictors.tar.gz \
    && tar -xvf pkg_predictors.tar.gz  -C /usr/local/ \
    && rm pkg_predictors.tar.gz

#set envirnomental variables for prediction methods
ENV NETCHOP /usr/local/predictors/netchop/netchop-3.1
ENV TMPDIR /tmp
ENV PATH /usr/local/predictors/bin:$PATH

#Fred2
RUN pip install git+https://github.com/FRED-2/Fred2@master


#Get Fred2 COMMANDLINE TOOLS
RUN git clone https://github.com/FRED-2/ImmunoNodes.git \
    && chmod -R 777 /ImmunoNodes/src/

ENV PATH /ImmunoNodes/src/:$PATH

