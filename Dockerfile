FROM ubuntu:latest

RUN apt-get update && \
apt-get install --no-install-recommends -y \
unzip wget git parallel  \
python3.8 python3-pip python3.8-dev \
&& rm -rf /var/lib/apt/lists/*

RUN pip3 install scipy h5py

CMD alias python3=/usr/bin/python3.8

RUN git clone https://github.com/getian107/PRScsx.git
WORKDIR /PRScsx