# syntax=docker/dockerfile:experimental
FROM ubuntu:focal as build
SHELL ["/bin/bash", "-lc"]

RUN apt-get -y update && \
DEBIAN_FRONTEND=noninteractive apt-get install -y git openssh-client build-essential \
  llvm clang libffi-dev curl libssl-dev libssl1.1 openssl pkg-config

RUN  mkdir -p /root/.ssh && chmod 700 /root/.ssh && ssh-keyscan -t rsa gitlab.com >>/root/.ssh/known_hosts
ARG project_id
ARG token
ADD "https://gitlab.com/api/v4/projects/${project_id}/repository/branches/master?private_token=${token}" .invalidateCache

RUN --mount=type=ssh cd /root && \
git clone --depth 1 "git@gitlab.com:bdgp/simsplice.git" && \
cd simsplice && \
curl https://sh.rustup.rs -sSf | sh -s -- -y && \
source $HOME/.cargo/env && \
rustup default stable && \
cargo build --release

FROM ubuntu:focal as install
SHELL ["/bin/bash", "-lc"]

RUN apt-get -y update && \
DEBIAN_FRONTEND=noninteractive apt-get install -y rsync parallel pigz bcftools htop locales

ENV LC_ALL=en_US.UTF-8
ENV LANG=en_US.UTF-8
RUN locale-gen en_US.UTF-8

RUN rsync -avP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/{bedGraphToBigWig,bigWigToBedGraph,faToTwoBit,bedToBigBed,bigBedToBed} /usr/bin

RUN --mount=type=bind,target=/root/simsplice,source=/root/simsplice,from=build,rw \
cp -v /root/simsplice/target/release/{simsplice,genvcf,liftover} /root/simsplice/*.sh /usr/bin
