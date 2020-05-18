# syntax=docker/dockerfile:experimental
FROM boss:6000/archdev as git

RUN pacman -Sy --noconfirm --needed llvm llvm-libs clang libffi && pacman -Sc --noconfirm ||true

#RUN ln -s /usr/lib/libffi.so.7 /usr/lib/libffi.so.6

RUN  mkdir -p /root/.ssh && chmod 700 /root/.ssh && ssh-keyscan -t rsa gitlab.com >>/root/.ssh/known_hosts

RUN --mount=type=ssh cd /root && \
git clone -q "git@gitlab.com:bdgp/simsplice.git" && \
cd simsplice && \
curl https://sh.rustup.rs -sSf | sh -s -- -y && \
source $HOME/.cargo/env && \
rustup default stable && \
cargo build --release #test

FROM boss:6000/archdev as install

RUN --mount=type=bind,target=/root/simsplice,source=/root/simsplice,from=git,rw \
cp -v /root/simsplice/target/release/{simsplice,genvcf} /usr/bin
