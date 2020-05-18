# syntax=docker/dockerfile:experimental
FROM boss:6000/archdev as git

RUN pacman -Sy --noconfirm --needed llvm llvm-libs clang libffi git openssh && pacman -Sc --noconfirm ||true

RUN mkdir -p -m 0600 ~/.ssh && ssh-keyscan gitlab.com >> ~/.ssh/known_hosts

RUN --mount=type=ssh cd /root && \
git clone -q "git@gitlab.com:bdgp/simsplice.git" && \
cd simsplice && \
curl https://sh.rustup.rs -sSf | sh -s -- -y && \
source $HOME/.cargo/env && \
rustup default stable && \
cargo build --release #test3

FROM boss:6000/archdev as install

RUN --mount=type=bind,target=/root/simsplice,source=/root/simsplice,from=git,rw \
cp -v /root/simsplice/target/release/{simsplice,genvcf} /usr/bin
