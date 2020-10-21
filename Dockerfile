# syntax=docker/dockerfile:experimental
FROM archlinux/base as git

RUN pacman --noconfirm -Sy && \
pacman --noconfirm -S --needed archlinux-keyring && \
pacman --noconfirm --needed -S pacman pacman-mirrorlist sudo git base-devel \
    && rm -rf /etc/pacman.d/mirrorlist.pacnew \
    && pacman-db-upgrade

# set locale
RUN echo 'en_US.UTF-8 UTF-8' >>/etc/locale.gen && \
locale-gen && \
echo 'LANG=en_US.UTF-8' >/etc/locale.conf && \
echo 'export LANG=en_US.UTF-8' >>/etc/profile && \
echo 'export LC_ALL=en_US.UTF-8' >>/etc/profile

RUN useradd -m -s /bin/bash bdgp && \
chmod u+s /sbin/unix_chkpwd /sbin/su /sbin/sudo /sbin/passwd && \
echo 'ALL ALL=(ALL) NOPASSWD:ALL' >>/etc/sudoers

USER bdgp

RUN cd /home/bdgp && git clone https://aur.archlinux.org/yay-bin.git && \
cd /home/bdgp/yay-bin && makepkg -sri --noconfirm

USER bdgp
RUN yay -Sy --noconfirm --needed llvm llvm-libs clang libffi git openssh && yay -Sc --noconfirm ||true
USER root

RUN curl https://sh.rustup.rs -sSf | sh -s -- -y && \
source $HOME/.cargo/env && \
rustup default stable

RUN mkdir -p -m 0600 ~/.ssh && ssh-keyscan gitlab.com >> ~/.ssh/known_hosts
ARG project_id
ARG token
ADD "https://gitlab.com/api/v4/projects/${project_id}/repository/branches/master?private_token=${token}" .invalidateCache

RUN --mount=type=ssh cd /root && \
git clone -q "git@gitlab.com:bdgp/simsplice.git" && \
cd simsplice && \
source $HOME/.cargo/env && \
cargo build

FROM archlinux/base as install

RUN pacman --noconfirm -Sy && \
pacman --noconfirm -S --needed archlinux-keyring && \
pacman --noconfirm --needed -S pacman pacman-mirrorlist sudo git base-devel \
    && rm -rf /etc/pacman.d/mirrorlist.pacnew \
    && pacman-db-upgrade

# set locale
RUN echo 'en_US.UTF-8 UTF-8' >>/etc/locale.gen && \
locale-gen && \
echo 'LANG=en_US.UTF-8' >/etc/locale.conf && \
echo 'export LANG=en_US.UTF-8' >>/etc/profile && \
echo 'export LC_ALL=en_US.UTF-8' >>/etc/profile

RUN useradd -m -s /bin/bash bdgp && \
chmod u+s /sbin/unix_chkpwd /sbin/su /sbin/sudo /sbin/passwd && \
echo 'ALL ALL=(ALL) NOPASSWD:ALL' >>/etc/sudoers

USER bdgp
RUN yay -Sy --noconfirm --needed samtools bcftools pigz parallel htop rsync && yay -Sc --noconfirm ||true
USER root

RUN rsync -avP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/{bedGraphToBigWig,bigWigToBedGraph,faToTwoBit,bedToBigBed,bigBedToBed} /usr/bin

RUN --mount=type=bind,target=/root/simsplice,source=/root/simsplice,from=git,rw \
cp -v /root/simsplice/target/debug/{simsplice,genvcf,liftover} /root/simsplice/*.sh /usr/bin
