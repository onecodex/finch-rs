# This Dockerfile builds the contents of this repo in a Rust builder container, then copies the
# resulting binaries into a Debian container for use.

# BUILDER
FROM ekidd/rust-musl-builder:stable AS builder

RUN mkdir /home/rust/finch
ADD --chown=rust:rust src/ /home/rust/finch/src/
ADD --chown=rust:rust Cargo.toml /home/rust/finch/Cargo.toml
ADD --chown=rust:rust rust-toolchain /home/rust/finch/rust-toolchain

RUN cd /home/rust/finch \
    && cargo +stable build --release --target x86_64-unknown-linux-musl

# MAIN CONTAINER
FROM python:3.7-slim-stretch

# install system dependencies
RUN apt-get update && apt-get install -y \
    curl \
    ssh \
    vim \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# copy finch from builder
COPY --from=builder \
    /home/rust/finch/target/x86_64-unknown-linux-musl/release/finch* \
    /usr/local/bin/

CMD /bin/bash