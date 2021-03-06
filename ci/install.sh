set -euxo pipefail

main() {
    if [ $TARGET = cargo-fmt ]; then
        rustup component add rustfmt-preview
        return
    fi

    rustup target add x86_64-unknown-linux-musl

    if [ $TARGET != x86_64-unknown-linux-gnu ]; then
        rustup target add $TARGET
    fi
}

main
