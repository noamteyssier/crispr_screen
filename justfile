
install:
  cargo install --path .

build:
  cargo build --release

clean:
  cargo clean

test:
  cargo nextest run

lint: build
  cargo clippy -- \
    -W clippy::pedantic \
    -A clippy::module_name_repetitions
