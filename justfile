
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
    -A clippy::module_name_repetitions \
    -A clippy::cast_precision_loss \
    -A clippy::cast_sign_loss \
    -A clippy::cast_possible_truncation
