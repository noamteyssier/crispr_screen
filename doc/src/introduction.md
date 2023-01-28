# [ `crispr_screen` ]

[![MIT licensed](https://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/noamteyssier/crispr_screen/blob/main/LICENSE)
[![actions status](https://github.com/noamteyssier/crispr_screen/workflows/CI/badge.svg)](https://github.com/noamteyssier/crispr_screen/actions)
[![codecov](https://codecov.io/gh/noamteyssier/crispr_screen/branch/main/graph/badge.svg?token=9ALCE60W2T)](https://codecov.io/gh/noamteyssier/crispr_screen)
[![Crates.io](https://img.shields.io/crates/d/crispr_screen?color=orange&label=crates.io)](https://crates.io/crates/crispr_screen)

## Introduction

`crispr_screen` is a free, open-source command-line tool that enables easy and
efficient differential expression analysis for CRISPR screens.

It is a faithful recreation of the `MAGeCK` algorithm, which is a form of the
`DESeq` method specifically applied for CRISPR screens. This tool extends the
`MAGeCK` algorithm with the `MAGeCK-INC` method, but also implements the `αRRA`
algorithm described in the original `MAGeCK` paper.

It is a drop-in replacement of current differential expression analyses that is
faster, more stable, and easier to use.
