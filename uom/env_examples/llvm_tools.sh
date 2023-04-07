#!/bin/bash
# Copyright @2023 UOM project. All rights reserved.

fail() {
  >&2 echo "$@"
  exit 1
}

[[ -a "external/llvm_toolchain_llvm/bin/clangd" ]] || fail "bin/clangd not found"
[[ -a "external/llvm_toolchain_llvm/bin/clang-format" ]] || fail "bin/clang-format not found"

echo "SUCCESS!"
