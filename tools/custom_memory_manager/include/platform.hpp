#pragma once

#if defined(_WIN32) && defined(_MSC_VER)
#define WIN_MSVC_BENCHMARK
#elif (defined(__unix__) && defined(__GNUC__)) || \
    ((defined(__unix__) || defined(_WIN32)) && defined(__clang__))
#define LINUX_GNUC_BENCHMARK
#else
#error "Unsupported platform and/or compiler"
#endif