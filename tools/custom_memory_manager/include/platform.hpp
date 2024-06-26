#pragma once

#if defined(_WIN32) && defined(_MSC_VER)
#define WIN_MSVC_BENCHMARK
#elif (defined(__linux__) || defined(__APPLE__)) && defined(__GNUC__)
#define LINUX_OR_APPLE_GNUC_BENCHMARK
#else
#error "Unsupported platform and/or compiler"
#endif
