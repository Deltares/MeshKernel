#pragma once

#include <iostream>
#include <memory>

#if defined(_WIN32)
#include "memory_monitor_windows.hpp"
#elif defined(__linux__)
#include "memory_monitor_linux.hpp"
#else
#error "Unsupported platform"
#endif

inline static std::unique_ptr<MemoryMonitor> CreateMemoryMonitor()
{
#if defined(_WIN32)
    return std::make_unique<WindowsMemoryMonitor>();
#elif defined(__linux__)
    return std::make_unique<LinuxMemoryMonitor>();
#else
#error "Unsupported platform"
#endif
}