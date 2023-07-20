#ifdef __linux__

#pragma once

#include <filesystem>
#include <fstream>
#include <map>
#include <stdexcept>
#include <string>

#include <cmath>
#include <cstdint>

#include <unistd.h>

#include "memory_monitor.hpp"

class LinuxMemoryMonitor final : public MemoryMonitor
{
public:
    LinuxMemoryMonitor()
    {
        std::filesystem::path path("/proc");
        path /= std::to_string(getpid());
        path /= "status";
        m_path_str = path.string();
        std::cout << m_path_str << std::endl;
    }

    ~LinuxMemoryMonitor() = default;

    // non-copyable
    LinuxMemoryMonitor(const LinuxMemoryMonitor&) = delete;
    LinuxMemoryMonitor& operator=(const LinuxMemoryMonitor&) = delete;
    LinuxMemoryMonitor(LinuxMemoryMonitor&&) = delete;
    LinuxMemoryMonitor& operator=(LinuxMemoryMonitor&&) = delete;

    uint64_t Usage(MemoryType memory_type) const override
    {
        uint64_t usage = 0;
        std::ifstream file(m_path_str);
        std::string line;

        std::string const memory_type_str = memory_type_map.at(memory_type);
        while (std::getline(file, line))
        {
            if (line.find(memory_type_str) != std::string::npos)
            {
                usage = Bytes(line);
                /*if (memory_type == MemoryType::Physical)
                {
                    usage *= page_size;
                }*/
                break;
            }
        }
        return usage;
    }

private:
    std::string m_path_str;

    inline static std::map<MemoryType, std::string> const memory_type_map = {
        {MemoryType::Physical, "VmRSS"},
        {MemoryType::PhysicalPeak, "VmHWM"},
        {MemoryType::Virtual, "VmSize"},
        {MemoryType::VirtualPeak, "VmPeak"} //
    };

    inline static std::map<char, uint64_t> unit_prefix_map = {
        {'k', 1}, // kilo
        {'K', 1}, // Kilo
        {'m', 2}, // mega
        {'M', 2}, // Mega
        {'g', 3}, // giga
        {'G', 3}  // Giga
    };

    inline static uint64_t constexpr kilobyte = 1024;

    inline static uint64_t const page_size = sysconf(_SC_PAGE_SIZE); // in bytes

    static uint64_t ConversionFactor(std::string const& unit)
    {
        uint64_t conversion_factor = 1;
        char const byte = unit[1];
        if (byte == 'b' || byte == 'B')
        {
            char const unit_prefix = unit[0];
            conversion_factor = std::pow(kilobyte, unit_prefix_map.at(unit_prefix));
        }
        else
        {
            throw std::runtime_error("Memory is not measured bytes");
        }
        return conversion_factor;
    }

    // ex: VmPeak:     4296 kB
    //  KEY:            VALUE UNIT
    // >   < key
    //    > < colon
    //     >    < tab
    //         >        < spaces
    //                 >     < value
    //                      > < space
    //                        >  < unit (2 chars)
    // position of colon separating the key from the value

    static uint64_t Bytes(std::string& line)
    {
        size_t const key_end_pos = line.find(':');
        size_t const value_start_pos = line.find_first_not_of('\t', key_end_pos + 1);
        size_t const value_end_pos = line.size() - 2;
        std::string const value(line.begin() + value_start_pos, line.begin() + value_end_pos);
        uint64_t const usage = std::stoull(value.c_str());
        std::string const unit(line.begin() + value_end_pos, line.end());
        uint64_t const conversion_factor = ConversionFactor(unit);
        return usage * conversion_factor;
    }
};

#endif
