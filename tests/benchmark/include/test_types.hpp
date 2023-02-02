#pragma once

class Point
{
public:
    Point() = default;
    Point(int x, int y, int z) : m_x{x}, m_y{y}, m_z{z} {}
    ~Point() = default;
    Point(Point const&) = default;
    Point& operator=(Point const&) = default;

private:
    int m_x = -1; // 4 bytes
    int m_y = -2; // 4 bytes
    int m_z = -3; // 4 bytes
};

static size_t constexpr size_of_point = sizeof(Point); // 3 x 4 bytes = 12 bytes