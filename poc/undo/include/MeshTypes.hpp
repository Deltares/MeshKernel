#ifndef MESH_TYPES__HPP
#define MESH_TYPES__HPP

namespace core
{
    constexpr size_t nullValueId = -999;

    constexpr double nullValueLoc = 1.0e20;

} // namespace core

class Point
{
public:
    Point(const double x, const double y) : x_(x), y_(y) {}
    Point() = default;

    double x() const;

    double& x();

    double y() const;

    double& y();

    bool isValid() const;

    void setInvalid();

private:
    double x_ = core::nullValueLoc;
    double y_ = core::nullValueLoc;
};

class Edge
{
public:
    enum Location
    {
        Start,
        End,
        Unknown
    };

    Edge(const size_t start, const size_t end) : start_(start), end_(end) {}
    Edge() = default;

    size_t start() const;

    size_t& start();

    size_t end() const;

    size_t& end();

    bool isValid() const;

    Location location(const size_t id) const;

private:
    size_t start_ = -1;
    size_t end_ = -1;
};

//--------------------------------

inline double Point::x() const
{
    return x_;
}

inline double& Point::x()
{
    return x_;
}

inline double Point::y() const
{
    return y_;
}

inline double& Point::y()
{
    return y_;
}

inline bool Point::isValid() const
{
    return x_ != core::nullValueLoc or y_ != core::nullValueLoc;
}

inline void Point::setInvalid()
{
    x_ = core::nullValueLoc;
    y_ = core::nullValueLoc;
}

inline size_t Edge::start() const
{
    return start_;
}

inline size_t& Edge::start()
{
    return start_;
}

inline size_t Edge::end() const
{
    return end_;
}

inline size_t& Edge::end()
{
    return end_;
}

inline bool Edge::isValid() const
{
    return start_ != core::nullValueId or end_ != core::nullValueId;
}

inline Edge::Location Edge::location(const size_t id) const
{
    if (id == start_)
    {
        return Start;
    }
    else if (id == end_)
    {
        return End;
    }
    else
    {
        return Unknown;
    }
}

#endif // MESH_TYPES__HPP
