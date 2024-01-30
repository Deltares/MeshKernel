#ifndef MESH_TYPES__HPP
#define MESH_TYPES__HPP

namespace core
{
    constexpr int nullValueId = -999;

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

    Edge(const int start, const int end) : start_(start), end_(end) {}
    Edge() = default;

    int start() const;

    int& start();

    int end() const;

    int& end();

    bool isValid() const;

    Location location(const int id) const;

private:
    int start_ = -1;
    int end_ = -1;
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

inline int Edge::start() const
{
    return start_;
}

inline int& Edge::start()
{
    return start_;
}

inline int Edge::end() const
{
    return end_;
}

inline int& Edge::end()
{
    return end_;
}

inline bool Edge::isValid() const
{
    return start_ != core::nullValueId or end_ != core::nullValueId;
}

inline Edge::Location Edge::location(const int id) const
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
