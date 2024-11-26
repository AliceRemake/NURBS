/**
  ******************************************************************************
  * @file           : TestFindSpan.cpp
  * @author         : AliceRemake
  * @brief          : None
  * @attention      : None
  * @date           : 24-11-21
  ******************************************************************************
  */



#include <doctest/doctest.h>
#include <NURBS.h>
#include <tinynurbs/tinynurbs.h>

TEST_CASE("FindSpan")
{
    constexpr size_t degree = 2;
    const std::vector knots = {0.0f, 0.0f, 0.0f, 1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 5.0f, 5.0f};
    CHECK(NURBS::FindSpan(degree, knots, 0.0f) == tinynurbs::findSpan(degree, knots, 0.0f));
    CHECK(NURBS::FindSpan(degree, knots, 0.5f) == tinynurbs::findSpan(degree, knots, 0.5f));
    CHECK(NURBS::FindSpan(degree, knots, 1.5f) == tinynurbs::findSpan(degree, knots, 1.5f));
    CHECK(NURBS::FindSpan(degree, knots, 2.5f) == tinynurbs::findSpan(degree, knots, 2.5f));
    CHECK(NURBS::FindSpan(degree, knots, 3.5f) == tinynurbs::findSpan(degree, knots, 3.5f));
    CHECK(NURBS::FindSpan(degree, knots, 4.5f) == tinynurbs::findSpan(degree, knots, 4.5f));
    CHECK(NURBS::FindSpan(degree, knots, 5.0f) == tinynurbs::findSpan(degree, knots, 5.0f));
}
