/**
  ******************************************************************************
  * @file           : TestBSplineBasis.cpp
  * @author         : AliceRemake
  * @brief          : None
  * @attention      : None
  * @date           : 24-11-22
  ******************************************************************************
  */



#include <doctest/doctest.h>
#include <NURBS.h>
#include <tinynurbs/tinynurbs.h>

#define CHECK_FLOAT_VECTOR(lhs, rhs)                                               \
    do {                                                                           \
        auto l = lhs;                                                              \
        auto r = rhs;                                                              \
        CHECK(l.size() == r.size());                                               \
        for (size_t i = 0; i < l.size(); ++i)                                      \
        {                                                                          \
            CHECK(std::fabs(l[i] - r[i]) < std::numeric_limits<float>::epsilon()); \
        }                                                                          \
    } while(0)

TEST_CASE("BSplineBasis")
{
    constexpr size_t degree = 2;
    const std::vector knots = {0.0f, 0.0f, 0.0f, 1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 5.0f, 5.0f};
    CHECK_FLOAT_VECTOR(NURBS::BSplineBasis(degree, NURBS::FindSpan(degree, knots, 0.0f), knots, 0.0f), tinynurbs::bsplineBasis(degree, tinynurbs::findSpan(degree, knots, 0.0f), knots, 0.0f));
    CHECK_FLOAT_VECTOR(NURBS::BSplineBasis(degree, NURBS::FindSpan(degree, knots, 0.5f), knots, 0.5f), tinynurbs::bsplineBasis(degree, tinynurbs::findSpan(degree, knots, 0.5f), knots, 0.5f));
    CHECK_FLOAT_VECTOR(NURBS::BSplineBasis(degree, NURBS::FindSpan(degree, knots, 1.5f), knots, 1.5f), tinynurbs::bsplineBasis(degree, tinynurbs::findSpan(degree, knots, 1.5f), knots, 1.5f));
    CHECK_FLOAT_VECTOR(NURBS::BSplineBasis(degree, NURBS::FindSpan(degree, knots, 2.5f), knots, 2.5f), tinynurbs::bsplineBasis(degree, tinynurbs::findSpan(degree, knots, 2.5f), knots, 2.5f));
    CHECK_FLOAT_VECTOR(NURBS::BSplineBasis(degree, NURBS::FindSpan(degree, knots, 3.5f), knots, 3.5f), tinynurbs::bsplineBasis(degree, tinynurbs::findSpan(degree, knots, 3.5f), knots, 3.5f));
    CHECK_FLOAT_VECTOR(NURBS::BSplineBasis(degree, NURBS::FindSpan(degree, knots, 4.5f), knots, 4.5f), tinynurbs::bsplineBasis(degree, tinynurbs::findSpan(degree, knots, 4.5f), knots, 4.5f));
    CHECK_FLOAT_VECTOR(NURBS::BSplineBasis(degree, NURBS::FindSpan(degree, knots, 5.0f), knots, 5.0f), tinynurbs::bsplineBasis(degree, tinynurbs::findSpan(degree, knots, 5.0f), knots, 5.0f));
}
