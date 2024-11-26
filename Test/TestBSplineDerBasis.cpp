/**
  ******************************************************************************
  * @file           : TestBSplineDerBasis.cpp
  * @author         : AliceRemake
  * @brief          : None
  * @attention      : None
  * @date           : 24-11-22
  ******************************************************************************
  */



#include <doctest/doctest.h>
#include <NURBS.h>
#include <tinynurbs/tinynurbs.h>

#define CHECK_FLOAT_MATRIX(lhs, rhs)                                                        \
    do {                                                                                    \
        auto l = lhs;                                                                       \
        auto r = rhs;                                                                       \
        CHECK(l.size() == r.rows());                                                        \
        for (size_t i = 0; i < l.size(); ++i)                                               \
        {                                                                                   \
            CHECK(l[i].size() == r.cols());                                                 \
            for (size_t j = 0; j < l[i].size(); ++j)                                        \
            {                                                                               \
                CHECK(std::fabs(l[i][j] - r(i,j)) < std::numeric_limits<float>::epsilon()); \
            }                                                                               \
        }                                                                                   \
    } while(0)

TEST_CASE("BSplineDerBasis")
{
    constexpr size_t degree = 2;
    const std::vector knots = {0.0f, 0.0f, 0.0f, 1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 5.0f, 5.0f};
    CHECK_FLOAT_MATRIX(NURBS::BSplineDerBasis(degree, NURBS::FindSpan(degree, knots, 0.0f), knots, 0.0f, 0), tinynurbs::bsplineDerBasis(degree, tinynurbs::findSpan(degree, knots, 0.0f), knots, 0.0f, 0));
    CHECK_FLOAT_MATRIX(NURBS::BSplineDerBasis(degree, NURBS::FindSpan(degree, knots, 0.5f), knots, 0.5f, 0), tinynurbs::bsplineDerBasis(degree, tinynurbs::findSpan(degree, knots, 0.5f), knots, 0.5f, 0));
    CHECK_FLOAT_MATRIX(NURBS::BSplineDerBasis(degree, NURBS::FindSpan(degree, knots, 1.5f), knots, 1.5f, 0), tinynurbs::bsplineDerBasis(degree, tinynurbs::findSpan(degree, knots, 1.5f), knots, 1.5f, 0));
    CHECK_FLOAT_MATRIX(NURBS::BSplineDerBasis(degree, NURBS::FindSpan(degree, knots, 2.5f), knots, 2.5f, 0), tinynurbs::bsplineDerBasis(degree, tinynurbs::findSpan(degree, knots, 2.5f), knots, 2.5f, 0));
    CHECK_FLOAT_MATRIX(NURBS::BSplineDerBasis(degree, NURBS::FindSpan(degree, knots, 3.5f), knots, 3.5f, 0), tinynurbs::bsplineDerBasis(degree, tinynurbs::findSpan(degree, knots, 3.5f), knots, 3.5f, 0));
    CHECK_FLOAT_MATRIX(NURBS::BSplineDerBasis(degree, NURBS::FindSpan(degree, knots, 4.5f), knots, 4.5f, 0), tinynurbs::bsplineDerBasis(degree, tinynurbs::findSpan(degree, knots, 4.5f), knots, 4.5f, 0));
    CHECK_FLOAT_MATRIX(NURBS::BSplineDerBasis(degree, NURBS::FindSpan(degree, knots, 5.0f), knots, 5.0f, 0), tinynurbs::bsplineDerBasis(degree, tinynurbs::findSpan(degree, knots, 5.0f), knots, 5.0f, 0));

    CHECK_FLOAT_MATRIX(NURBS::BSplineDerBasis(degree, NURBS::FindSpan(degree, knots, 0.0f), knots, 0.0f, 1), tinynurbs::bsplineDerBasis(degree, tinynurbs::findSpan(degree, knots, 0.0f), knots, 0.0f, 1));
    CHECK_FLOAT_MATRIX(NURBS::BSplineDerBasis(degree, NURBS::FindSpan(degree, knots, 0.5f), knots, 0.5f, 1), tinynurbs::bsplineDerBasis(degree, tinynurbs::findSpan(degree, knots, 0.5f), knots, 0.5f, 1));
    CHECK_FLOAT_MATRIX(NURBS::BSplineDerBasis(degree, NURBS::FindSpan(degree, knots, 1.5f), knots, 1.5f, 1), tinynurbs::bsplineDerBasis(degree, tinynurbs::findSpan(degree, knots, 1.5f), knots, 1.5f, 1));
    CHECK_FLOAT_MATRIX(NURBS::BSplineDerBasis(degree, NURBS::FindSpan(degree, knots, 2.5f), knots, 2.5f, 1), tinynurbs::bsplineDerBasis(degree, tinynurbs::findSpan(degree, knots, 2.5f), knots, 2.5f, 1));
    CHECK_FLOAT_MATRIX(NURBS::BSplineDerBasis(degree, NURBS::FindSpan(degree, knots, 3.5f), knots, 3.5f, 1), tinynurbs::bsplineDerBasis(degree, tinynurbs::findSpan(degree, knots, 3.5f), knots, 3.5f, 1));
    CHECK_FLOAT_MATRIX(NURBS::BSplineDerBasis(degree, NURBS::FindSpan(degree, knots, 4.5f), knots, 4.5f, 1), tinynurbs::bsplineDerBasis(degree, tinynurbs::findSpan(degree, knots, 4.5f), knots, 4.5f, 1));
    CHECK_FLOAT_MATRIX(NURBS::BSplineDerBasis(degree, NURBS::FindSpan(degree, knots, 5.0f), knots, 5.0f, 1), tinynurbs::bsplineDerBasis(degree, tinynurbs::findSpan(degree, knots, 5.0f), knots, 5.0f, 1));

    CHECK_FLOAT_MATRIX(NURBS::BSplineDerBasis(degree, NURBS::FindSpan(degree, knots, 0.0f), knots, 0.0f, 2), tinynurbs::bsplineDerBasis(degree, tinynurbs::findSpan(degree, knots, 0.0f), knots, 0.0f, 2));
    CHECK_FLOAT_MATRIX(NURBS::BSplineDerBasis(degree, NURBS::FindSpan(degree, knots, 0.5f), knots, 0.5f, 2), tinynurbs::bsplineDerBasis(degree, tinynurbs::findSpan(degree, knots, 0.5f), knots, 0.5f, 2));
    CHECK_FLOAT_MATRIX(NURBS::BSplineDerBasis(degree, NURBS::FindSpan(degree, knots, 1.5f), knots, 1.5f, 2), tinynurbs::bsplineDerBasis(degree, tinynurbs::findSpan(degree, knots, 1.5f), knots, 1.5f, 2));
    CHECK_FLOAT_MATRIX(NURBS::BSplineDerBasis(degree, NURBS::FindSpan(degree, knots, 2.5f), knots, 2.5f, 2), tinynurbs::bsplineDerBasis(degree, tinynurbs::findSpan(degree, knots, 2.5f), knots, 2.5f, 2));
    CHECK_FLOAT_MATRIX(NURBS::BSplineDerBasis(degree, NURBS::FindSpan(degree, knots, 3.5f), knots, 3.5f, 2), tinynurbs::bsplineDerBasis(degree, tinynurbs::findSpan(degree, knots, 3.5f), knots, 3.5f, 2));
    CHECK_FLOAT_MATRIX(NURBS::BSplineDerBasis(degree, NURBS::FindSpan(degree, knots, 4.5f), knots, 4.5f, 2), tinynurbs::bsplineDerBasis(degree, tinynurbs::findSpan(degree, knots, 4.5f), knots, 4.5f, 2));
    CHECK_FLOAT_MATRIX(NURBS::BSplineDerBasis(degree, NURBS::FindSpan(degree, knots, 5.0f), knots, 5.0f, 2), tinynurbs::bsplineDerBasis(degree, tinynurbs::findSpan(degree, knots, 5.0f), knots, 5.0f, 2));

    CHECK_FLOAT_MATRIX(NURBS::BSplineDerBasis(degree, NURBS::FindSpan(degree, knots, 0.0f), knots, 0.0f, 3), tinynurbs::bsplineDerBasis(degree, tinynurbs::findSpan(degree, knots, 0.0f), knots, 0.0f, 3));
    CHECK_FLOAT_MATRIX(NURBS::BSplineDerBasis(degree, NURBS::FindSpan(degree, knots, 0.5f), knots, 0.5f, 3), tinynurbs::bsplineDerBasis(degree, tinynurbs::findSpan(degree, knots, 0.5f), knots, 0.5f, 3));
    CHECK_FLOAT_MATRIX(NURBS::BSplineDerBasis(degree, NURBS::FindSpan(degree, knots, 1.5f), knots, 1.5f, 3), tinynurbs::bsplineDerBasis(degree, tinynurbs::findSpan(degree, knots, 1.5f), knots, 1.5f, 3));
    CHECK_FLOAT_MATRIX(NURBS::BSplineDerBasis(degree, NURBS::FindSpan(degree, knots, 2.5f), knots, 2.5f, 3), tinynurbs::bsplineDerBasis(degree, tinynurbs::findSpan(degree, knots, 2.5f), knots, 2.5f, 3));
    CHECK_FLOAT_MATRIX(NURBS::BSplineDerBasis(degree, NURBS::FindSpan(degree, knots, 3.5f), knots, 3.5f, 3), tinynurbs::bsplineDerBasis(degree, tinynurbs::findSpan(degree, knots, 3.5f), knots, 3.5f, 3));
    CHECK_FLOAT_MATRIX(NURBS::BSplineDerBasis(degree, NURBS::FindSpan(degree, knots, 4.5f), knots, 4.5f, 3), tinynurbs::bsplineDerBasis(degree, tinynurbs::findSpan(degree, knots, 4.5f), knots, 4.5f, 3));
    CHECK_FLOAT_MATRIX(NURBS::BSplineDerBasis(degree, NURBS::FindSpan(degree, knots, 5.0f), knots, 5.0f, 3), tinynurbs::bsplineDerBasis(degree, tinynurbs::findSpan(degree, knots, 5.0f), knots, 5.0f, 3));

    CHECK_FLOAT_MATRIX(NURBS::BSplineDerBasis(degree, NURBS::FindSpan(degree, knots, 0.0f), knots, 0.0f, 4), tinynurbs::bsplineDerBasis(degree, tinynurbs::findSpan(degree, knots, 0.0f), knots, 0.0f, 4));
    CHECK_FLOAT_MATRIX(NURBS::BSplineDerBasis(degree, NURBS::FindSpan(degree, knots, 0.5f), knots, 0.5f, 4), tinynurbs::bsplineDerBasis(degree, tinynurbs::findSpan(degree, knots, 0.5f), knots, 0.5f, 4));
    CHECK_FLOAT_MATRIX(NURBS::BSplineDerBasis(degree, NURBS::FindSpan(degree, knots, 1.5f), knots, 1.5f, 4), tinynurbs::bsplineDerBasis(degree, tinynurbs::findSpan(degree, knots, 1.5f), knots, 1.5f, 4));
    CHECK_FLOAT_MATRIX(NURBS::BSplineDerBasis(degree, NURBS::FindSpan(degree, knots, 2.5f), knots, 2.5f, 4), tinynurbs::bsplineDerBasis(degree, tinynurbs::findSpan(degree, knots, 2.5f), knots, 2.5f, 4));
    CHECK_FLOAT_MATRIX(NURBS::BSplineDerBasis(degree, NURBS::FindSpan(degree, knots, 3.5f), knots, 3.5f, 4), tinynurbs::bsplineDerBasis(degree, tinynurbs::findSpan(degree, knots, 3.5f), knots, 3.5f, 4));
    CHECK_FLOAT_MATRIX(NURBS::BSplineDerBasis(degree, NURBS::FindSpan(degree, knots, 4.5f), knots, 4.5f, 4), tinynurbs::bsplineDerBasis(degree, tinynurbs::findSpan(degree, knots, 4.5f), knots, 4.5f, 4));
    CHECK_FLOAT_MATRIX(NURBS::BSplineDerBasis(degree, NURBS::FindSpan(degree, knots, 5.0f), knots, 5.0f, 4), tinynurbs::bsplineDerBasis(degree, tinynurbs::findSpan(degree, knots, 5.0f), knots, 5.0f, 4));

}
