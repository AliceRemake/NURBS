/**
  ******************************************************************************
  * @file           : TestCurveDerivative.cpp
  * @author         : AliceRemake
  * @brief          : None
  * @attention      : None
  * @date           : 24-11-26
  ******************************************************************************
  */



#include <doctest/doctest.h>
#include <NURBS.h>
#include <tinynurbs/tinynurbs.h>

#define CHECK_VERTEX_VECTOR(lhs, rhs)                                                     \
    do {                                                                                  \
        auto l = lhs;                                                                     \
        auto r = rhs;                                                                     \
        CHECK(l.size() == r.size());                                                      \
        for (size_t i = 0; i < l.size(); ++i)                                             \
        {                                                                                 \
            CHECK(glm::distance(l[i], r[i]) < 2 * std::numeric_limits<float>::epsilon()); \
        }                                                                                 \
    } while(0)

TEST_CASE("CurveDerivatives")
{
    constexpr size_t degree = 2;
    const std::vector knots = { 0.0f, 0.0f, 0.0f, 1.0f, 1.0f, 1.0f };
    const std::vector control_points = {
        glm::vec3(-1, 0, 0),
        glm::vec3( 0, 1, 0),
        glm::vec3( 1, 0, 0),
    };
    const std::vector weights = { 1.0f, 2.0f, 3.0f };
    tinynurbs::RationalCurve crv(
        degree,
        knots,
        control_points,
        weights
    );

    CHECK_VERTEX_VECTOR(NURBS::CurveDerivatives(crv, 0, 0.0f), tinynurbs::curveDerivatives(crv, 0, 0.0f));
    CHECK_VERTEX_VECTOR(NURBS::CurveDerivatives(crv, 0, 0.1f), tinynurbs::curveDerivatives(crv, 0, 0.1f));
    CHECK_VERTEX_VECTOR(NURBS::CurveDerivatives(crv, 0, 0.2f), tinynurbs::curveDerivatives(crv, 0, 0.2f));
    CHECK_VERTEX_VECTOR(NURBS::CurveDerivatives(crv, 0, 0.3f), tinynurbs::curveDerivatives(crv, 0, 0.3f));
    CHECK_VERTEX_VECTOR(NURBS::CurveDerivatives(crv, 0, 0.4f), tinynurbs::curveDerivatives(crv, 0, 0.4f));
    CHECK_VERTEX_VECTOR(NURBS::CurveDerivatives(crv, 0, 0.5f), tinynurbs::curveDerivatives(crv, 0, 0.5f));
    CHECK_VERTEX_VECTOR(NURBS::CurveDerivatives(crv, 0, 0.6f), tinynurbs::curveDerivatives(crv, 0, 0.6f));
    CHECK_VERTEX_VECTOR(NURBS::CurveDerivatives(crv, 0, 0.7f), tinynurbs::curveDerivatives(crv, 0, 0.7f));
    CHECK_VERTEX_VECTOR(NURBS::CurveDerivatives(crv, 0, 0.8f), tinynurbs::curveDerivatives(crv, 0, 0.8f));
    CHECK_VERTEX_VECTOR(NURBS::CurveDerivatives(crv, 0, 0.9f), tinynurbs::curveDerivatives(crv, 0, 0.9f));
    CHECK_VERTEX_VECTOR(NURBS::CurveDerivatives(crv, 0, 1.0f), tinynurbs::curveDerivatives(crv, 0, 1.0f));

    CHECK_VERTEX_VECTOR(NURBS::CurveDerivatives(crv, 1, 0.0f), tinynurbs::curveDerivatives(crv, 1, 0.0f));
    CHECK_VERTEX_VECTOR(NURBS::CurveDerivatives(crv, 1, 0.1f), tinynurbs::curveDerivatives(crv, 1, 0.1f));
    CHECK_VERTEX_VECTOR(NURBS::CurveDerivatives(crv, 1, 0.2f), tinynurbs::curveDerivatives(crv, 1, 0.2f));
    CHECK_VERTEX_VECTOR(NURBS::CurveDerivatives(crv, 1, 0.3f), tinynurbs::curveDerivatives(crv, 1, 0.3f));
    CHECK_VERTEX_VECTOR(NURBS::CurveDerivatives(crv, 1, 0.4f), tinynurbs::curveDerivatives(crv, 1, 0.4f));
    CHECK_VERTEX_VECTOR(NURBS::CurveDerivatives(crv, 1, 0.5f), tinynurbs::curveDerivatives(crv, 1, 0.5f));
    CHECK_VERTEX_VECTOR(NURBS::CurveDerivatives(crv, 1, 0.6f), tinynurbs::curveDerivatives(crv, 1, 0.6f));
    CHECK_VERTEX_VECTOR(NURBS::CurveDerivatives(crv, 1, 0.7f), tinynurbs::curveDerivatives(crv, 1, 0.7f));
    CHECK_VERTEX_VECTOR(NURBS::CurveDerivatives(crv, 1, 0.8f), tinynurbs::curveDerivatives(crv, 1, 0.8f));
    CHECK_VERTEX_VECTOR(NURBS::CurveDerivatives(crv, 1, 0.9f), tinynurbs::curveDerivatives(crv, 1, 0.9f));
    CHECK_VERTEX_VECTOR(NURBS::CurveDerivatives(crv, 1, 1.0f), tinynurbs::curveDerivatives(crv, 1, 1.0f));

    CHECK_VERTEX_VECTOR(NURBS::CurveDerivatives(crv, 2, 0.0f), tinynurbs::curveDerivatives(crv, 2, 0.0f));
    CHECK_VERTEX_VECTOR(NURBS::CurveDerivatives(crv, 2, 0.1f), tinynurbs::curveDerivatives(crv, 2, 0.1f));
    CHECK_VERTEX_VECTOR(NURBS::CurveDerivatives(crv, 2, 0.2f), tinynurbs::curveDerivatives(crv, 2, 0.2f));
    CHECK_VERTEX_VECTOR(NURBS::CurveDerivatives(crv, 2, 0.3f), tinynurbs::curveDerivatives(crv, 2, 0.3f));
    CHECK_VERTEX_VECTOR(NURBS::CurveDerivatives(crv, 2, 0.4f), tinynurbs::curveDerivatives(crv, 2, 0.4f));
    CHECK_VERTEX_VECTOR(NURBS::CurveDerivatives(crv, 2, 0.5f), tinynurbs::curveDerivatives(crv, 2, 0.5f));
    CHECK_VERTEX_VECTOR(NURBS::CurveDerivatives(crv, 2, 0.6f), tinynurbs::curveDerivatives(crv, 2, 0.6f));
    CHECK_VERTEX_VECTOR(NURBS::CurveDerivatives(crv, 2, 0.7f), tinynurbs::curveDerivatives(crv, 2, 0.7f));
    CHECK_VERTEX_VECTOR(NURBS::CurveDerivatives(crv, 2, 0.8f), tinynurbs::curveDerivatives(crv, 2, 0.8f));
    CHECK_VERTEX_VECTOR(NURBS::CurveDerivatives(crv, 2, 0.9f), tinynurbs::curveDerivatives(crv, 2, 0.9f));
    CHECK_VERTEX_VECTOR(NURBS::CurveDerivatives(crv, 2, 1.0f), tinynurbs::curveDerivatives(crv, 2, 1.0f));

    CHECK_VERTEX_VECTOR(NURBS::CurveDerivatives(crv, 3, 0.0f), tinynurbs::curveDerivatives(crv, 3, 0.0f));
    CHECK_VERTEX_VECTOR(NURBS::CurveDerivatives(crv, 3, 0.1f), tinynurbs::curveDerivatives(crv, 3, 0.1f));
    CHECK_VERTEX_VECTOR(NURBS::CurveDerivatives(crv, 3, 0.2f), tinynurbs::curveDerivatives(crv, 3, 0.2f));
    CHECK_VERTEX_VECTOR(NURBS::CurveDerivatives(crv, 3, 0.3f), tinynurbs::curveDerivatives(crv, 3, 0.3f));
    CHECK_VERTEX_VECTOR(NURBS::CurveDerivatives(crv, 3, 0.4f), tinynurbs::curveDerivatives(crv, 3, 0.4f));
    CHECK_VERTEX_VECTOR(NURBS::CurveDerivatives(crv, 3, 0.5f), tinynurbs::curveDerivatives(crv, 3, 0.5f));
    CHECK_VERTEX_VECTOR(NURBS::CurveDerivatives(crv, 3, 0.6f), tinynurbs::curveDerivatives(crv, 3, 0.6f));
    CHECK_VERTEX_VECTOR(NURBS::CurveDerivatives(crv, 3, 0.7f), tinynurbs::curveDerivatives(crv, 3, 0.7f));
    CHECK_VERTEX_VECTOR(NURBS::CurveDerivatives(crv, 3, 0.8f), tinynurbs::curveDerivatives(crv, 3, 0.8f));
    CHECK_VERTEX_VECTOR(NURBS::CurveDerivatives(crv, 3, 0.9f), tinynurbs::curveDerivatives(crv, 3, 0.9f));
    CHECK_VERTEX_VECTOR(NURBS::CurveDerivatives(crv, 3, 1.0f), tinynurbs::curveDerivatives(crv, 3, 1.0f));
}
