/**
  ******************************************************************************
  * @file           : TestCurvePoint.cpp
  * @author         : AliceRemake
  * @brief          : None
  * @attention      : None
  * @date           : 24-11-22
  ******************************************************************************
  */



#include <doctest/doctest.h>
#include <NURBS.h>
#include <tinynurbs/tinynurbs.h>

#define CHECK_GLM_VERTEX(lhs, rhs) CHECK((glm::distance(lhs, rhs) < 2 * std::numeric_limits<float>::epsilon()))

TEST_CASE("CurvePoint")
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

    CHECK_GLM_VERTEX(NURBS::CurvePoint(crv, 0.0f), tinynurbs::curvePoint(crv, 0.0f));
    CHECK_GLM_VERTEX(NURBS::CurvePoint(crv, 0.1f), tinynurbs::curvePoint(crv, 0.1f));
    CHECK_GLM_VERTEX(NURBS::CurvePoint(crv, 0.2f), tinynurbs::curvePoint(crv, 0.2f));
    CHECK_GLM_VERTEX(NURBS::CurvePoint(crv, 0.3f), tinynurbs::curvePoint(crv, 0.3f));
    CHECK_GLM_VERTEX(NURBS::CurvePoint(crv, 0.4f), tinynurbs::curvePoint(crv, 0.4f));
    CHECK_GLM_VERTEX(NURBS::CurvePoint(crv, 0.5f), tinynurbs::curvePoint(crv, 0.5f));
    CHECK_GLM_VERTEX(NURBS::CurvePoint(crv, 0.6f), tinynurbs::curvePoint(crv, 0.6f));
    CHECK_GLM_VERTEX(NURBS::CurvePoint(crv, 0.7f), tinynurbs::curvePoint(crv, 0.7f));
    CHECK_GLM_VERTEX(NURBS::CurvePoint(crv, 0.8f), tinynurbs::curvePoint(crv, 0.8f));
    CHECK_GLM_VERTEX(NURBS::CurvePoint(crv, 0.9f), tinynurbs::curvePoint(crv, 0.9f));
    CHECK_GLM_VERTEX(NURBS::CurvePoint(crv, 1.0f), tinynurbs::curvePoint(crv, 1.0f));
    
}
