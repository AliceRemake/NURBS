/**
  ******************************************************************************
  * @file           : TestSurfaceDerivatives.cpp
  * @author         : AliceRemake
  * @brief          : None
  * @attention      : None
  * @date           : 24-11-26
  ******************************************************************************
  */



#include <doctest/doctest.h>
#include <NURBS.h>
#include <tinynurbs/tinynurbs.h>

#define CHECK_VERTEX_MATRIX(lhs, rhs) \
    do {                                                                                           \
        auto l = lhs;                                                                              \
        auto r = rhs;                                                                              \
        CHECK(l.size() == r.rows());                                                               \
        for (size_t i = 0; i < l.size(); ++i)                                                      \
        {                                                                                          \
            CHECK(l[i].size() == r.cols());                                                        \
            for (size_t j = 0; j < l[i].size(); ++j)                                               \
            {                                                                                      \
                CHECK(glm::distance(l[i][j], r(i,j)) < 2 * std::numeric_limits<float>::epsilon()); \
            }                                                                                      \
        }                                                                                          \
    } while(0)

TEST_CASE("SurfaceDerivatives")
{
    tinynurbs::RationalSurface3f srf;
    srf.degree_u = 3;
    srf.degree_v = 3;
    srf.knots_u = {0, 0, 0, 0, 1, 1, 1, 1};
    srf.knots_v = {0, 0, 0, 0, 1, 1, 1, 1};
    // 4x4 grid (tinynurbs::array2) of control points and weights
    // https://www.geometrictools.com/Documentation/NURBSCircleSphere.pdf
    srf.control_points = {4, 4, 
                          {glm::vec3(0, 0, 1), glm::vec3(0, 0, 1), glm::vec3(0, 0, 1), glm::vec3(0, 0, 1),
                           glm::vec3(2, 0, 1), glm::vec3(2, 4, 1),  glm::vec3(-2, 4, 1),  glm::vec3(-2, 0, 1),
                           glm::vec3(2, 0, -1), glm::vec3(2, 4, -1), glm::vec3(-2, 4, -1), glm::vec3(-2, 0, -1),
                           glm::vec3(0, 0, -1), glm::vec3(0, 0, -1), glm::vec3(0, 0, -1), glm::vec3(0, 0, -1)
                          }
    };
    srf.weights = {4, 4,
                   {1,       1.f/3.f, 1.f/3.f, 1,
                    1.f/3.f, 1.f/9.f, 1.f/9.f, 1.f/3.f,
                    1.f/3.f, 1.f/9.f, 1.f/9.f, 1.f/3.f,
                    1,       1.f/3.f, 1.f/3.f, 1
                   }
    };

    for (const auto u : { 0.0f,0.1f,0.2f,0.3f,0.4f,0.5f,0.6f,0.7f,0.8f,0.9f,1.0f })
    {
        for (const auto v : { 0.0f,0.1f,0.2f,0.3f,0.4f,0.5f,0.6f,0.7f,0.8f,0.9f,1.0f })
        {
            CHECK_VERTEX_MATRIX(NURBS::SurfaceDerivatives(srf, 0, u, v), tinynurbs::surfaceDerivatives(srf, 0, u, v));
        }
    }

    for (const auto u : { 0.0f,0.1f,0.2f,0.3f,0.4f,0.5f,0.6f,0.7f,0.8f,0.9f,1.0f })
    {
        for (const auto v : { 0.0f,0.1f,0.2f,0.3f,0.4f,0.5f,0.6f,0.7f,0.8f,0.9f,1.0f })
        {
            CHECK_VERTEX_MATRIX(NURBS::SurfaceDerivatives(srf, 1, u, v), tinynurbs::surfaceDerivatives(srf, 1, u, v));
        }
    }
    
    for (const auto u : { 0.0f,0.1f,0.2f,0.3f,0.4f,0.5f,0.6f,0.7f,0.8f,0.9f,1.0f })
    {
        for (const auto v : { 0.0f,0.1f,0.2f,0.3f,0.4f,0.5f,0.6f,0.7f,0.8f,0.9f,1.0f })
        {
            CHECK_VERTEX_MATRIX(NURBS::SurfaceDerivatives(srf, 2, u, v), tinynurbs::surfaceDerivatives(srf, 2, u, v));
        }
    }
    
    for (const auto u : { 0.0f,0.1f,0.2f,0.3f,0.4f,0.5f,0.6f,0.7f,0.8f,0.9f,1.0f })
    {
        for (const auto v : { 0.0f,0.1f,0.2f,0.3f,0.4f,0.5f,0.6f,0.7f,0.8f,0.9f,1.0f })
        {
            CHECK_VERTEX_MATRIX(NURBS::SurfaceDerivatives(srf, 3, u, v), tinynurbs::surfaceDerivatives(srf, 3, u, v));
        }
    }
    
    for (const auto u : { 0.0f,0.1f,0.2f,0.3f,0.4f,0.5f,0.6f,0.7f,0.8f,0.9f,1.0f })
    {
        for (const auto v : { 0.0f,0.1f,0.2f,0.3f,0.4f,0.5f,0.6f,0.7f,0.8f,0.9f,1.0f })
        {
            CHECK_VERTEX_MATRIX(NURBS::SurfaceDerivatives(srf, 4, u, v), tinynurbs::surfaceDerivatives(srf, 4, u, v));
        }
    }
    
    for (const auto u : { 0.0f,0.1f,0.2f,0.3f,0.4f,0.5f,0.6f,0.7f,0.8f,0.9f,1.0f })
    {
        for (const auto v : { 0.0f,0.1f,0.2f,0.3f,0.4f,0.5f,0.6f,0.7f,0.8f,0.9f,1.0f })
        {
            CHECK_VERTEX_MATRIX(NURBS::SurfaceDerivatives(srf, 5, u, v), tinynurbs::surfaceDerivatives(srf, 5, u, v));
        }
    }
    
    for (const auto u : { 0.0f,0.1f,0.2f,0.3f,0.4f,0.5f,0.6f,0.7f,0.8f,0.9f,1.0f })
    {
        for (const auto v : { 0.0f,0.1f,0.2f,0.3f,0.4f,0.5f,0.6f,0.7f,0.8f,0.9f,1.0f })
        {
            CHECK_VERTEX_MATRIX(NURBS::SurfaceDerivatives(srf, 6, u, v), tinynurbs::surfaceDerivatives(srf, 6, u, v));
        }
    }
    
    for (const auto u : { 0.0f,0.1f,0.2f,0.3f,0.4f,0.5f,0.6f,0.7f,0.8f,0.9f,1.0f })
    {
        for (const auto v : { 0.0f,0.1f,0.2f,0.3f,0.4f,0.5f,0.6f,0.7f,0.8f,0.9f,1.0f })
        {
            CHECK_VERTEX_MATRIX(NURBS::SurfaceDerivatives(srf, 7, u, v), tinynurbs::surfaceDerivatives(srf, 7, u, v));
        }
    }
}
