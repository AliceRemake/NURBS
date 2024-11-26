/**
  ******************************************************************************
  * @file           : TestSurfaceNormal.cpp
  * @author         : AliceRemake
  * @brief          : None
  * @attention      : None
  * @date           : 24-11-26
  ******************************************************************************
  */



#include <doctest/doctest.h>
#include <NURBS.h>
#include <tinynurbs/tinynurbs.h>

#define CHECK_GLM_VERTEX(lhs, rhs) CHECK((glm::distance(lhs, rhs) < 2 * std::numeric_limits<float>::epsilon()))

TEST_CASE("SurfaceNormal")
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
            CHECK_GLM_VERTEX(NURBS::SurfaceNormal(srf, u, v), tinynurbs::surfaceNormal(srf, u, v));
        }
    }
}
