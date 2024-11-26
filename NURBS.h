/**
  ******************************************************************************
  * @file           : Nurbs.h
  * @author         : AliceRemake
  * @brief          : None
  * @attention      : None
  * @date           : 24-11-21
  ******************************************************************************
  */



#ifndef NURBS_H
#define NURBS_H

#include <bits/stdc++.h>
#include <tinynurbs/tinynurbs.h>

namespace NURBS
{

/// @brief Binary Search For The First u_i s.t. u < u_i In Interval [degree + 1, knots.size() - degree - 1].
/// Then span Should Be i - 1. Different From A2.1 In The NURBS Book.
inline size_t FindSpan(const size_t degree, const std::vector<float>& knots, const float u) noexcept
{
    assert(!knots.empty() && knots.front() - std::numeric_limits<float>::epsilon() <= u && u <= knots.back() + std::numeric_limits<float>::epsilon());
    return (size_t)(std::upper_bound(knots.begin() + (long long)degree + 1, knots.end() - (long long)degree - 1, u) - knots.begin() - 1);
}

/// @brief Compute Nonzero B-Spline Basis Functions.
///
///    0         1            d     <--Index In b_spline_basis
/// [span  ]                        degree=0
/// [span-1] [span    ]             degree=1
/// [span-d] [span-d+1] ... [span]  degree=d
///
/// LET: b_spline_basis[i] = N_{span-d+i,d}.
/// LET: left[i] = u - u_{span+1-i}.
/// LET: right[i] = u_{span+i} - u
/// 
///              u - u_i                   u_{i+p+1} - u
/// N_{i,p} = ------------- N_{i,p-1} + ------------------- N_{i+1,p-1}
///           u_{i+p} - u_i             u_{i+p+1} - u_{i+1}
///           
///                  u - u_{i+1}                     u_{i+p+2} - u
/// N_{i+1,p} = ------------------- N_{i+1,p-1} + ------------------- N_{i+2,p-1}
///             u_{i+p+1} - u_{i+1}               u_{i+p+2} - u_{i+2}
///
///                            N_{i+1,p-1}
/// FOR SAME `p`. THE TERM ------------------- CAN BE REUSE.
///                        u_{i+p+1} - u_{i+1}
///
inline std::vector<float> BSplineBasis(const size_t degree, const size_t span, const std::vector<float>& knots, const float u) noexcept
{
    std::vector<float> b_spline_basis(degree + 1);
    std::vector<float> left(degree + 1);
    std::vector<float> right(degree + 1);

    b_spline_basis[0] = 1.0f; // N_{span,0} = 1.0f.
    
    for (size_t d = 1; d <= degree; ++d) // Loop Over Degree From `1` To `degree`.
    {
        // First Calc left And Right.
        left[d] = u - knots[span+1-d];
        right[d] = knots[span+d] - u;
        
        // The First Term Of The First Nonzero B-Spline Basis Function Is 0.0f.
        float first_term = 0.0f, reused_term = 0.0f;

        for (size_t i = 0; i < d; ++i)
        {
            reused_term = b_spline_basis[i] / (left[d-i] + right[i+1]);
            b_spline_basis[i] = first_term + right[i+1] * reused_term;
            first_term = left[d-i] * reused_term;
        }

        // The Second Term Of The Last Nonzero B-Spline Basis Function Is 0.0f.
        b_spline_basis[d] = first_term;
    }

    return b_spline_basis;
}

/// @brief Compute Nonzero Derivatives Of B-Spline Basis Functions.
///
/// LET: ndu[i][j] = N_{span+i-j,j}.                i <= j.
/// LET: ndu[i][j] = u_{span+1+j} - u_{span+1-i+j}. i >  j.
/// 
inline std::vector<std::vector<float>> BSplineDerBasis(const size_t degree, const size_t span, const std::vector<float>& knots, const float u, const size_t num_ders)
{
    std::vector ndu(degree + 1, std::vector<float>(degree + 1));
    std::vector<float> left(degree + 1);
    std::vector<float> right(degree + 1);

    ndu[0][0] = 1.0f; // N_{span,0} = 1.0f.
    
    for (size_t d = 1; d <= degree; ++d) // Loop Over Degree From `1` To `degree`.
    {
        // First Calc left And Right.
        left[d] = u - knots[span+1-d];
        right[d] = knots[span+d] - u;
        
        // The First Term Of The First Nonzero B-Spline Basis Function Is 0.0f.
        float first_term = 0.0f, reused_term = 0.0f;

        for (size_t i = 0; i < d; ++i)
        {
            ndu[d][i] = right[i+1] + left[d-i];
            reused_term = ndu[i][d-1] / (left[d-i] + right[i+1]);
            ndu[i][d] = first_term + right[i+1] * reused_term;
            first_term = left[d-i] * reused_term;
        }

        // The Second Term Of The Last Nonzero B-Spline Basis Function Is 0.0f.
        ndu[d][d] = first_term;
    }

    std::vector b_spline_der_basis(num_ders + 1, std::vector(degree + 1, 0.0f));

    // 0-th Derivatives Is Basis Functions.
    for (size_t i = 0; i <= degree; ++i)
    {
        b_spline_der_basis[0][i] = ndu[i][degree];
    }

    for (size_t d = 0; d <= degree; ++d) // Loop Over Degree.
    {
        int s1 = 0;
        int s2 = 1;
        std::vector a(2, std::vector<float>(num_ders + 1));

        a[s1][0] = 1.0f; // a_{k-1,0} = a_{0,0} = 1.0f;
        
        for (size_t k = 1; k <= std::min(degree, num_ders); ++k) // Loop Over k-th Derivative.
        {
            float derivative = 0.0f;
            
            // Calc a_{k, 0}
            if (d >= k)
            {
                a[s2][0] = a[s1][0] / ndu[degree-k+1][d-k];
                derivative = a[s2][0] * ndu[d-k][degree-k];
            }
            // Calc a_{k,j}
            for (size_t j = std::max(1ull, k-d); j <= std::min(k-1, degree-d); ++j)
            {
                a[s2][j] = (a[s1][j] - a[s1][j-1]) / ndu[degree-k+1][d-k+j];
                derivative += a[s2][j] * ndu[d-k+j][degree-k];
            }
            // Calc a_{k,k}
            if (d <= degree - k)
            {
                a[s2][k] = -a[s1][k-1] / ndu[degree-k+1][d];
                derivative += a[s2][k] * ndu[d][degree-k];
            }
            b_spline_der_basis[k][d] = derivative;
            // Flip.
            std::swap(s1, s2);
        }
    }

    for (size_t k = 1, factor = degree; k <= num_ders; factor *= (degree - k), ++k)
    {
        // degree * (degree-1) * ... * (degree-k+1)
        for (size_t d = 0; d <= degree; ++d)
        {
            b_spline_der_basis[k][d] *= (float)factor;
        }
    }
    
    return b_spline_der_basis;
}

inline glm::vec3 CurvePoint(const tinynurbs::RationalCurve<float>& crv, const float u)
{
    const size_t span = FindSpan(crv.degree, crv.knots, u);

    const auto b_spline_basis = BSplineBasis(crv.degree, span, crv.knots, u);

    std::vector<glm::vec4> homo_control_points(crv.control_points.size());

    for (size_t i = 0; i < crv.control_points.size(); ++i)
    {
        homo_control_points[i] = glm::vec4(crv.control_points[i] * crv.weights[i], crv.weights[i]);
    }

    glm::vec4 point(0.0f);

    for (size_t i = 0; i < b_spline_basis.size(); ++i)
    {
        point += b_spline_basis[i] * homo_control_points[span-crv.degree+i];
    }

    return glm::vec3(point) / point.w;
}

inline glm::vec3 SurfacePoint(const tinynurbs::RationalSurface<float>& srf, const float u, const float v)
{
    const size_t u_span = FindSpan(srf.degree_u, srf.knots_u, u);    
    const size_t v_span = FindSpan(srf.degree_v, srf.knots_v, v);

    const auto u_b_spline_basis = BSplineBasis(srf.degree_u, u_span, srf.knots_u, u);
    const auto v_b_spline_basis = BSplineBasis(srf.degree_v, v_span, srf.knots_v, v);

    std::vector homo_control_points(srf.control_points.rows(), std::vector<glm::vec4>(srf.control_points.cols()));

    for (size_t i = 0; i < srf.control_points.rows(); ++i)
    {
        for (size_t j = 0; j < srf.control_points.cols(); ++j)
        {
            homo_control_points[i][j] = glm::vec4(srf.control_points(i, j) * srf.weights(i, j), srf.weights(i, j));
        }
    }

    glm::vec4 point(0.0f);

    for (size_t i = 0; i < v_b_spline_basis.size(); ++i)
    {
        glm::vec4 tmp(0.0f);
        for (size_t j = 0; j < u_b_spline_basis.size(); ++j)
        {
            tmp += u_b_spline_basis[j] * homo_control_points[u_span-srf.degree_u+j][v_span-srf.degree_v+i];
        }
        point += v_b_spline_basis[i] * tmp;
    }

    return glm::vec3(point) / point.w;
}

// C_n^i
inline size_t Binomial(const size_t i, const size_t n)
{
    size_t ans = 1;
    if (i > n)
    {
        return 0;
    }
    for (size_t j = 1; j <= i; ++j)
    {
        ans *= n + 1 - j;
        ans /= j;
    }
    return ans;
}

inline std::vector<glm::vec3> CurveDerivatives(const tinynurbs::RationalCurve<float>& crv, const size_t num_ders, const float u)
{
    const size_t span = FindSpan(crv.degree, crv.knots, u);

    std::vector<glm::vec4> homo_control_points(crv.control_points.size());

    for (size_t i = 0; i < crv.control_points.size(); ++i)
    {
        homo_control_points[i] = glm::vec4(crv.control_points[i] * crv.weights[i], crv.weights[i]);
    }
    
    const auto b_spline_der_basis = BSplineDerBasis(crv.degree, span, crv.knots, u, num_ders);
    
    std::vector homo_curve_derivative(num_ders + 1, glm::vec4(0.0f));
    
    const size_t du = std::min(num_ders, (size_t)crv.degree);

    for (size_t k = 0; k <= du; ++k)
    {
        for (size_t j = 0; j <= crv.degree; ++j)
        {
            homo_curve_derivative[k] += b_spline_der_basis[k][j] *  homo_control_points[span-crv.degree+j];         
        }
    }

    std::vector<glm::vec3> Aders(homo_curve_derivative.size());
    std::vector<float> wders(homo_curve_derivative.size());

    for (size_t i = 0; i < homo_curve_derivative.size(); ++i)
    {
        Aders[i] = homo_curve_derivative[i];
        wders[i] = homo_curve_derivative[i].w;
    }

    std::vector<glm::vec3> ders(num_ders+1);

    for (size_t d = 0; d <= num_ders; ++d)
    {
        ders[d] = Aders[d];
        for (size_t i = 1; i <= d; ++i)
        {
            ders[d] -= Binomial(i, d) * wders[i] * ders[d-i];
        }
        ders[d] /= wders[0];
    }

    return ders;
}

inline std::vector<std::vector<glm::vec3>> SurfaceDerivatives(const tinynurbs::RationalSurface<float>& srf, const size_t num_ders, const float u, const float v)
{
    const size_t u_span = FindSpan(srf.degree_u, srf.knots_u, u);
    const size_t v_span = FindSpan(srf.degree_v, srf.knots_v, v);

    std::vector homo_control_points(srf.control_points.rows(), std::vector<glm::vec4>(srf.control_points.cols()));
    tinynurbs::array2<glm::vec4> hhomo_control_points(srf.control_points.rows(), srf.control_points.cols());

    for (size_t i = 0; i < srf.control_points.rows(); ++i)
    {
        for (size_t j = 0; j < srf.control_points.cols(); ++j)
        {
            homo_control_points[i][j] = glm::vec4(srf.control_points(i, j) * srf.weights(i, j), srf.weights(i, j));
            hhomo_control_points(i, j) = homo_control_points[i][j];
        }
    }

    auto u_b_spline_der_basis = BSplineDerBasis(srf.degree_u, u_span, srf.knots_u, u, num_ders);
    auto v_b_spline_der_basis = BSplineDerBasis(srf.degree_v, v_span, srf.knots_v, v, num_ders);

    const size_t du = std::min(num_ders, (size_t)srf.degree_u);
    const size_t dv = std::min(num_ders, (size_t)srf.degree_v);

    std::vector homo_surface_derivatives(num_ders + 1, std::vector(num_ders + 1, glm::vec4(0.0f)));
    
    for (size_t k = 0; k <= du; ++k)
    {
        std::vector<glm::vec4> temp(srf.degree_v + 1);
        for (size_t s = 0; s <= srf.degree_v; ++s)
        {
            temp[s] = glm::vec4(0.0);
            for (size_t r = 0; r <= srf.degree_u; ++r)
            {
                temp[s] += u_b_spline_der_basis[k][r] * homo_control_points[u_span+r-srf.degree_u][v_span+s-srf.degree_v];
            }
        }
        const size_t dd = std::min(num_ders-k, dv);
        for (size_t l = 0; l <= dd; ++l)
        {
            homo_surface_derivatives[k][l] = glm::vec4(0.0);
            for (size_t s = 0; s <= srf.degree_v; ++s)
            {
                homo_surface_derivatives[k][l] += v_b_spline_der_basis[l][s] * temp[s];
            }
        }
    }

    auto hhomo_surface_derivatives = tinynurbs::internal::surfaceDerivatives(srf.degree_u, srf.degree_v, srf.knots_u, srf.knots_v, hhomo_control_points, num_ders, u, v);
    
    for (size_t i = 0; i < homo_surface_derivatives.size(); ++i)
    {
        for (size_t j = 0; j < homo_surface_derivatives[i].size(); ++j)
        {
            assert(glm::distance(homo_surface_derivatives[i][j], hhomo_surface_derivatives(i, j)) < 2 * std::numeric_limits<float>::epsilon());
        }
    }

    std::vector Aders(homo_surface_derivatives.size(), std::vector<glm::vec3>(homo_surface_derivatives[0].size()));
    std::vector wders(homo_surface_derivatives.size(), std::vector<float>(homo_surface_derivatives[0].size()));

    for (size_t i = 0; i < homo_surface_derivatives.size(); ++i)
    {
        for (size_t j = 0; j < homo_surface_derivatives[i].size(); ++j)
        {
            Aders[i][j] = homo_surface_derivatives[i][j];
            wders[i][j] = homo_surface_derivatives[i][j].w;
        }
    }

    std::vector ders(num_ders + 1, std::vector(num_ders + 1, glm::vec3(0.0f)));

    for (int k = 0; k <= num_ders; ++k)
    {
        for (int l = 0; l <= num_ders - k; ++l)
        {
            auto v0 = Aders[k][l];

            for (int j = 1; j <= l; ++j)
            {
                v0 -= (float)Binomial(j, l) * wders[0][j] * ders[k][l - j];
            }

            for (int i = 1; i <= k; ++i)
            {
                v0 -= (float)Binomial(i, k) * wders[i][0] * ders[k - i][l];

                glm::vec3 v1(0.0f);
                for (int j = 1; j <= l; ++j)
                {
                    v1 -= (float)Binomial(j, l) * wders[i][j] * ders[k - 1][l - j];
                }

                v0 -= (float)Binomial(i, k) * v1;
            }

            v0 *= 1 / wders[0][0];
            ders[k][l] = v0;
        }
    }

    return ders;
}

inline glm::vec3 SurfaceNormal(const tinynurbs::RationalSurface<float>& srf, const float u, const float v)
{
    const auto surface_derivatives = SurfaceDerivatives(srf, 1, u, v);
    const auto n = glm::cross(surface_derivatives[0][1], surface_derivatives[1][0]);
    if (glm::length(n) <= std::numeric_limits<float>::epsilon())
    {
        return glm::vec3(0.0f);
    }
    return glm::normalize(n);
}

}

#endif //NURBS_H
