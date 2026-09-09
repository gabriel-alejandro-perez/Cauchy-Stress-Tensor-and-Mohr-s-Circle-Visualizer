#ifndef UTILS_H
#define UTILS_H

#include <Eigen/Dense>
#include <vector>

/**
    * @file  Utils.hpp
    * @author      Gabriel Alejandro Perez Yanez
    * @date        2026-09-08
    * @brief Utility data structures including constants, solvers, and parsers.
    *
**/

/// @brief Pi constant (single precision).
inline constexpr float pi = 3.14159265358979323846;

/**
    * @brief Layout and positioning constants for rendering the matrix on screen.
    * 
    * Defines the dimensions, font scaling, and origin coordinates used 
    * by the graphics engine to draw the visual representation of the Cauchy 
    * Stress Tensor (3x3 matrix) and its brackets in the user interface.
    * 
**/
namespace matrixLayout {

    constexpr int scale      = 2;                ///< Scale factor for text and UI elements.
    constexpr int widthChar  = 8 * scale;        ///< Width in pixels of an individually rendered character (8px base * scale).
    constexpr int baseX      = 390;              ///< Origin X coordinate (top-left corner) for drawing the matrix.
    constexpr int baseY      = 60;               ///< Origin Y coordinate for drawing the matrix.
    constexpr int heightRow  = 25;               ///< Vertical distance (height) in pixels between each row of the matrix.
    constexpr int marginCol  = 20;               ///< Horizontal separation margin in pixels between columns.
    constexpr int width      = (7 * widthChar);  ///< Width in pixels reserved for the value string (max 7 characters).
    constexpr int widthCol   = width + 20;       ///< Total column width including the formatted value and its margin.
    constexpr int heightMat  = 70;               ///< Total height in pixels of the matrix's visual brackets.
    constexpr int widthMat   = 5;                ///< Thickness/Width in pixels of the brackets delimiting the matrix.
    constexpr int filCol     = 3;                ///< Matrix dimension (3 rows x 3 columns).
}

namespace mathUtils {   
    /**
        * @brief Calculates the eigenvalues and eigenvectors of a 3x3 symmetric matrix in real-time.
        * 
        * Implements an iterative solver based on the Jacobi algorithm. It is highly
        * optimized for the specific case of 3x3 matrices (such as the Cauchy Tensor), 
        * avoiding dynamic loops and using direct search for the maximum off-diagonal element.
        * 
        * To maximize performance, the angles of the Givens rotations (sine and cosine) 
        * are deduced purely algebraically from the matrix entries, completely bypassing 
        * expensive trigonometric functions (atan, sin, cos).
        * 
        * @param[in,out] A 3x3 Eigen symmetric matrix. Upon completion, this matrix is 
        *                  modified in-place (A = P^T * A * P) and will contain the principal 
        *                  stresses (eigenvalues) on its diagonal.
        * @param[out] V    3x3 Eigen matrix initialized internally as identity. 
        *                  Upon completion, its columns will contain the eigenvectors 
        *                  (principal directions).
        * @return true if the algorithm converges before reaching the maximum iteration limit.
        * @return false if the limit (maxit) is reached without converging to the tolerance threshold.
        *
        * @note The theoretical basis and algorithmic optimization techniques for this 
        * implementation were adapted from Numerical Recipes.
        *
        * @see Press, W. H., Teukolsky, S. A., Vetterling, W. T., & Flannery, B. P. (2007). 
        * "Numerical Recipes: The Art of Scientific Computing" (3rd ed.). 
        * Cambridge University Press. (Chapter 11: Eigensystems).
        *
    **/
    bool eigenVV( Eigen::Matrix3f& A, Eigen::Matrix3f& V )
    {   
        assert(A.isApprox(A.transpose()) && "A must be symmetric");

        float norm = A.norm();
        V = Eigen::Matrix3f::Identity();

        /// @note If the matrix is null, it is already diagonalized.
        if(norm == 0) return true;

        /// @brief Dynamic sweep threshold based on machine epsilon and matrix norm.
        float toler = std::numeric_limits<float>::epsilon()*norm;

        /// @brief Maximum number of iterations.
        int maxit  = 100;

        float tau  = 0.0f;
        float c    = 0.0f; // Cosine of the angle
        float s    = 0.0f; // Sine of the angle
        float t    = 0.0f; // Tangent of the angle

        int p, q, r;

        for(int x=0;x<maxit;x++)
        {

            // Hardcoded search for the maximum off-diagonal value
            // Indices of rotated axes (p,q) and the invariant axis (r)
            p = 0; q = 1; r = 2;
            float max_offdiag = std::abs(A(1,0));

            if( std::abs(A(2,0)) > max_offdiag ) 
            {
                max_offdiag = std::abs(A(2,0));
                p = 0; q = 2; r = 1;
            }
            if(std::abs(A(2,1)) > max_offdiag)
            {
                max_offdiag = std::abs(A(2,1));
                p = 1; q = 2; r = 0;
            }

            // Exit condition: Off-diagonal elements are negligible
            if (max_offdiag < toler) return true;

            // Algebraic deduction of the rotation without trigonometric functions
            // Givens rotation (without math::sin/cos/atan)
            tau = (A(p,p)-A(q,q)) / (2.0f*A(p,q));
            t = std::copysign(1.0f,tau) / (std::abs(tau) + std::sqrtf((tau*tau) +1));
            c = 1.0f / std::sqrtf((t*t) + 1.0f); // cos(θ) = 1 / √ (t^2 + 1)
            s = t*c;                             // tan(θ) = sin(θ) / cos(θ)

            // Temporary variables for the update
            float pp = A(p,p) + A(p,q) * t;
            float qq = A(q,q) - A(p,q) * t;

            float pr =  A(p,r)*c  + s * A(q,r); 
            float qr =  A(p,r)*-s + c * A(q,r);  

            /// @brief Applies similarity transformation: A = P^T * A * P
            A(p,p) = pp;
            A(q,q) = qq;
            
            A(p,r) = pr;
            A(q,r) = qr;

            A(r,q) = qr;
            A(r,p) = pr;

            A(p,q) = 0.0f; 
            A(q,p) = 0.0f;

            /// @brief Accumulates the rotation in the eigenvector matrix: V = V * P
            for(int y = 0; y < 3; y ++)
            {
                float yp =  V(y,p)*c + V(y,q)*s;
                float yq = -V(y,p)*s + V(y,q)*c;
                
                V(y,p) = yp;
                V(y,q) = yq;
            }
        }
        return false;
    }
}

/**
    * @brief Safely converts a text string to a floating-point number (float).
    *
    * Designed to process user inputs intended for the Cauchy Stress Tensor matrix. 
    * Automatically handles edge cases and standard C++ exceptions (empty strings, incomplete 
    * signs, or out-of-range values) to prevent the application from crashing. Additionally, 
    * normalizes the value '-0.0f' to '0.0f' to avoid display issues.
    *
    * @param[in] str Text string (std::string) entered by the user in the interface.
    *
    * @return The numeric value converted to float if the format is valid.
    * @return 0.0f if an invalid input is detected, is empty, or throws an exception.
    *
 */
inline float tryStof(const std::string& str) {
    if (str.empty() || str == "-" || str == "." || str == "-.") return 0.0f;
    try {
        float n = std::stof(str);
        return (n == -0.0f) ? 0.0f : n;
    } catch (const std::invalid_argument&) {
        return 0.0f;
    } catch (const std::out_of_range&) {
        return 0.0f; 
    }
}

/**
    * @brief Converts a floating-point number to a formatted string with exactly two decimal places.
    *
    * Uses a string stream (ostringstream) alongside format manipulators 
    * (std::fixed and std::setprecision) to ensure uniform representation. 
    * This is essential for rendering clean numerical data in the interface (such as 
    * principal stresses, tensor components, and Mohr's Circle data).
    *
    * @param[in] val Numeric float value to format.
    *
    * @return Text string (std::string) with the value represented to two decimal places (e.g., 15.00).
    * 
 **/
inline std::string Format2Dec(float val) {
    std::ostringstream stream;
    stream << std::fixed << std::setprecision(2) << val;
    return stream.str();
}

/**
    * @brief Precalculated mapping structure for the CRT monitor visual filter.
    * 
    * Stores screen coordinates transformed by barrel distortion (Brown-Conrady) 
    * along with their corresponding vignette attenuation factor. Allows 
    * applying the CRT effect via a Lookup Table (LUT), avoiding expensive 
    * geometric function calculations per pixel on every frame.
    * 
 **/
struct crtMap 
{
    int   ux; ///< Remapped X coordinate in the source buffer.
    int   uy; ///< Remapped Y coordinate in the source buffer.
    float vm; ///< Vignette attenuation factor for the monitor edges.
};

// Windows Compilation //
//g++ -o Mohr main.cpp -luser32 -lgdi32 -lopengl32 -lgdiplus -lShlwapi -ldwmapi -lstdc++fs -static -std=c++17 -I "C:\msys64\mingw64\include\eigen3" -O3//

#endif