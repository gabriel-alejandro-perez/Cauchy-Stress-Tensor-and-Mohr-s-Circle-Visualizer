#ifndef CUBE_H
#define CUBE_H

#include "Pipeline.hpp"

/**
    * @file Cube.hpp
    * @author      Gabriel Alejandro Perez Yanez
    * @date        2026-09-08
    * @brief Data structures for computing stresses on the faces of the differential element 
    *        and managing the 3D cube geometry.
**/

/**
    * @struct FaceStress
    * @brief Stores the position vectors and stress components for a single face of the cube.
    * 
    * Used to decompose the traction vector into its normal (axial) component 
    * and two tangential (shear) components over the surface of a face.
**/
struct FaceStress {
    Eigen::Vector3f center;    ///< Geometric center of the face in 3D space.
    Eigen::Vector3f axialVec;  ///< Normal stress vector to the face (perpendicular).
    Eigen::Vector3f shearVec1; ///< First shear stress vector in the plane of the face.
    Eigen::Vector3f shearVec2; ///< Second shear stress vector in the plane of the face.
};

/**
    * @brief Represents the Cauchy Stress Tensor and handles principal stress resolution.
    * 
    * Contains the 3x3 stress matrix (`sigma`), calculates the principal stresses 
    * ($\sigma_1 \ge \sigma_2 \ge \sigma_3$) by sorting the eigenvalues, and computes the 
    * maximum shear stress ($\tau_{max}$). It also allows computing the traction vectors 
    * on individual faces of the cube.
    *
**/
struct tensor{

    /**
        * @brief Diagonalizes the stress tensor to find principal stresses and directions.
        * 
        * Invokes the Jacobi iterative solver (`mathUtils::eigenVV`). Saves the eigenvector 
        * matrix in `solut` before sorting to preserve spatial orientation 
        * (preventing negative determinants), and then sorts the eigenvalues in descending order:
        * \f$ \sigma_1 \ge \sigma_2 \ge \sigma_3 \f$.
        *
    **/
    void solve(){
        Eigen::Matrix3f Vec;
        Eigen::Matrix3f Mat = sigma;

        mathUtils::eigenVV( Mat, Vec );

        //@note Check before sorting or the determinant could equal -1 
        solut = Vec; 

        if (Mat(0, 0) < Mat(1, 1)) {
            std::swap(Mat(0, 0), Mat(1, 1));
            Vec.col(0).swap(Vec.col(1));
        }
        if (Mat(0, 0) < Mat(2, 2)) {
            std::swap(Mat(0, 0), Mat(2, 2));
            Vec.col(0).swap(Vec.col(2));
        }
        if (Mat(1, 1) < Mat(2, 2)) {
            std::swap(Mat(1, 1), Mat(2, 2));
            Vec.col(1).swap(Vec.col(2));
        }

        // Principal stresses
        sigma1 = Mat(0, 0);
        sigma2 = Mat(1, 1);
        sigma3 = Mat(2, 2);

        // Absolute maximum shear (Radius of the main Mohr's Circle)
        tauMax = (sigma1 - sigma3) / 2.0;
    }

    /**
        * @brief Computes the stress decomposition (traction) on a face of the differential element.
        * 
        * Given the normal vector to a face and two orthogonal tangent vectors, it computes the 
        * traction vector using \f$ \mathbf{T} = \mathbf{\Sigma} \cdot \mathbf{n} \f$ and projects it 
        * into its axial (normal) and shear (tangential) components.
        * 
        * @param[in] globalSigma Stress tensor matrix in global coordinates.
        * @param[in] normal      Unit normal vector to the face.
        * @param[in] tangent1    First unit tangent vector to the face.
        * @param[in] tangent2    Second unit tangent vector to the face (orthogonal to tangent1).
        * @param[in] faceCenter  Position of the geometric center of the face.
        * 
        * @return FaceStress struct loaded with the projected components.
    **/
    FaceStress calculatefaceStress(const Eigen::Matrix3f& globalSigma, const Eigen::Vector3f& normal, const Eigen::Vector3f& tangent1, const Eigen::Vector3f& tangent2, const Eigen::Vector3f& faceCenter) {
        FaceStress face;
        face.center = faceCenter;
        
        // Traction vector T = Sigma * n
        Eigen::Vector3f traction = globalSigma * normal;
        
        // Projection of the traction vector into its axial and shear components
        face.axialVec  = traction.dot(normal)   * normal; 
        face.shearVec1 = traction.dot(tangent1) * tangent1;
        face.shearVec2 = traction.dot(tangent2) * tangent2;
        
        return face;
    }

    Eigen::Matrix3f sigma; ///< 3x3 Matrix of the Cauchy Stress Tensor.
    Eigen::Matrix3f solut; ///< Eigenvector matrix (orientation towards the principal basis).

    float sigma1 = 0.0f; ///< Maximum principal stress (Sigma 1).
    float sigma2 = 0.0f; ///< Intermediate principal stress (Sigma 2).
    float sigma3 = 0.0f; ///< Minimum principal stress (Sigma 3).
    float tauMax = 0.0f; ///< Maximum shear stress (Tau Max).

};

/**
    * @brief Vectorial geometric representation (wireframe) of the cubic differential element.
    * 
    * Defines the local coordinates of the 8 vertices of the cube centered at the origin 
    * and the index pairs that make up its 12 edges for line rendering.
**/
struct cube {

    Eigen::Vector4f tvrts[8]; ///< Vertices transformed to screen/clip coordinates after applying MVP.

    /// @brief Original vertices of the origin-centered cube in homogeneous coordinates (x, y, z, w).
    Eigen::Vector4f verts[8] = {
        Eigen::Vector4f(-1.0f, -1.0f, -1.0f, 1.0f), //0
        Eigen::Vector4f( 1.0f, -1.0f, -1.0f, 1.0f), //1
        Eigen::Vector4f(-1.0f,  1.0f, -1.0f, 1.0f), //2
        Eigen::Vector4f( 1.0f,  1.0f, -1.0f, 1.0f), //3
        Eigen::Vector4f(-1.0f, -1.0f,  1.0f, 1.0f), //4
        Eigen::Vector4f( 1.0f, -1.0f,  1.0f, 1.0f), //5
        Eigen::Vector4f(-1.0f,  1.0f,  1.0f, 1.0f), //6
        Eigen::Vector4f( 1.0f,  1.0f,  1.0f, 1.0f)  //7
    };

    /// @brief Pairs of vertex indices that define the 12 edges of the wireframe cube.
    Eigen::Vector2i lines[12] = {
        Eigen::Vector2i(0, 2), Eigen::Vector2i(2, 3), Eigen::Vector2i(3, 1), Eigen::Vector2i(1, 0),
        Eigen::Vector2i(4, 6), Eigen::Vector2i(6, 7), Eigen::Vector2i(7, 5), Eigen::Vector2i(5, 4),
        Eigen::Vector2i(3, 7), Eigen::Vector2i(1, 5), Eigen::Vector2i(0, 4), Eigen::Vector2i(2, 6)
    };

};

#endif