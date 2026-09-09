#ifndef PIPELINE_H
#define PIPELINE_H

#include "Utils.hpp"

/**
    * @file Pipeline.hpp
    * @author      Gabriel Alejandro Perez Yanez
    * @date        2026-09-08
    * @brief Data structures for handling the Graphics Pipeline.
    *
**/

/** 
    * @class pipeLine
    * @brief Manages the order of 3D transformations (Pipeline) and the camera system.
    * 
    * This class manually constructs the Projection (Perspective), View (LookAt),
    * and Model (rotation and scale) matrices. It uses quaternions to manipulate the 
    * orientation of the 3D differential volume element free of Gimbal Lock, and implements 
    * Spherical Linear Interpolation (SLERP) to smoothly animate the transition from 
    * any orientation to the principal stress basis (eigenvectors).
    * 
**/
class pipeLine {
    private:
        bool solving = false;           ///< Indicates if a SLERP animation towards the principal axes is currently running.
        float time   = 0.0f;            ///< Accumulator for the interpolation time 't' in the range [0.0, 1.0].
        float scale;                    ///< Uniform scale factor applied to the 3D differential element.

        Eigen::Matrix4f projection;     ///< Perspective projection matrix.
        Eigen::Matrix4f projView;       ///< Combined Projection * View matrix (precomputed for optimization).
        Eigen::Matrix4f view;           ///< View matrix (camera in 3D space).

        Eigen::Matrix4f model;          ///< Model matrix (Scale and Rotation of the object).
        Eigen::Matrix4f mvp;            ///< Final Model-View-Projection (MVP) matrix.

        Eigen::Quaternionf quat;        ///< Current orientation quaternion of the object.
        Eigen::Quaternionf quatStart;   ///< Source quaternion before starting the SLERP interpolation.
        Eigen::Quaternionf quatSolved;  ///< Target quaternion derived from the eigenvector matrix.

    public:

        /// @brief Gets the final combined MVP (Model-View-Projection) matrix.
        const Eigen::Matrix4f& getMVP() const { return mvp; }

        /// @brief Gets the precalculated Projection * View matrix.
        const Eigen::Matrix4f& getprojView() const { return projView; }

        /// @brief Gets the current orientation quaternion.
        const Eigen::Quaternionf& getQuat () const{return quat;}

        /// @brief Gets the current scale of the object.
        const float& getScale() const { return scale; }

        /**
            * @brief Prepares the rotation animation towards a new eigenvector basis.
            * 
            * Normalizes the columns of the rotation matrix to prevent scale distortion,
            * extracts the target quaternion, and triggers the SLERP interpolation flag.
            * 
            * @param[in] mat 3x3 Rotation matrix (usually the eigenvector matrix V).
        **/
        void setrotMat (Eigen::Matrix3f &mat){

            Eigen::Matrix3f rot = mat;

            // Column vector normalization to ensure an orthogonal rotation matrix
            rot.col(0).normalize();
            rot.col(1).normalize();
            rot.col(2).normalize();

            quatStart = quat;
            quatSolved = Eigen::Quaternionf(rot);

            solving = true;
        }

        /**
            * @brief Applies an incremental rotation to the object based on an axis and an angle.
            * 
            * Accumulates the rotation delta via quaternion multiplication and 
            * re-normalizes to prevent the accumulation of numerical floating-point errors.
            * 
            * @param[in] ag   Base rotation angle in radians.
            * @param[in] axis 3D axis around which the rotation is performed.
            * @param[in] flt  Delta time.
        **/
        void setRot (float ag, Eigen::Vector3f axis, float flt) {
            Eigen::Quaternionf dq(Eigen::AngleAxisf( ag*flt, axis) );
            quat = quat*dq; quat.normalize();
            update();
        }

        /**
            * @brief Updates the SLERP interpolation step frame by frame.
            * 
            * Smoothly and non-linearly transitions the current orientation (`quatStart`)
            * towards alignment with the principal stresses (`quatSolved`).
            * 
            * @param[in] t Delta time (dt) elapsed since the last frame.
        **/
        void solve(float t) {
            if (solving) {
                time += t;

                if (time >= 1.0f) {
                    time = 1.0f;
                    quat = quatSolved;
                    
                    solving = false;
                    time = 0.0f;
                    quatSolved = Eigen::Quaternionf::Identity(); 
                } else {
                    // Spherical Linear Interpolation of quaternions (SLERP)
                    quat = quatStart.slerp(time, quatSolved);
                }

                quat.normalize(); 
                
                update();
            }
        }

        /**
            * @brief Recalculates the model matrix and the final MVP matrix.
            * 
            * Transforms the current quaternion into its matrix representation, applies the scale 
            * using Eigen::Affine3f, and multiplies it with the previously calculated projView matrix.
        **/
        void update () {

            Eigen::Affine3f modelAffine = Eigen::Affine3f::Identity();
            modelAffine.rotate(quat);
            modelAffine.scale(scale);
            model = modelAffine.matrix();
            mvp = projView * model;
        }

        /**
            * @brief 3D Pipeline Constructor.
            * 
            * Manually constructs the Perspective Projection matrix and the LookAt View matrix.
            * 
            * @param[in] sw   Screen Width in pixels.
            * @param[in] sh   Screen Height in pixels.
            * @param[in] angle Field of View (FOV) angle in radians.
            * @param[in] np   Near clipping plane.
            * @param[in] fp   Far clipping plane.
            * @param[in] cp   Camera Position in world coordinates.
            * @param[in] tp   Target Position the camera is looking at.
            * @param[in] scl  Initial scale of the object in world space.
        **/
        pipeLine(int sw, int sh, float angle, float np, float fp, Eigen::Vector3f cp, Eigen::Vector3f tp, float scl){
            // Manual construction of the Perspective Projection Matrix
            projection = Eigen::Matrix4f::Zero();
            float tn = std::tanf(angle*0.5f);
            float as =  static_cast<float>(sw) / static_cast<float>(sh);

            projection(0,0) = 1.0f / (tn*as);
            projection(1,1) = 1.0f / tn;
            projection(2,2) =  fp/(fp-np);
            projection(2,3) = -(fp*np)/(fp-np);
            projection(3,2) =  1.0f;
            projection(3,3) =  0.0f;

            // Manual construction of the View Matrix (LookAt System)
            Eigen::Vector3f worldUp = Eigen::Vector3f(0.0f,1.0f,0.0f);
            Eigen::Vector3f cameraForward = (tp - cp).normalized();
            Eigen::Vector3f cameraRight   = cameraForward.cross( worldUp).normalized();
            Eigen::Vector3f cameraUp      = cameraRight.cross(cameraForward).normalized();

            view = Eigen::Matrix4f::Identity();

            view.row(0).head<3>() =   cameraRight.transpose();  
            view.row(1).head<3>() =      cameraUp.transpose();  
            view.row(2).head<3>() = cameraForward.transpose();  
            view(0,3) = -cp.dot(cameraRight);
            view(1,3) = -cp.dot(cameraUp);
            view(2,3) = -cp.dot(cameraForward);

            scale = scl;
            projView = projection*view;
            quat = Eigen::Quaternionf::Identity();

            update ();
        }   
    
    /// @brief Default constructor.
    pipeLine() = default;
};

#endif