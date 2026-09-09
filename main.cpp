#define OLC_PGE_APPLICATION
#include "olcPixelGameEngine.h"

#include "Cube.hpp"

/**
    * @file main.cpp
    * @author      Gabriel Alejandro Perez Yanez
    * @date        2026-09-08
    * @brief Entry file of the program where all the 3D rendering and UI logic is located.
**/

/**
    * @class visualizer3D
    * @brief Interactive application for the visualization of Mohr's Circle and the 3D stress state.
    * 
    * Extends `olc::PixelGameEngine` to render a graphical interface that includes:
    * - A matrix editor for the Cauchy Stress Tensor.
    * - A 3D representation of the differential element subjected to stresses.
    * - A 2D diagram of Mohr's Circle.
    * - A retro CRT post-processing effect (barrel distortion).
    *
**/
class visualizer3D : public olc::PixelGameEngine {
private:
    std::vector<std::vector<std::string>> sigmaS;   ///< Text representation of the stress matrix for the interface.
    std::unique_ptr<olc::Sprite> canvaSprite;       ///< Intermediate canvas to render the scene before software post-processing.
    std::vector<crtMap> crtLut;                     ///< Look-Up Table (LUT) for CRT lens distortion mapping.

    int   screenHeight;    int screenWidth;         ///< Full screen dimensions.
    int   hscreenHeight;   int hscreenWidth;        ///< Half of the screen dimensions (center).
    int   syPos;           int sxPos;               ///< Base position (offset) for 3D projection on screen.
    float inscreenWidth2; float inscreenHeight2;    ///< Precalculated multipliers for coordinate normalization.

    float k = 0.15f;                                ///< Barrel distortion coefficient for the CRT effect.

    cube differentialCube;                          ///< Geometric structure of the cubic differential element.
    pipeLine rendpipeLine;                          ///< Camera system and matrix transformation (MVP).

    tensor tensorState;                             ///< Mathematical engine that stores the tensor and solves the principal stresses.

    int cellSelX = 0;                               ///< X index of the currently selected stress matrix cell.
    int cellSelY = 0;                               ///< Y index of the currently selected stress matrix cell.
    
    /**
     * @brief Synchronizes the user's text input with the mathematical tensor.
     * 
     * Converts the modified text of the cell to float, forces tensor symmetry
     * ($\sigma_{ij} = \sigma_{ji}$) and recalculates the eigenvalues and eigenvectors.
     * Solves the stress tensor.
     *
    **/
    void hasChange(){
        char buffer[16];
        std::string& cellActs = sigmaS[cellSelY][cellSelX];
        
        tensorState.sigma(cellSelY,cellSelX) = tryStof(cellActs);
        std::snprintf(buffer, sizeof(buffer), "%7.2f", tensorState.sigma(cellSelY,cellSelX));
        sigmaS[cellSelY][cellSelX] = buffer;

        // Maintain the symmetry of the stress tensor
        if(cellSelY != cellSelX){ sigmaS[cellSelX][cellSelY] = buffer; tensorState.sigma(cellSelX,cellSelY) = tryStof(cellActs);}

        tensorState.solve();
    }

    /**
        * @brief Draws a 3D stress vector projected on the screen.
        * 
        * @param[in] mp Transformation matrix.
        * @param[in] startPos 3D start coordinate of the vector (center of the face).
        * @param[in] forceVec Force/stress vector to represent.
        * @param[in] color Base color of the vector.
        * @param[in] label Optional text label to identify the vector (e.g., "G11").
        *
    **/
    void DrawStressVector(const Eigen::Matrix4f& mp, const Eigen::Vector3f& startPos, const Eigen::Vector3f& forceVec, olc::Pixel color, const std::string& label = "") {
        if (forceVec.squaredNorm() < 1e-4f) return; 

        // Darkening based on depth (Z Axis)
        float zNorm = (-startPos.z() + 1.0f) * 0.5f; 
        float depthFactor = 0.30f + (0.70f * zNorm); 

        olc::Pixel dimmedColor(
            static_cast<uint8_t>(color.r * depthFactor),
            static_cast<uint8_t>(color.g * depthFactor),
            static_cast<uint8_t>(color.b * depthFactor)
        );

        float vectorScale = 0.005f; 
        float maxLength   =  0.95f; 

        Eigen::Vector3f displacement = forceVec * vectorScale;

        // Limit the maximum length of the vector in the visualization
        if (displacement.squaredNorm() > (maxLength * maxLength)) {
            displacement = displacement.normalized() * maxLength;
        }

        Eigen::Vector3f endPos = startPos + displacement;

        // Projection to Clip Space
        Eigen::Vector4f p0Clip = mp * Eigen::Vector4f(startPos.x()*rendpipeLine.getScale(), startPos.y()*rendpipeLine.getScale(), startPos.z()*rendpipeLine.getScale(), 1.0f);
        Eigen::Vector4f p1Clip = mp * Eigen::Vector4f(endPos.x()*rendpipeLine.getScale(), endPos.y()*rendpipeLine.getScale(), endPos.z()*rendpipeLine.getScale(), 1.0f);

        if (p0Clip.w() < 0.1f || p1Clip.w() < 0.1f) return;

        // Perspective Division
        p0Clip.head<3>() /= p0Clip.w();
        p1Clip.head<3>() /= p1Clip.w();

        // Mapping to Screen Coordinates
        int x0 = static_cast<int>((p0Clip.x() * hscreenWidth) + sxPos);
        int y0 = static_cast<int>((-p0Clip.y() * hscreenHeight) + syPos);
        int x1 = static_cast<int>((p1Clip.x() * hscreenWidth) + sxPos);
        int y1 = static_cast<int>((-p1Clip.y() * hscreenHeight) + syPos);

        // Vector rendering (thicker line and end point)
        if (std::abs(x1 - x0) <= 2 && std::abs(y1 - y0) <= 2) {
            DrawCircle(x1, y1, 4, dimmedColor); 
        } else {
            DrawLine(x0, y0, x1, y1, dimmedColor);
            DrawLine(x0 + 1, y0, x1 + 1, y1, dimmedColor);
            DrawLine(x0, y0 + 1, x1, y1 + 1, dimmedColor);
            FillCircle(x1, y1, 2, dimmedColor); 
        }

        // Text label rendering
        if (!label.empty()) {
            float dx = static_cast<float>(x1 - x0);
            float dy = static_cast<float>(y1 - y0);
            float len = std::sqrt(dx * dx + dy * dy);
            
            if (len > 0.1f) {
                dx /= len;
                dy /= len;
            } else {
                dx = 1.0f; dy = -1.0f;
            }

            int textWidth = static_cast<int>(label.length()) * 8;
            int textHeight = 8;
 
            float pushOut = 10.0f; 
            int labelX = x1 + static_cast<int>(dx * pushOut);
            int labelY = y1 + static_cast<int>(dy * pushOut);

            if (dx < 0.0f) labelX -= textWidth;
            if (dy < 0.0f) labelY -= textHeight;
            
            // Keep text within screen boundaries
            labelX = std::max(0, std::min(labelX, screenWidth - textWidth));
            labelY = std::max(0, std::min(labelY, screenHeight - textHeight));

            DrawString(labelX, labelY, label, dimmedColor, 1);
        }
    }

    /**
        * @brief Processes keyboard inputs to interact with the application.
        * 
        * Supports matrix navigation (Arrows), number input, negatives and decimals,
        * manual cube rotation (Q, W, E, A, S, D) and alignment with principal axes (R).
        * 
        * @param[in] fElapsedTime Time elapsed since the last frame (Delta Time).
        *
    **/
    void handleInput(float fElapsedTime) {
        std::string& cellActs = sigmaS[cellSelY][cellSelX];

        // Input of numeric digits
        for (int i = 0; i <= 9; i++) {
            if (GetKey(static_cast<olc::Key>(static_cast<int>(olc::Key::K0) + i)).bPressed ||
                GetKey(static_cast<olc::Key>(static_cast<int>(olc::Key::NP0) + i)).bPressed) {
                std::string n = cellActs + std::to_string(i);
                size_t dot = n.find('.');
                bool addDeci = true;

                if (dot != std::string::npos) {
                    int len = n.length() - dot - 1;
                    if (len > 2) addDeci = false;
                }

                if (addDeci) {
                    float value = tryStof(n);
                    if (value <= 1000.0f && value >= -999.0f) {
                        cellActs += std::to_string(i);
                    }
                }
            }
        }

        // Signs, decimals and backspace
        if (GetKey(olc::Key::MINUS).bPressed || GetKey(olc::Key::NP_SUB).bPressed) {
            if (cellActs.empty()) cellActs += "-";
            else if(tryStof(cellActs) == 0.0f) cellActs = "-";
        }
        if (GetKey(olc::Key::PERIOD).bPressed || GetKey(olc::Key::NP_DECIMAL).bPressed) {
            if (cellActs.find('.') == std::string::npos) cellActs += ".";
        }
        if (GetKey(olc::Key::BACK).bPressed && !cellActs.empty()) {
            cellActs.pop_back();
            if(tryStof(cellActs) == 0.0f) cellActs = "";
        }

        // Matrix navigation
        if (GetKey(olc::Key::UP).bPressed)    { hasChange(); cellSelY = std::max(0, cellSelY - 1); }
        if (GetKey(olc::Key::DOWN).bPressed)  { hasChange(); cellSelY = std::min(matrixLayout::filCol - 1, cellSelY + 1); }
        if (GetKey(olc::Key::LEFT).bPressed)  { hasChange(); cellSelX = std::max(0, cellSelX - 1); }
        if (GetKey(olc::Key::RIGHT).bPressed) { hasChange(); cellSelX = std::min(matrixLayout::filCol - 1, cellSelX + 1); }

        // Manual rotation of the differential cube
        if (GetKey(olc::Key::Q).bHeld) rendpipeLine.setRot ( 0.5f, Eigen::Vector3f::UnitX(), fElapsedTime);
        if (GetKey(olc::Key::W).bHeld) rendpipeLine.setRot ( 0.5f, Eigen::Vector3f::UnitY(), fElapsedTime);
        if (GetKey(olc::Key::E).bHeld) rendpipeLine.setRot ( 0.5f, Eigen::Vector3f::UnitZ(), fElapsedTime);
        if (GetKey(olc::Key::A).bHeld) rendpipeLine.setRot (-0.5f, Eigen::Vector3f::UnitX(), fElapsedTime);
        if (GetKey(olc::Key::S).bHeld) rendpipeLine.setRot (-0.5f, Eigen::Vector3f::UnitY(), fElapsedTime);
        if (GetKey(olc::Key::D).bHeld) rendpipeLine.setRot (-0.5f, Eigen::Vector3f::UnitZ(), fElapsedTime);

        // Align cube with principal stresses
        if (GetKey(olc::Key::R).bPressed) rendpipeLine.setrotMat(tensorState.solut);

        // Update SLERP animation interpolation (if active)
        rendpipeLine.solve(fElapsedTime);
    }

    /**
        * @brief Applies mathematical transformations to the vertices of the 3D model.
        * 
        * @param[in] fElapsedTime Time elapsed since the last frame (Delta Time).
        *
    **/
    void updateLogic(float fElapsedTime) {
        const Eigen::Matrix4f& mvp = rendpipeLine.getMVP();

        for(int i = 0; i < 8; i++){
            Eigen::Vector4f n = mvp * differentialCube.verts[i];
            if (n(3) != 0.0f) {
                n.head<3>() /= n(3); 
            }
            differentialCube.tvrts[i] = n;
        }
    }

    /**
        * @brief Builds and draws all graphical elements on the secondary buffer (canvaSprite).
    **/
    void renderScene() {
        SetDrawTarget(canvaSprite.get());
        Clear(olc::Pixel(10, 20, 10));

        // Draw the side panel with Mohr's Circle
        DrawMohrDiagram(380, 250, 260, 200);

        // Main texts and labels
        FillRect(20, 10, 220, 40, olc::Pixel(150, 255, 150));
        DrawString(25, 22,  "MOHR'S CIRCLE", olc::BLACK, 2);
        DrawString(20, 60,  "BY GABRIEL PEREZ", olc::Pixel(100, 255, 100), 2);
        DrawString(450, 20, "CAUCHY TENSOR", olc::Pixel(100, 200, 100), 2);

        int offsetRotadoY = 140; 
        Eigen::Matrix3f T = rendpipeLine.getQuat().toRotationMatrix();
        Eigen::Matrix3f sigmaRotated = T.transpose() * tensorState.sigma * T;

        // Rendering of the Stress Tensor editor
        for (int c = 0; c < 3; ++c) {
            const int posX = matrixLayout::baseX + (c * matrixLayout::widthCol);
            for (int f = 0; f < 3; ++f) {
                DrawString(posX, matrixLayout::baseY + (f * matrixLayout::heightRow), 
                           sigmaS[f][c], olc::Pixel(100, 200, 100), matrixLayout::scale);

                char buffer[16];

                float val = std::abs(sigmaRotated(f, c)) < 1e-4f ? 0.0f : sigmaRotated(f, c);
                std::snprintf(buffer, sizeof(buffer), "%7.2f", val);
                
                DrawString(posX, (offsetRotadoY + 40) + (f * matrixLayout::heightRow), buffer, olc::Pixel(255, 255, 100), matrixLayout::scale);
            }
        }   

        DrawString(410, offsetRotadoY+15, "ROTATED TENSOR (LOCAL)", olc::Pixel(100, 200, 100), 2);

        for (int c = 0; c < 3; ++c) {
            const int posX = matrixLayout::baseX + (c * matrixLayout::widthCol);
            for (int f = 0; f < 3; ++f) {
                char buffer[16];
                
                // Avoid mathematical noise and "-0.00" from floats
                float val = std::abs(sigmaRotated(f, c)) < 1e-4f ? 0.0f : sigmaRotated(f, c);
                std::snprintf(buffer, sizeof(buffer), "%7.2f", val);
                
                DrawString(posX, (offsetRotadoY + 40) + (f * matrixLayout::heightRow), 
                           buffer, olc::Pixel(255, 255, 100), matrixLayout::scale);
            }
        }
        
        // Matrix brackets and selection indicator
        FillRect(matrixLayout::baseX - (matrixLayout::widthMat*2), matrixLayout::baseY, matrixLayout::widthMat, matrixLayout::heightMat, olc::RED);
        FillRect(matrixLayout::baseX*2, matrixLayout::baseY, matrixLayout::widthMat, matrixLayout::heightMat, olc::RED);

        FillRect(matrixLayout::baseX - (matrixLayout::widthMat*2), matrixLayout::baseY+(matrixLayout::heightMat*1.65f), matrixLayout::widthMat, matrixLayout::heightMat, olc::RED);
        FillRect(matrixLayout::baseX*2, matrixLayout::baseY+(matrixLayout::heightMat*1.65f), matrixLayout::widthMat, matrixLayout::heightMat, olc::RED);
        DrawRect(matrixLayout::baseX + (cellSelX * matrixLayout::widthCol), matrixLayout::heightRow*2 + (cellSelY * matrixLayout::heightRow), matrixLayout::width, matrixLayout::heightRow, olc::RED);

        // Depth calculation
        float minZ = std::numeric_limits<float>::max();
        float maxZ = std::numeric_limits<float>::lowest();

        for (int v = 0; v < 8; v++) { 
            float z = differentialCube.tvrts[v](2);
            minZ = std::min(minZ, z);
            maxZ = std::max(maxZ, z);
        }

        float zRange = (maxZ - minZ < 1e-5f) ? 1.0f : (maxZ - minZ); 

        // Rendering of cube edges
        for(int i = 0; i < 12; i++){
            int li0 = differentialCube.lines[i](0); 
            int li1 = differentialCube.lines[i](1);

            float x0 = differentialCube.tvrts[li0](0); float y0 = differentialCube.tvrts[li0](1);
            float x1 = differentialCube.tvrts[li1](0); float y1 = differentialCube.tvrts[li1](1);
            float z1 = differentialCube.tvrts[li0](2); float z2 = differentialCube.tvrts[li1](2);

            float midZ = (z1 + z2) * 0.5f; 

            float zNorm = (midZ - minZ) / zRange; 
            
            float depthFactor = std::clamp(0.15f + (0.85f * (1.0f - zNorm)), 0.0f, 1.0f); 

            DrawLine(
                static_cast<int>((x0 * hscreenWidth) + sxPos),
                static_cast<int>((-y0 * hscreenHeight) + syPos),
                static_cast<int>((x1 * hscreenWidth) + sxPos),
                static_cast<int>((-y1 * hscreenHeight) + syPos),
                olc::Pixel(0,static_cast<uint8_t>(255 * depthFactor),0)
            );
        }

        // Extraction of normal vectors from the object's quaternion
        const Eigen::Matrix4f& mvp = rendpipeLine.getprojView();
        
        Eigen::Matrix3f R = rendpipeLine.getQuat().toRotationMatrix();
        Eigen::Vector3f nX = R.col(0); 
        Eigen::Vector3f nY = R.col(1); 
        Eigen::Vector3f nZ = R.col(2); 

        // Calculation of traction and axial/shear stresses on each face
        FaceStress stressX   = tensorState.calculatefaceStress(tensorState.sigma, nX, nY, nZ, nX);
        FaceStress stressY   = tensorState.calculatefaceStress(tensorState.sigma, nY, nX, nZ, nY);
        FaceStress stressZ   = tensorState.calculatefaceStress(tensorState.sigma, nZ, nX, nY, nZ);
        FaceStress stress_nX = tensorState.calculatefaceStress(tensorState.sigma, -nX, nY, nZ, -nX);
        FaceStress stress_nY = tensorState.calculatefaceStress(tensorState.sigma, -nY, nX, nZ, -nY);
        FaceStress stress_nZ = tensorState.calculatefaceStress(tensorState.sigma, -nZ, nX, nY, -nZ);

        // Drawing of axial components (Normals)
        DrawStressVector(mvp, stressX.center, stressX.axialVec, olc::GREEN, "G11'");
        DrawStressVector(mvp, stressY.center, stressY.axialVec, olc::GREEN, "G22");
        DrawStressVector(mvp, stressZ.center, stressZ.axialVec, olc::GREEN, "G33'");

        DrawStressVector(mvp, stress_nX.center, stress_nX.axialVec, olc::GREEN, "G11");
        DrawStressVector(mvp, stress_nY.center, stress_nY.axialVec, olc::GREEN, "G22'");
        DrawStressVector(mvp, stress_nZ.center, stress_nZ.axialVec, olc::GREEN, "G33");

        // Drawing of shear components
        DrawStressVector(mvp, stressX.center, stressX.shearVec1, olc::GREEN, "G12");
        DrawStressVector(mvp, stressX.center, stressX.shearVec2, olc::GREEN, "G13");
        DrawStressVector(mvp, stress_nX.center, stress_nX.shearVec1, olc::GREEN, "G12'");
        DrawStressVector(mvp, stress_nX.center, stress_nX.shearVec2, olc::GREEN, "G13'");

        DrawStressVector(mvp, stressY.center, stressY.shearVec1, olc::GREEN, "G21");
        DrawStressVector(mvp, stressY.center, stressY.shearVec2, olc::GREEN, "G23");
        DrawStressVector(mvp, stress_nY.center, stress_nY.shearVec1, olc::GREEN, "G21'");
        DrawStressVector(mvp, stress_nY.center, stress_nY.shearVec2, olc::GREEN, "G23'");

        DrawStressVector(mvp, stressZ.center, stressZ.shearVec1, olc::GREEN, "G31");
        DrawStressVector(mvp, stressZ.center, stressZ.shearVec2, olc::GREEN, "G32");
        DrawStressVector(mvp, stress_nZ.center, stress_nZ.shearVec1, olc::GREEN, "G31'");
        DrawStressVector(mvp, stress_nZ.center, stress_nZ.shearVec2, olc::GREEN, "G32'");

        // Finish drawing on the canvas and apply filter
        SetDrawTarget(nullptr);

        drawCrt();
    }

    /**
        * @brief Draws the 2D parametric representation of Mohr's Circle.
        * 
        * Represents the principal stresses ($\sigma_1, \sigma_2, \sigma_3$), calculates maximum shears
        * and prints the Von Mises criterion based on:
        * $$ \sigma_v = \sqrt{\frac{(\sigma_1 - \sigma_2)^2 + (\sigma_2 - \sigma_3)^2 + (\sigma_3 - \sigma_1)^2}{2}} $$
        * 
        * @param[in] originX Origin position on screen (X).
        * @param[in] originY Origin position on screen (Y).
        * @param[in] width   Width of the assigned area.
        * @param[in] height  Height of the assigned area.
        *
    **/
    void DrawMohrDiagram(int originX, int originY, int width, int height) {

        float offset = 430; 
        float centerY = 400;

        float s1 = tensorState.sigma1;
        float s2 = tensorState.sigma2;
        float s3 = tensorState.sigma3;

        // (Circle of radius 0)
        if (std::abs(s1 - s3) < 1e-5f) return;

        float maxradiusPx = 130.0f; 

        float screencenterX = offset + maxradiusPx + 30.0f; 

        float centerStress = (s1 + s3) / 2.0f;       
        float maxRadiusStress = std::abs(s1 - s3) / 2.0f; 

        // Dynamic scale factor from Stress to Pixels
        float scale = maxradiusPx / maxRadiusStress;

        auto toscreenX = [&](float sigma) {
            return screencenterX + (sigma - centerStress) * scale;
        };
  
        float p1X = toscreenX(s1);
        float p2X = toscreenX(s2);
        float p3X = toscreenX(s3);

        float c2X = toscreenX((s1 + s2) / 2.0f); 
        float c3X = toscreenX((s2 + s3) / 2.0f);

        float r1Px = maxradiusPx; 
        float r2Px = std::abs(s1 - s2) / 2.0f * scale;
        float r3Px = std::abs(s2 - s3) / 2.0f * scale;


        float zero_x = std::max(offset, toscreenX(0.0f));
        
        float startX = std::max(offset, toscreenX(std::min(s3, 0.0f)) - 10.0f);
        float endX   = toscreenX(std::max(s1, 0.0f)) + 30.0f;

        // The three Mohr's Circles (Triaxial State) (semi...)
        DrawCircle(screencenterX, centerY, r1Px, olc::Pixel(0, 255, 0), 0xC3);
        DrawCircle(c2X, centerY, r2Px, olc::Pixel(0, 150, 0), 0xC3);
        DrawCircle(c3X, centerY, r3Px, olc::Pixel(0, 150, 0), 0xC3);

        // X and Y axes
        DrawLine(startX, centerY, endX, centerY, olc::Pixel(100, 255, 100), 0xF0000FFF);
        DrawLine(zero_x, centerY - maxradiusPx + 5, zero_x, centerY + maxradiusPx-50, olc::Pixel(100, 255, 100), 0xF0000FFF);

        DrawString(zero_x + 3, centerY + 5, "0", olc::Pixel(150, 255, 150), 1);

        // Principal Stresses Marks
        DrawString(p1X - 15, centerY + 10, Format2Dec(s1), olc::Pixel(100, 200, 100), 1);
        DrawString(p2X - 15, centerY + 10, Format2Dec(s2), olc::Pixel(100, 200, 100), 1);
        DrawString(p3X - 15, centerY + 10, Format2Dec(s3), olc::Pixel(100, 200, 100), 1);

        FillCircle(p1X, centerY, 3, olc::GREEN);
        FillCircle(p2X, centerY, 3, olc::GREEN);
        FillCircle(p3X, centerY, 3, olc::GREEN);

        // Shear Stress Calculations
        float tau12 = std::abs(s1 - s2) / 2.0f;
        float tau23 = std::abs(s2 - s3) / 2.0f;
        float tau13 = std::abs(s1 - s3) / 2.0f; 

        // Von Mises equivalent stress
        float vonMises = std::sqrt((std::pow(s1 - s2, 2) + std::pow(s2 - s3, 2) + std::pow(s3 - s1, 2)) / 2.0f);

        int startY = 485;
        int col1X = 30;
        int col2X = 410;
        int spacing = 20;

        // Summary Printout
        DrawString(col1X, startY, "=== STRESS STATE SUMMARY ===", olc::Pixel(0, 200, 0), 2);
        startY += spacing + 5;

        DrawString(col1X, startY,               "Sigma 1 (Max) : " + Format2Dec(s1), olc::Pixel(100, 255, 100), 2);
        DrawString(col1X, startY + spacing,     "Sigma 2 (Med) : " + Format2Dec(s2), olc::Pixel(100, 255, 100), 2);
        DrawString(col1X, startY + spacing * 2, "Sigma 3 (Min) : " + Format2Dec(s3), olc::Pixel(100, 255, 100), 2);

        DrawString(col1X, startY + spacing * 3, "Von Mises     : " + Format2Dec(vonMises), olc::Pixel(180, 255, 50), 2);

        DrawString(col2X, startY,               "Tau (1-2)     : " + Format2Dec(tau12), olc::Pixel(100, 255, 100), 2);
        DrawString(col2X, startY + spacing,     "Tau (2-3)     : " + Format2Dec(tau23), olc::Pixel(100, 255, 100), 2);
        DrawString(col2X, startY + spacing * 2, "Tau MAX (1-3) : " + Format2Dec(tau13), olc::Pixel(180, 255, 50),  2);
    }

public:
    visualizer3D() {
        sAppName = "Mohr's Circle 3D";
    }

    /**
        * @brief Initializes the main components of the program.
        * 
        * Defines the initial tensor, initializes the visual matrix, configures the 3D pipeline (camera)
        * and precalculates the look-up table (LUT) to optimize the CRT filter.
        * 
        * @return boolean indicating the success of the creation.
        *
    **/
    bool OnUserCreate() override {
        // Initial tensor with default values
        tensorState.sigma << 150.0f,  100.0f,  -55.0f,
                             100.0f,  150.0f,  -95.0f,
                             -55.0f,  -95.0f,  150.0f;

        sigmaS.resize(3);
        for (int f = 0; f < 3; ++f) {
            sigmaS[f].resize(3);
        }
        
        for (int f = 0; f < 3; ++f) {
            for (int c = 0; c < 3; ++c) {
                char buffer[16];
                std::snprintf(buffer, sizeof(buffer), "%7.2f", tensorState.sigma(f, c));
                sigmaS[f][c] = buffer;

            }
        }

        screenHeight    = ScreenHeight();  screenWidth = ScreenWidth();
        hscreenHeight   = screenHeight >> 1; hscreenWidth = screenWidth >> 1;
        inscreenWidth2  = 2.0f / static_cast<float>(screenWidth);
        inscreenHeight2 = 2.0f / static_cast<float>(screenHeight);
        syPos = hscreenHeight; sxPos = (hscreenWidth >> 1);

        // Camera initialization (Pipeline)
        rendpipeLine = pipeLine(screenWidth, screenHeight, pi/2, 0.1f, 50.0f, Eigen::Vector3f(0.0f, 0.0f, -25.0f), Eigen::Vector3f(0.0f, 0.0f, 0.0f), 6.5f);

        canvaSprite = std::make_unique<olc::Sprite>(screenWidth,screenHeight);

        crtLut.resize(screenWidth*screenHeight);

        // Pre-calculation of the CRT distortion map
        for(int y = 0; y < screenHeight; y++){
            float ny = (static_cast<float>(y)*inscreenHeight2) - 1.0f; // normalized y coordinate //
            float sm = (y % 2 != 0) ? 0.55f : 1.0f;
            float ny2 = ny * ny;

            for(int x = 0; x < screenWidth; x++){
                float nx = (static_cast<float>(x)*inscreenWidth2) - 1.0f; // normalized x coordinate //

                float r2 = (nx*nx)+ny2;

                float dt = (1 + k*r2);
                float ux = ((nx*dt) + 1) * (hscreenWidth);
                float uy = ((ny*dt) + 1) * (hscreenHeight);

                if (ux >= 0 && ux < screenWidth && uy >= 0 && uy < screenHeight){
                    crtLut[y*screenWidth+x] = {static_cast<int>(ux), static_cast<int>(uy), std::max(0.0f, 1.0f - (r2 * 0.45f)) * sm};
                } 
                else{
                    crtLut[y*screenWidth+x] = {-1, -1, 0.0f}; 
                }
            }
        }

        tensorState.solve();

        return true;
    }
    
    /**
        * @brief Main engine loop, called once per frame.
        * 
        * @param[in] fElapsedTime Time elapsed since the last frame (Delta Time).
        * @return boolean true.
        *
    **/
    bool OnUserUpdate(float fElapsedTime) override {
        
        handleInput(fElapsedTime);
        updateLogic(fElapsedTime);
        renderScene();

        return true;
    }

    /**
        * @brief Applies a CRT filter based on the precalculated LUT.
        * 
        * Extracts the rendered pixels in `canvaSprite`, calculates grayscale luminance
        * with a retro terminal green tint, and applies spatial mapping (Brown-Conrady distortion).
    **/
    void drawCrt(){
        olc::Pixel* pSrc = canvaSprite->GetData();
        olc::Pixel* pDest = GetDrawTarget()->GetData();

        for(int y = 0; y < screenHeight; y++){
            for(int x = 0; x < screenWidth; x++){
                crtMap map = crtLut[y*screenWidth+x];
            
                if (map.ux >= 0) {   
                    olc::Pixel c = pSrc[map.uy*screenWidth+map.ux];
                    
                    // Fast calculation of relative luminance
                    float luma = (0.2126f * c.r) + (0.7152f * c.g) + (0.0722f * c.b);
                    int final_g = std::min(255, static_cast<int>(map.vm * luma * 1.3f));
                    
                    pDest[y*screenWidth+x] = olc::Pixel(0, final_g, 0);
                }
            }
        }
    }

};

/**
    * @brief Main entry point of the application.
**/
int main() {
    visualizer3D visualizador;
    // Start with 800, 600 resolution, 1x1 pixels
    if (visualizador.Construct(800, 600, 1, 1)) {
        visualizador.Start();
    }
    return 0;
}

//      |\__/,|   (`\_
//    *.|o o  |*   ) )
//---(((---(((------------------
//|                            |
//|       Mohr Visualizer      |
//|____________________________|