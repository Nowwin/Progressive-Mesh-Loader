#pragma once

#include <iostream>
#include <vector>
#include <filesystem>

#include <SDL.h>
#include <glad/glad.h>
#include <glm/vec3.hpp>
#include <glm/vec4.hpp>
#include <glm/mat4x4.hpp>
#include <glm/ext/matrix_transform.hpp>
#include <glm/ext/matrix_clip_space.hpp>
#include <glm/ext/scalar_constants.hpp>


class GraphicsApp {
public:
    GraphicsApp(int width, int height, float offset, float rotation) {
        m_screenWidth = width;
        m_screenHeight = height;
        m_uOffset = offset;
        m_uRotation = rotation;
        gQuit = false;
        gGraphicsApplicationWindow = nullptr;
        gOpenGlContext = nullptr;
        gVertexArrayObject = 0;
        gVertexBufferObject = 0;
        gIndexBufferObject = 0;
        gGraphicsPipelineShaderProgram = 0;
    };
    ~GraphicsApp() {
    };

    void Initialize();
    void MainLoop();
    void CleanUp();

private:
    // Private helper methods
    void GetOpenGLVersionInfo();
    void IntializeProgram();
    void VertexSpecification();
    void LoadShaders();
    void Input();
    void PreDraw();
    void Draw();

    // Private member variables
    bool gQuit;
    SDL_Window* gGraphicsApplicationWindow;
    SDL_GLContext gOpenGlContext;
    GLuint gVertexArrayObject, gVertexBufferObject, gIndexBufferObject, gGraphicsPipelineShaderProgram;
    float m_uOffset, m_uRotation;
    int m_screenWidth, m_screenHeight;
};
