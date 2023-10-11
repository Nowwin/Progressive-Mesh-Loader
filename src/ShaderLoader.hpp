#pragma once
#include <string>

#include <glad/glad.h>

class ShaderLoader {
public:
    static std::string LoadShaderAsString(const std::string& filename);
    static GLuint CompileShader(GLuint type, const std::string &source);
    static GLuint CreateShaderProgram(const std::string &vertexShaderSource, const std::string &fragmentShaderSource);
};
