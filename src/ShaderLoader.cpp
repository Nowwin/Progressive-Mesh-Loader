#include <string>
#include <fstream>

#include <glad/glad.h>

#include "ShaderLoader.hpp"


std::string ShaderLoader::LoadShaderAsString(const std::string& filename) {
    std::string result = "";

    std::string line = "";
    std::ifstream myFile(filename.c_str());

    if (myFile.is_open())
    {
        while (std::getline(myFile, line))
        {
            result += line + '\n';
        }
        myFile.close();   
    }

    

    return result;
    
}


GLuint ShaderLoader::CompileShader(GLuint type, const std::string &source) {
    
    GLuint shaderPbject;
    
    if (type == GL_VERTEX_SHADER)
    {
        shaderPbject = glCreateShader(GL_VERTEX_SHADER);
    } else if(type == GL_FRAGMENT_SHADER) {
        shaderPbject = glCreateShader(GL_FRAGMENT_SHADER);
    }

    const char* src = source.c_str();
    glShaderSource(shaderPbject, 1, &src, nullptr);
    glCompileShader(shaderPbject);

    return shaderPbject;
    
}

GLuint ShaderLoader::CreateShaderProgram(const std::string &vertexshadersource, const std::string &fragmentshadersource) {
    GLuint programObject = glCreateProgram();

    GLuint myVertexShader = CompileShader(GL_VERTEX_SHADER, vertexshadersource);
    GLuint myFragmentShader = CompileShader(GL_FRAGMENT_SHADER, fragmentshadersource);

    glAttachShader(programObject, myVertexShader);
    glAttachShader(programObject, myFragmentShader);
    glLinkProgram(programObject);

    //Validation
    glValidateProgram(programObject);

    return programObject;
}
