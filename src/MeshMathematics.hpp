#pragma once

#include <vector>

#include "glad/glad.h"
#include "glm/glm.hpp"

//Step1: Calculate the QEM using Quadratice Error

//Calculates the equation of plane
glm::vec4 calculatePlane(const glm::vec3& v1, const glm::vec3& v2, const glm::vec3& v3) {
    glm::vec3 vec1 = v2 - v1;
    glm::vec3 vec2 = v3 - v1;
    
    glm::vec3 normal = glm::cross(vec1, vec2);
    
    normal = glm::normalize(normal);
    float D = -glm::dot(normal, v1);
    
    return glm::vec4(normal, D);
}

//Calculates the Quadric of a plane
glm::mat4 computeQuadric(const glm::vec4& planeEquation) {
    float A = planeEquation.x;
    float B = planeEquation.y;
    float C = planeEquation.z;
    float D = planeEquation.w;

    glm::mat4 quadric(
        A * A, A * B, A * C, A * D,
        A * B, B * B, B * C, B * D,
        A * C, B * C, C * C, C * D,
        A * D, B * D, C * D, D * D
    );

    return quadric;
}

