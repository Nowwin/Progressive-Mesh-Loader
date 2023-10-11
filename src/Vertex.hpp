#pragma once

#include "glad/glad.h"
#include "glm/glm.hpp"

//A vertex can have a position and normal
struct Vertex {
    glm::vec3 position_;
    glm::vec3 normal_;

    Vertex(const glm::vec3& position, const glm::vec3& normal)
        : position_(position), normal_(normal){};

    Vertex(const glm::vec3& position)
        : position_(position), normal_(glm::vec3(0.0f)){};

    void updateNormal(const glm::vec3& normal) {
        normal_ = normal;
    } 

    void updatePosition(const glm::vec3& position) {
        position_ = position;
    }    
  };