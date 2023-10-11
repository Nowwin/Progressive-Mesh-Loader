#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <regex>

#include "glad/glad.h"
#include "glm/glm.hpp"

#include "Vertex.hpp"

class Object {
 public:
  Object() {};
  Object(const std::string& object_file_path);
  ~Object() {};

  int getFacePattern();
  const std::vector<GLfloat> vertices();
  const std::vector<GLuint> indices();

 private:
  std::vector<Vertex> vertices_;
  std::vector<GLuint> indices_;
  int facePattern;
};
