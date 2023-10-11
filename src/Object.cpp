
#include "Object.hpp"

Object::Object(const std::string& object_file_path) {
    

    std::string line = "";
    std::vector<glm::vec3> normal;
    this->facePattern = -1;
    std::ifstream myFile(object_file_path.c_str());
   
    std::cout <<"File reading has started" << std::endl;

    if (myFile.is_open())
    {
        while (std::getline(myFile, line))
        {
            
            //std::cout << line << std::endl;

            if (line[0] != '#')
            {
                if (line[0] == 'v' && line[1] != 'n')
                {
                    std::stringstream ss(line);

                    char v;
                    GLfloat num1, num2, num3;

                    ss >> v >> num1 >> num2 >> num3;

                    this->vertices_.push_back(Vertex(glm::vec3(num1, num2, num3)));
                   
                } else if (line[0] == 'v' && line[1] == 'n')
                {
                    
                    std::stringstream ss(line);

                    char v, n;
                    GLfloat num1, num2, num3;

                    ss >> v >> n >> num1 >> num2 >> num3;

                    normal.push_back(glm::vec3(num1, num2, num3));

                } else if (line[0] == 'f')
                {
                    
                    //determine the pattern
                    if (this->facePattern == -1)
                    {
                        std::regex vtn_pattern(R"(\d+/\d+/\d+)");

                        if (std::regex_search(line, std::regex("//")))
                        {
                            //Then the f is like v//n
                            this->facePattern = 1;
                        } else if (std::regex_search(line, vtn_pattern))
                        {
                            //Then the f is like v/t/n
                            this->facePattern = 2;
                        } else {
                            this->facePattern = 0;
                        }

                        std::cout << "Face pattern: " << this->facePattern << std::endl;    
                    }

                    //Now use the pattern to read the face values
                    if (this->facePattern == 0)
                    {
                        std::stringstream ss(line);

                        char v;
                        GLuint num1, num2, num3;

                        ss >> v >> num1 >> num2 >> num3;

                        this->indices_.push_back(num1 - 1);
                        this->indices_.push_back(num2 - 1);
                        this->indices_.push_back(num3 - 1);
                    } else if (this->facePattern == 1)
                    {
                        std::stringstream ss(line);

                        char f;
                        GLuint v1, n1, v2, n2, v3, n3;
                        char delim; // Used to skip the //

                        ss >> f; // Read the 'f' character

                        // Read the first vertex/normal pair
                        ss >> v1 >> delim >> delim >> n1;

                        // Read the second vertex/normal pair
                        ss >> v2 >> delim >> delim >> n2;

                        // Read the third vertex/normal pair
                        ss >> v3 >> delim >> delim >> n3;
 
                        this->indices_.push_back(v1 - 1);
                        this->vertices_[v1 - 1].updateNormal(normal[n1 -1]);
                        this->indices_.push_back(v2 - 1);
                        this->vertices_[v2 - 1].updateNormal(normal[n2 -1]);
                        this->indices_.push_back(v3 - 1);
                        this->vertices_[v3 - 1].updateNormal(normal[n3 -1]);   
                    } else {
                        std::cout << "Error: Format of face not recognized" << std::endl;
                        //Write a better way to handle it
                    }
                    

                }
                
                
                
            }
            
        }
        myFile.close();

        std::cout <<"File reading has finished" << std::endl;   
    }

    
}

int Object::getFacePattern() {
    return this->facePattern;
}

const std::vector<GLfloat> Object::vertices() {
    std::vector<GLfloat> vertexInfo;
    
    for (auto vertex : this->vertices_)
    {
        vertexInfo.push_back(vertex.position_.x);
        vertexInfo.push_back(vertex.position_.y);
        vertexInfo.push_back(vertex.position_.z);

        if (facePattern == 1)
        {
            vertexInfo.push_back(vertex.normal_.x);
            vertexInfo.push_back(vertex.normal_.y);
            vertexInfo.push_back(vertex.normal_.z);
        }
        
    }

    return vertexInfo;    
}

const std::vector<GLuint> Object::indices() {
    return this->indices_;  
}