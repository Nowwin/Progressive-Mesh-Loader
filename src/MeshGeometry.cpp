#include "MeshGeometry.hpp"

bool Mesh::ReadOBJFile(char *filename)
{
    FILE *fp;

    if((fp = fopen(filename, "r")) == NULL ){
        std::cerr << "file cannot be read.\n";
        return false;
    }

    char buf[512];

    std::vector<VertexIter> vertex_iterator;
    int v_id = 0, f_id = 0;
    n_vertices = 0; // Will compute this
    n_faces = 0;    // Will compute this

    while(fgets(buf, 512, fp)) {
        if(buf[0] == 'v' && buf[1] == ' ') { // Geometric vertex
            glm::vec3 coord_in;
            sscanf(buf, "v %f %f %f", &coord_in.x, &coord_in.y, &coord_in.z);

            vertices.push_back(Vertex(coord_in, v_id++));
            vertex_iterator.push_back(--(vertices.end()));
            n_vertices++;
        } 
        else if(buf[0] == 'f') { // Face definition
            int v_id[3];
            if(strstr(buf, "//")) { // format: f v1//vn1 v2//vn2 v3//vn3
                sscanf(buf, "f %d//%*d %d//%*d %d//%*d", &v_id[0], &v_id[1], &v_id[2]);
            } else { // format: f v1 v2 v3
                sscanf(buf, "f %d %d %d", &v_id[0], &v_id[1], &v_id[2]);
            }

            // Update the following line according to the Face constructor signature.
            faces.push_back(Face(vertex_iterator[v_id[0]-1], vertex_iterator[v_id[1]-1], vertex_iterator[v_id[2]-1], f_id++));
            n_faces++;
        }
    }


    std::cerr << "Reading OBJ file done...\n";
    fclose(fp);

    std::cerr << "# of vertices  " << n_vertices << std::endl;
    std::cerr << "# of faces     " << n_faces    << std::endl;

    glm::vec3 range_min(1.0e6f, 1.0e6f, 1.0e6f);
    glm::vec3 range_max(-1.0e6f, -1.0e6f, -1.0e6f);
    glm::vec3 center;

    for (VertexIter vi = vertices.begin(); vi != vertices.end(); vi++) {
        for (int i = 0; i < 3; i++) {
            if (vi->position_[i] < range_min[i])	range_min[i] = vi->position_[i];
            if (vi->position_[i] > range_max[i])	range_max[i] = vi->position_[i];
        }
    }

    // Set input model within (-1,-1,-1) and (1,1,1) and the centroid is at the origin
    center = (range_min + range_max) * 0.5f;

    float largest_range = -1.0f;
    for (int i = 0; i < 3; i++) {
        if (largest_range < range_max[i] - range_min[i]) largest_range = range_max[i] - range_min[i];
    }

    float scale_factor = 2.0f / largest_range;

    for (VertexIter vi = vertices.begin(); vi != vertices.end(); vi++) {
        vi->position_ = (vi->position_ - center) * scale_factor;
    }


    return true;
}

void Mesh::MakeCircularList(FaceIter &fi)
{
    fi->halfedge[0].next = &(fi->halfedge[1]);
    fi->halfedge[1].next = &(fi->halfedge[2]);
    fi->halfedge[2].next = &(fi->halfedge[0]);

    fi->halfedge[0].prev = &(fi->halfedge[2]);
    fi->halfedge[1].prev = &(fi->halfedge[0]);
    fi->halfedge[2].prev = &(fi->halfedge[1]);
}

void Mesh::AssignFaceNormal(FaceIter &fi) {
    glm::vec3 vec1 = fi->halfedge[1].vertex->position_ - fi->halfedge[0].vertex->position_;
    glm::vec3 vec2 = fi->halfedge[2].vertex->position_ - fi->halfedge[0].vertex->position_;

    fi->normal_ = glm::cross(vec1, vec2); // Notice the underscore in normal_
    fi->area = glm::length(fi->normal_) * 0.5f;  // Assuming area is half the length of the cross product
    fi->normal_ = glm::normalize(fi->normal_);   // Notice the underscore in normal_
}

void Mesh::AssignVertexNormal(VertexIter &vi) {
    bool isBoundaryVertex = false;

    vi->normal_ = glm::vec3(0.0f, 0.0f, 0.0f);
    double cumulativeArea = 0.0;

    // Traverse faces incident to "vi" in CCW
    HalfEdge* hep = vi->neighborHe;
    do {
        FaceIter fi = hep->face;

        vi->normal_ += fi->normal_ * static_cast<float>(fi->area);
        cumulativeArea += fi->area;

        hep = hep->prev->mate;

        if (hep == nullptr) {
            isBoundaryVertex = true;
            break;
        }
    } while (hep != vi->neighborHe);

    // When we cannot traverse all incident faces since "vi" is on boundary
    // we traverse faces incident to "vi" in CW to check all incident faces
    if (isBoundaryVertex) {
        HalfEdge* hep = vi->neighborHe->mate;
        while (hep != nullptr) {
            FaceIter fi = hep->face;

            vi->normal_ += fi->normal_ * static_cast<float>(fi->area);
            cumulativeArea += fi->area;

            hep = hep->next->mate;
        }
    }

    vi->normal_ *= (1.0 / cumulativeArea);
}


void Mesh::AddEdgeInfo() {
    // Store face iterators incident to each vertex
    std::vector<std::vector<FaceIter>> Ring(n_vertices);

    for (FaceIter fi = faces.begin(); fi != faces.end(); fi++) {
        // construct circular list
        MakeCircularList(fi);

        for (int i = 0; i < 3; i++) {
            fi->halfedge[i].face = fi;
            fi->halfedge[i].vertex->neighborHe = &(fi->halfedge[i]);

            Ring[fi->halfedge[i].vertex->id].push_back(fi);
        }
    }

    std::cerr << "Halfedges are set\n";

    // Construct mates of halfedge
    for (size_t i = 0; i < Ring.size(); i++) { // For each vertex
        for (size_t j = 0; j < Ring[i].size(); j++) {
            HalfEdge* candidate_he;
            VertexIter candidate_vertex;

            for (int m = 0; m < 3; m++) {
                if (Ring[i][j]->halfedge[m].vertex->id == i) {
                    candidate_he = &(Ring[i][j]->halfedge[m]);
                    candidate_vertex = Ring[i][j]->halfedge[m].next->vertex;
                    break;
                }
            }

            for (size_t k = 0; k < Ring[i].size(); k++) {
                if (j == k) continue;

                for (int m = 0; m < 3; m++) {
                    if (Ring[i][k]->halfedge[m].vertex == candidate_vertex && 
                        Ring[i][k]->halfedge[m].next->vertex->id == i) {
                        candidate_he->mate = &(Ring[i][k]->halfedge[m]);
                        Ring[i][k]->halfedge[m].mate = candidate_he;
                        break;
                    }
                }
            }
        }
    }

    std::cerr << "Halfedge mates are set\n";

    // Add edge information
    for (FaceIter fi = faces.begin(); fi != faces.end(); fi++) {
        for (int i = 0; i < 3; i++) {
            if (fi->halfedge[i].mate == nullptr || 
                fi->halfedge[i].vertex->id < fi->halfedge[i].mate->vertex->id) {
                edges.push_back(Edge(&(fi->halfedge[i]), fi->halfedge[i].mate, n_edges++));
            }
            if (fi->halfedge[i].mate == nullptr) 
                fi->halfedge[i].vertex->isBoundary = true;
        }
    }

    std::cerr << "Edges are set\n";

    // Construct link from halfedge to the corresponding edge
    for (EdgeIter ei = edges.begin(); ei != edges.end(); ei++) {
        ei->halfedge[0]->edge = ei;
        if (ei->halfedge[1] != nullptr) ei->halfedge[1]->edge = ei;
    }

    std::cerr << "# of edges " << n_edges << std::endl;

    for (FaceIter fi = faces.begin(); fi != faces.end(); fi++) 
        AssignFaceNormal(fi);
    for (VertexIter vi = vertices.begin(); vi != vertices.end(); vi++) 
        AssignVertexNormal(vi);
}

bool Mesh::ConstructMeshDataStructure(char *filename)
{
    if( ReadOBJFile(filename) == false ) return false;
    AddEdgeInfo();

    return true;
}

std::vector<GLfloat> Mesh::GetVertexData() {
    std::vector<GLfloat> vertexData;
    vertexIndexMap.clear();
    int newIndex = 0; // New index for each active vertex
    // Iterate over each vertex in the mesh.
    for (const Vertex& vertex : vertices) {
        if (vertex.isActive) {
            vertexIndexMap[vertex.id] = newIndex++;
            // Add vertex position
            vertexData.push_back(vertex.position_.x);
            vertexData.push_back(vertex.position_.y);
            vertexData.push_back(vertex.position_.z);

            // Add vertex normal
            vertexData.push_back(vertex.normal_.x);
            vertexData.push_back(vertex.normal_.y);
            vertexData.push_back(vertex.normal_.z);
        }
    }
    
    return vertexData;
}

std::vector<GLuint> Mesh::GetIndexData() {
    std::vector<GLuint> indexData;
    
    // Iterate over each face in the mesh.
    for (const Face& face : faces) {
        // For each half-edge of the face, retrieve the corresponding vertex's id.
        if (face.isActive) {
        // Assuming triangular faces
            indexData.push_back(vertexIndexMap[face.halfedge[0].vertex->id]);
            indexData.push_back(vertexIndexMap[face.halfedge[1].vertex->id]);
            indexData.push_back(vertexIndexMap[face.halfedge[2].vertex->id]);
        }
    }
    
    return indexData;
}
