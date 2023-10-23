#include <MeshGeometry.hpp>

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

void Mesh::AddEdgeInfo() {

}