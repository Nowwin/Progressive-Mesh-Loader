#version 410 core


layout(location=0) in vec3 position;
layout(location = 1) in vec3 inNormal;

out vec3 fragNormal;

uniform mat4 u_ModelMatrix;
uniform mat4 u_PerspectiveMatrix;

void main()
{
    fragNormal = inNormal;

    vec4 newPosition = u_PerspectiveMatrix * u_ModelMatrix * vec4(position, 1.0f);
    
    gl_Position = vec4(newPosition.x, newPosition.y, newPosition.z, newPosition.w);
}