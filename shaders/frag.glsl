#version 410 core
in vec3 fragNormal;

out vec4 color;

uniform bool isWireframeMode;

void main()
{
    
    if(isWireframeMode) {
        color = vec4(1.0, 1.0, 1.0, 1.0);
    } else {
        vec3 normalizedNormal = normalize(fragNormal);
        color = vec4(normalizedNormal * 0.5 + 0.5, 1.0);
    }
    
    
}