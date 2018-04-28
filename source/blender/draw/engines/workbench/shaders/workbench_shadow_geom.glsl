layout(lines) in;
layout(triangle_strip, max_vertices=4) out;

uniform vec4 lightDirection = vec4(-4.0, -25.0, 0.0, 0.0);
uniform vec3 lightDirection2 = vec3(0.57, 0.57, 0.57);
in vec3 faceNormal1[]; 
in vec3 faceNormal2[];

void main()
{
	bool front1 = dot(faceNormal1[0], lightDirection2) >= 0;
	bool front2 = dot(faceNormal2[1], lightDirection2) >= 0;

	if (front1 != front2) {
		gl_Position = gl_in[0].gl_Position;
		EmitVertex();
		gl_Position = gl_in[1].gl_Position;
		EmitVertex();
		gl_Position = gl_in[0].gl_Position + lightDirection;
		EmitVertex();
		gl_Position = gl_in[1].gl_Position + lightDirection;
		EmitVertex();
		EndPrimitive();
	}
}