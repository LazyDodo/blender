#define EPSILON 0.0000001

layout(lines) in;
layout(triangle_strip, max_vertices=4) out;

uniform mat4 ViewProjectionMatrix;

uniform vec3 lightDirection = vec3(0.57, -0.57, -0.57);
in vec3 faceNormal1[]; 
in vec3 faceNormal2[];
in vec3 localPos[];

void emit_quad(vec3 start_vertex, vec3 end_vertex, vec3 light_direction) {
	gl_Position = ViewProjectionMatrix * vec4((start_vertex + light_direction * EPSILON), 1.0);
	EmitVertex();
	gl_Position = ViewProjectionMatrix * vec4((start_vertex + light_direction * 25.0), 1.0);
	EmitVertex();
	gl_Position = ViewProjectionMatrix * vec4((end_vertex + light_direction * EPSILON), 1.0);
	EmitVertex();
	gl_Position = ViewProjectionMatrix * vec4((end_vertex + light_direction * 25.0), 1.0);
	EmitVertex();
	EndPrimitive();
}
void main()
{
	/* 
		TODO: light is currently connected to the object.
		multiply with inverse object matrix to get to counteracti this.
		or fo not use the modelviewprojectionmatrix (but the ViewProjectionMatrix)
		
		Should the faceNormal and local pos not be done in the vertex shader (conversion to world coordinates...)
	*/
	vec3 light_direction = lightDirection;
	bool front1 = dot(faceNormal1[0], normalize(light_direction)) >= 0;
	bool front2 = dot(faceNormal2[1], normalize(light_direction)) >= 0;

	if (!front1 && front2) {
		emit_quad(gl_in[0].gl_Position.xyz, gl_in[1].gl_Position.xyz, light_direction);
	}
	else if (front1 && !front2) {
		emit_quad(gl_in[1].gl_Position.xyz, gl_in[0].gl_Position.xyz, light_direction);
	}
}