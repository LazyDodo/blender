#define EPSILON 0.0000001
#define INFINITE 1000.0
layout(lines) in;
layout(triangle_strip, max_vertices=4) out;

uniform mat4 ModelMatrix;
uniform mat4 ModelViewProjectionMatrix;

uniform vec3 lightDirection = vec3(0.57, 0.57, -0.57);
flat in vec3 faceNormal1[]; 
flat in vec3 faceNormal2[];

void emit_quad(vec4 start_vertex, vec4 end_vertex, vec4 light_direction) {
	gl_Position = ModelViewProjectionMatrix * (start_vertex + light_direction * EPSILON);
	EmitVertex();
	gl_Position = ModelViewProjectionMatrix * (start_vertex + light_direction * INFINITE);
	EmitVertex();
	gl_Position = ModelViewProjectionMatrix * (end_vertex + light_direction * EPSILON);
	EmitVertex();
	gl_Position = ModelViewProjectionMatrix * (end_vertex + light_direction * INFINITE);
	EmitVertex();
	EndPrimitive();
}
void main()
{
	/* light_direction: Light direction in object space TODO: Move to CPU */
	vec3 light_direction = normalize((vec4(lightDirection, 0.0) * ModelMatrix).xyz);
	bool front1 = dot(faceNormal1[0], light_direction) >= 0;
	bool front2 = dot(faceNormal2[0], light_direction) >= 0;

	if (!front1 && front2) {
		emit_quad(gl_in[0].gl_Position, gl_in[1].gl_Position, vec4(light_direction, 0.0));
	}
	else if (front1 && !front2) {
		emit_quad(gl_in[1].gl_Position, gl_in[0].gl_Position, vec4(light_direction, 0.0));
	}
}