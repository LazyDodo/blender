#define EPSILON 0.0000001
#define INFINITE 100.0

layout(triangles) in;
layout(triangle_strip, max_vertices=14) out;

uniform mat4 ModelMatrix;
uniform mat4 ModelViewProjectionMatrix;

uniform vec3 lightDirection = vec3(0.57, 0.57, -0.57);
vec3 face_normal(vec3 v1, vec3 v2, vec3 v3) {
	return normalize(cross(v2 - v1, v3 - v1));
}
void emit_side_quad(vec4 v1, vec4 v2, vec4 light_direction) {
	gl_Position = ModelViewProjectionMatrix * (v1 + light_direction * EPSILON);
	EmitVertex();
	gl_Position = ModelViewProjectionMatrix * (v1 + light_direction * INFINITE);
	EmitVertex();
	gl_Position = ModelViewProjectionMatrix * (v2 + light_direction * EPSILON);
	EmitVertex();
	gl_Position = ModelViewProjectionMatrix * (v2 + light_direction * INFINITE);
	EmitVertex();
	EndPrimitive();
}
void emit_front_cap(vec4 v1, vec4 v2, vec4 v3, vec4 light_direction) {
	gl_Position = ModelViewProjectionMatrix * (v1 + light_direction * EPSILON);
	EmitVertex();
	gl_Position = ModelViewProjectionMatrix * (v2 + light_direction * EPSILON);
	EmitVertex();
	gl_Position = ModelViewProjectionMatrix * (v3 + light_direction * EPSILON);
	EmitVertex();
	EndPrimitive();
}

void emit_back_cap(vec4 v1, vec4 v2, vec4 v3, vec4 light_direction) {
	gl_Position = ModelViewProjectionMatrix * (v3 + light_direction * INFINITE);
	EmitVertex();
	gl_Position = ModelViewProjectionMatrix * (v2 + light_direction * INFINITE);
	EmitVertex();
	gl_Position = ModelViewProjectionMatrix * (v1 + light_direction * INFINITE);
	EmitVertex();
	EndPrimitive();
}

void main()
{
	/* light_direction: Light direction in object space TODO: Move to world space */
	vec3 light_direction = normalize((vec4(lightDirection, 0.0) * ModelMatrix).xyz);
	vec4 v1 = gl_in[0].gl_Position;
	vec4 v2 = gl_in[1].gl_Position;
	vec4 v3 = gl_in[2].gl_Position;
	bool backface = dot(face_normal(v1.xyz, v2.xyz, v3.xyz), light_direction) > 0.0;

	if (backface) {
		vec4 light_direction4 = vec4(light_direction, 0.0);
		/* emit side faces */
		emit_side_quad(v2, v1, light_direction4);
		emit_side_quad(v3, v2, light_direction4);
		emit_side_quad(v1, v3, light_direction4);
		/* emit front and back cap */
		emit_front_cap(v3, v2, v1, light_direction4);
		emit_back_cap(v3, v2, v1, light_direction4);
	}
}
