layout(triangles) in;
layout(triangle_strip, max_vertices=12) out;

uniform vec4 lightDirection = vec4(-4.0, -25.0, 0.0, 0.0);

void main()
{
	gl_Position = gl_in[0].gl_Position;
	EmitVertex();
	gl_Position = gl_in[1].gl_Position;
	EmitVertex();
	gl_Position = gl_in[0].gl_Position + lightDirection;
	EmitVertex();
	gl_Position = gl_in[1].gl_Position + lightDirection;
	EmitVertex();
	EndPrimitive();

	gl_Position = gl_in[1].gl_Position;
	EmitVertex();
	gl_Position = gl_in[2].gl_Position;
	EmitVertex();
	gl_Position = gl_in[1].gl_Position + lightDirection;
	EmitVertex();
	gl_Position = gl_in[2].gl_Position + lightDirection;
	EmitVertex();
	EndPrimitive();

	gl_Position = gl_in[2].gl_Position;
	EmitVertex();
	gl_Position = gl_in[0].gl_Position;
	EmitVertex();
	gl_Position = gl_in[2].gl_Position + lightDirection;
	EmitVertex();
	gl_Position = gl_in[0].gl_Position + lightDirection;
	EmitVertex();
	EndPrimitive();
}