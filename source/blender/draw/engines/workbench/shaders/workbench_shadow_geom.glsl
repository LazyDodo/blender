layout(triangles) in;
layout(triangle_strip, max_vertices=3) out;

uniform vec4 direction = vec4(0.57, 0.57, 0.0, 0.0);

void main()
{
	for(int i = 0; i < gl_in.length(); i++)
	{
		vec4 new_pos = gl_in[i].gl_Position;

		new_pos += direction;
		gl_Position = new_pos;
		EmitVertex();
	}
}