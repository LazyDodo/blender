in vec3 pos;
in vec3 N1;
in vec3 N2;
flat out vec3 faceNormal1;
flat out vec3 faceNormal2;
void main()
{
	/* Pass through, MVP calculation happens in geometry shader */
	faceNormal1 = N1;
	faceNormal2 = N2;
	gl_Position = vec4(pos, 1.0);
}
