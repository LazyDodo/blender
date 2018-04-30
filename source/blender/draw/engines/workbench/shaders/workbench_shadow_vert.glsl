in vec3 pos;
in vec3 N1;
in vec3 N2;
out vec3 faceNormal1;
out vec3 faceNormal2;
void main()
{
	faceNormal1 = N1;
	faceNormal2 = N2;
	gl_Position = vec4(pos, 1.0);
}
