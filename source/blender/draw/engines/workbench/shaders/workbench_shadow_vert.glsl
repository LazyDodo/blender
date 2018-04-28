uniform mat4 ModelViewProjectionMatrix;
uniform mat3 NormalMatrix;
in vec3 pos;
in vec3 N1;
in vec3 N2;
out vec3 faceNormal1;
out vec3 faceNormal2;
void main()
{
	gl_Position = ModelViewProjectionMatrix * vec4(pos, 1.0);
	faceNormal1 = NormalMatrix * N1;
	faceNormal2 = NormalMatrix * N2;
}
