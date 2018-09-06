in vec4 pos;
in vec2 uvs;
in vec3 normal;
in int type;
in int level;

uniform vec4 output_viewport;
uniform vec4 preview_viewport;

out vec2 gOffset;
out int gType;
out int gLevel;
out vec3 gNormal;

void main(){
	vec4 p = pos;
	// move to geo shader
	//float asp1 = output_viewport.z / output_viewport.w;
	//float asp2 = preview_viewport.z / preview_viewport.w;
	//p.x = pos.x / asp2 * asp1;

	gOffset = uvs;
	gType = type;
	gLevel = level;
	gNormal = normal;
	gl_Position = vec4(vec3(p), 1);
}
