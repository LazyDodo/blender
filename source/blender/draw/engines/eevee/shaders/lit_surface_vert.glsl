
uniform mat4 ModelViewProjectionMatrix;
uniform mat4 ModelMatrix;
uniform mat4 ModelViewMatrix;
uniform mat4 ModelViewMatrixInverse;
uniform mat3 WorldNormalMatrix;
#ifndef ATTRIB
uniform mat3 NormalMatrix;
#endif

#ifdef HAIR_SHADER

#ifdef HAIR_SHADER_FIBERS
in int fiber_index;
in float curve_param;
#endif

#else
in vec3 pos;
in vec3 nor;
#endif

out vec3 worldPosition;
out vec3 viewPosition;

/* Used for planar reflections */
/* keep in sync with EEVEE_ClipPlanesUniformBuffer */
layout(std140) uniform clip_block {
	vec4 ClipPlanes[1];
};

#ifdef USE_FLAT_NORMAL
flat out vec3 worldNormal;
flat out vec3 viewNormal;
#else
out vec3 worldNormal;
out vec3 viewNormal;
#endif

#ifdef HAIR_SHADER
out vec3 hairTangent;
out float hairThickTime;
out float hairThickness;
out float hairTime;
flat out int hairStrandID;
#endif

void main()
{
#ifdef HAIR_SHADER
	bool is_persp = (ProjectionMatrix[3][3] == 0.0);

#ifdef HAIR_SHADER_FIBERS
	vec3 pos, tang, binor;
	hair_fiber_get_vertex(
	        fiber_index, curve_param,
	        is_persp, ModelViewMatrixInverse[3].xyz, ModelViewMatrixInverse[2].xyz,
	        pos, tang, binor,
	        hairTime, hairThickness, hairThickTime);
	vec3 nor = cross(binor, tang);

	gl_Position = ModelViewProjectionMatrix * vec4(pos, 1.0);
	viewPosition = (ModelViewMatrix * vec4(pos, 1.0)).xyz;
	worldPosition = (ModelMatrix * vec4(pos, 1.0)).xyz;
	hairTangent = (ModelMatrix * vec4(tang, 0.0)).xyz;
	worldNormal = (ModelMatrix * vec4(nor, 0.0)).xyz;
	viewNormal = normalize(mat3(ViewMatrix) * worldNormal);
#else
	hairStrandID = hair_get_strand_id();
	vec3 pos, binor;
	hair_get_pos_tan_binor_time(
	        is_persp, ViewMatrixInverse[3].xyz, ViewMatrixInverse[2].xyz,
	        pos, hairTangent, binor, hairTime, hairThickness, hairThickTime);

	gl_Position = ViewProjectionMatrix * vec4(pos, 1.0);
	viewPosition = (ViewMatrix * vec4(pos, 1.0)).xyz;
	worldPosition = pos;
	hairTangent = normalize(hairTangent);
	worldNormal = cross(binor, hairTangent);
	viewNormal = normalize(mat3(ViewMatrix) * worldNormal);
#endif

#else
	gl_Position = ModelViewProjectionMatrix * vec4(pos, 1.0);
	viewPosition = (ModelViewMatrix * vec4(pos, 1.0)).xyz;
	worldPosition = (ModelMatrix * vec4(pos, 1.0)).xyz;
	worldNormal = normalize(WorldNormalMatrix * nor);
	viewNormal = normalize(NormalMatrix * nor);
#endif

	/* Used for planar reflections */
	gl_ClipDistance[0] = dot(vec4(worldPosition, 1.0), ClipPlanes[0]);

#ifdef ATTRIB
	pass_attrib(pos);
#endif
}
