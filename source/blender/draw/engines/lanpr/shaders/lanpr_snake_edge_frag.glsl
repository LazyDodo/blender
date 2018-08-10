in vec4 uvcoordsvar;
uniform sampler2DMS tex_sampe_0;//depth
uniform sampler2DMS tex_sample_1;//color
uniform sampler2DMS tex_sample_2;//normal
uniform float normal_clamp;//normal clamp
uniform float normal_strength;//normal strength
uniform float depth_clamp;//depth clamp
uniform float depth_strength;//depth strength
uniform float z_near;//z_near
uniform float z_far;//z_far

mat3 sx = mat3(
	1.0, 2.0, 1.0,
	0.0, 0.0, 0.0,
	-1.0, -2.0, -1.0
	);
mat3 sy = mat3(
	1.0, 0.0, -1.0,
	2.0, 0.0, -2.0,
	1.0, 0.0, -1.0
	);
vec3 rgb2hsv(vec3 c)
{
	vec4 K = vec4(0.0, -1.0 / 3.0, 2.0 / 3.0, -1.0);
	vec4 p = mix(vec4(c.bg, K.wz), vec4(c.gb, K.xy), step(c.b, c.g));
	vec4 q = mix(vec4(p.xyw, c.r), vec4(c.r, p.yzx), step(p.x, c.r));
	float d = q.x - min(q.w, q.y);
	float e = 1.0e-10;
	return vec3(abs(q.z + (q.w - q.y) / (6.0 * d + e)), d / (q.x + e), q.x);
}
float linearDepth(float depthSample){
	float d = 2.0 * depthSample - 1.0;
	float zLinear = 2.0 * z_near * z_far / (z_far + z_near - d * (z_far - z_near));
	return zLinear;
}

vec4 DetectEdge(sampler2DMS tex, float clamp, float strength){
	mat3 I = mat3(0);
	mat3 J = mat3(0);
	mat3 K = mat3(0);

	ivec2 texSize = textureSize(tex);
	ivec2 sp = ivec2(uvcoordsvar.xy * texSize);
	vec4 cs = vec4(0);

	//sample hardcoded (8) now
	for (int s = 0; s < 8; s++) {
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				vec4 col = texelFetch(tex, sp + ivec2(i - 1, j - 1), s);
				vec3 sample1 = vec3((col.r));
				vec3 sample2 = vec3((col.g));
				vec3 sample3 = vec3((col.b));
				I[i][j] += length(sample1) / 8;
				J[i][j] += length(sample2) / 8;
				K[i][j] += length(sample3) / 8;
			}
		}
		cs += texelFetch(tex_sample_1, sp, s) / 8;
	}

	float gx1 = dot(sx[0], I[0]) + dot(sx[1], I[1]) + dot(sx[2], I[2]);
	float gy1 = dot(sy[0], I[0]) + dot(sy[1], I[1]) + dot(sy[2], I[2]);
	float g1 = sqrt(pow(gx1, 2.0) + pow(gy1, 2.0));

	float gx2 = dot(sx[0], J[0]) + dot(sx[1], J[1]) + dot(sx[2], J[2]);
	float gy2 = dot(sy[0], J[0]) + dot(sy[1], J[1]) + dot(sy[2], J[2]);
	float g2 = sqrt(pow(gx2, 2.0) + pow(gy2, 2.0));

	float gx3 = dot(sx[0], K[0]) + dot(sx[1], K[1]) + dot(sx[2], K[2]);
	float gy3 = dot(sy[0], K[0]) + dot(sy[1], K[1]) + dot(sy[2], K[2]);
	float g3 = sqrt(pow(gx3, 2.0) + pow(gy3, 2.0));

	float value = max(max(g1, g2), g3);

	value = value > clamp ? value : 0;
	return vec4(vec3(value * strength), 1);

	//if(value<clamp) value=0;

	//return vec4(pow(value,strength));
}

void main()
{
	float nc = normal_clamp;//(normal_clamp==0? 0.01:normal_clamp);
	float ns = normal_strength;//(normal_strength==0? 5:   normal_strength);
	float dc = depth_clamp;//(depth_clamp==0? 0.2: depth_clamp);
	float ds = depth_strength;//(depth_strength==0? 2.5: depth_strength);

	vec4 diffuse = vec4(1, 1, 1, 1);
	vec4 color = (DetectEdge(tex_sampe_0, dc, ds) + DetectEdge(tex_sample_2, nc, ns));
	gl_FragColor = color;
};