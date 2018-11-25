in vec4 uvcoordsvar;

out vec4 FragColor;

uniform sampler2D strokeColor;
uniform sampler2D strokeDepth;
uniform sampler2D blendColor;
uniform sampler2D blendDepth;
uniform int mode;
uniform int disable_mask;

#define THRESHOLD 0.001f
#define ON 1
#define OFF 0

#define MODE_NORMAL   0
#define MODE_OVERLAY  1
#define MODE_ADD      2
#define MODE_SUB      3
#define MODE_MULTIPLY 4
#define MODE_DIVIDE   5

float overlay_color(float a, float b)
{
	float rtn;
		if (a < 0.5) {
			rtn = 2.0 * a * b;
		}
		else {
			rtn = 1.0 - 2.0 * (1.0 - a) * (1.0 - b);
		}

	return rtn;
}

vec4 get_blend_color(int mode, vec4 src_color, vec4 blend_color)
{
	vec4 mix_color = blend_color;
	vec4 outcolor;

	if (mix_color.a == 0) {
		outcolor = src_color;
	}
	else if (mode == MODE_OVERLAY) {
		mix_color.rgb = mix_color.rgb * mix_color.a;
		outcolor.r = overlay_color(src_color.r, mix_color.r);
		outcolor.g = overlay_color(src_color.g, mix_color.g);
		outcolor.b = overlay_color(src_color.b, mix_color.b);
		outcolor.a = src_color.a;
	}
	else if (mode == MODE_ADD){
		mix_color.rgb = mix_color.rgb * mix_color.a;
		outcolor = src_color + mix_color;
		outcolor.a = src_color.a;
	}
	else if (mode == MODE_SUB){
		outcolor = src_color - mix_color;
		outcolor.a = clamp(src_color.a - mix_color.a, 0.0, 1.0);
	}
	else if (mode == MODE_MULTIPLY)	{
		mix_color.rgb = mix_color.rgb * mix_color.a;
		outcolor = src_color * mix_color;
		outcolor.a = src_color.a;
	}
	else if (mode == MODE_DIVIDE) {
		mix_color.rgb = mix_color.rgb * mix_color.a;
		outcolor = src_color / mix_color;
		outcolor.a = src_color.a;
	}
	else {
		outcolor = mix_color;
		outcolor.a = src_color.a;
	}
	
	return outcolor;
}

vec4 get_blurcolor(ivec2 uv, int limit)
{
	int pixels = ((limit * 2) + 1) * ((limit * 2) + 1);
	vec4 blend_color = vec4(0, 0, 0, 0);
	vec4 color;
	int hit = 0;
	for (int x = -limit; x < limit + 1; x++) {
		for (int y = -limit; y < limit + 1; y++) {
			color = texelFetch(blendColor, uv + ivec2(x, y), 0).rgba;
			if (color.a > THRESHOLD) {
				hit++;
			}
			blend_color += color;
		}
		
	}
	/* color is the result of the visible pixel */
	if (hit > 0) {
		blend_color.rgb = blend_color.rgb / hit;
	}
	/* alpha is divided by numberof pixels */
	blend_color.a = blend_color.a / pixels;
	

	return blend_color;
}

void main()
{
	vec4 outcolor;
	ivec2 uv = ivec2(gl_FragCoord.xy);
	vec4 stroke_color =  texelFetch(strokeColor, uv, 0).rgba;
	float stroke_depth = texelFetch(strokeDepth, uv, 0).r;
	if ((stroke_color.a == 0) && (disable_mask == ON)) {
		discard;
	}
	
	vec4 mix_color =  texelFetch(blendColor, uv, 0).rgba;
	float mix_depth = texelFetch(blendDepth, uv, 0).r;

	/* premult alpha factor to remove double blend effects */
	if (stroke_color.a > 0) {
		stroke_color = vec4(vec3(stroke_color.rgb / stroke_color.a), stroke_color.a);
	}
	if (mix_color.a > 0) {
		mix_color = vec4(vec3(mix_color.rgb / mix_color.a), mix_color.a);
	}
	
	/* Normal mode */
	if (mode == MODE_NORMAL) {
		if (mix_color.a > THRESHOLD) {
			outcolor =  mix_color;
			gl_FragDepth = mix_depth;
		}
		else {
			/* blur edges with box blur */
			vec4 blur = get_blurcolor(uv, 1);
			/* interpolate color using the alpha factor */
			outcolor = vec4(mix(stroke_color.rgb, blur.rgb, blur.a), stroke_color.a); 
			gl_FragDepth = stroke_depth;
		}
		FragColor = outcolor;
		return;
	}
	
	/* if not using mask, return mix color */
	if ((stroke_color.a == 0) && (disable_mask == OFF)) {
		FragColor = mix_color;
		gl_FragDepth = mix_depth;
		return;
	}

	/* apply blend mode */
	FragColor = get_blend_color(mode, stroke_color, mix_color);
	gl_FragDepth = stroke_depth;
}
