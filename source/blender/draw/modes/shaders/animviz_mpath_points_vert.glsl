
uniform mat4 ViewProjectionMatrix;

uniform int pointSize = 2;
uniform int frameCurrent;
uniform int cacheStart;
uniform bool showKeyFrames = true;
uniform bool useCustomColor;
uniform vec3 customColor;

in vec3 pos;
in int flag;

out vec4 finalColor;

void main()
{
	gl_Position = ViewProjectionMatrix * vec4(pos, 1.0);
	gl_PointSize = float(pointSize + 3);

	int frame = gl_VertexID + cacheStart;
	finalColor = (useCustomColor) ? vec4(customColor, 1.0) : vec4(1.0);

	/* Draw big green dot where the current frame is.
	 * NOTE: this is only done when keyframes are shown, since this adds similar types of clutter
	 */
	if (frame == frameCurrent) {
		gl_PointSize = float(pointSize + 7);
		finalColor = colorCurrentFrame;
	}
}
