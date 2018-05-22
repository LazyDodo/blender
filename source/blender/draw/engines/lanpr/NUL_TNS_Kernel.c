#include <gl/glew.h>
#include "NUL4.h"

/*

Ported from NUL4.0

Author(s):WuYiming - xp8110@outlook.com

*/

#include <math.h>


char TNS_VERTEX_SHADER_SRC_ATLAS_PREVIEW_130[] = "#version 330\n"

"\n"
"in vec4 vVertex;\n"
"in vec4 vColor;\n"
"uniform float uValue0; // buffer_w\n"
"out vec4 gColor;\n"
"\n"
"void main(){\n"
"    gl_Position = ((vVertex+1)/2*uValue0);\n"
"    gColor = vColor;\n"
"}";

char TNS_VERTEX_SHADER_SRC_UI_DEFALUT_130[] = "#version 330\n"
"\n"
"uniform mat4 mProjection;\n"
"uniform mat4 mModel;\n"
"uniform mat4 mView;\n"
"\n"
"in vec4 vVertex;\n"
"in vec4 vColor;\n"
"out vec4 fColor;\n"
"\n"
"void main(){\n"
"    gl_Position = mProjection * mView * mModel * vVertex;\n"
"    fColor = vColor;\n"
"}";

char TNS_VERTEX_SHADER_SRC_LINE_CONNECTION_130[] = "#version 330\n"
"\n"
"uniform mat4 mProjection;\n"
"uniform mat4 mModel;\n"
"uniform mat4 mView;\n"
"\n"
"in vec4 vVertex;\n"
"in vec4 vColor;\n"
"in float vUV;\n"
"out vec4 gColor;\n"
"out float gOffset;\n"
"\n"
"void main(){\n"
"    gl_Position = mProjection * mView * mModel * vVertex;\n"
"    gColor = vColor;\n"
"    gOffset = vUV;\n"
"}";

char TNS_FRAGMENT_SHADER_SRC_UI_DEFALUT_130[] = "#version 330\n"
"\n"
"in vec4 fColor;\n"
"\n"
"void main(){\n"
"    gl_FragColor = fColor;\n"
"}";

char TNS_VERTEX_SHADER_SRC_UI_TEXT_130[] = "#version 130\n"
"uniform mat4 mProjection;\n"
"uniform mat4 mModel;\n"
"uniform mat4 mView;\n"
"\n"
"in vec4 vVertex;\n"
"in vec4 vColor;\n"
"in vec2 vUV;\n"
"out vec4 fColor;\n"
"out vec2 fUV;\n"
"\n"
"void main(){\n"
"    gl_Position = mProjection * mView  * mModel * vVertex;\n"
"    fUV = vUV;\n"
"    fColor = vColor;\n"
"}";

char TNS_FRAGMENT_SHADER_SRC_UI_TEXT_130[] = "#version 130\n"
"\n"
"uniform sampler2D TexSample0;\n"
"in vec4 fColor;\n"
"in vec2 fUV;\n"
"\n"
"void main(){\n"
"    gl_FragColor = vec4(fColor.rgb, pow(texture2D(TexSample0,fUV.st).a,1/2.2));\n"
"}";
char TNS_FRAGMENT_SHADER_SRC_SINGLE_TEXTURE_130[] = "#version 130\n"
"\n"
"uniform sampler2D TexSample0;\n"
"in vec4 fColor;\n"
"in vec2 fUV;\n"
"\n"
"void main(){"
"    vec4 Color = texture2D(TexSample0,fUV.st);\n"
"    Color.a = pow(Color.a,1/2.2);\n"
"    gl_FragColor = Color;"
"}";
char TNS_FRAGMENT_SHADER_SRC_RECTANGLE_TEXTURE_130[] = "#version 130\n"
"\n"
"uniform sampler2DRect TexSample0;\n"
"uniform vec2 uValue0;"
"in vec4 fColor;\n"
"in vec2 fUV;\n"
"\n"
"void main(){"
"    vec4 Color = texture2DRect(TexSample0,vec2(fUV.s*uValue0.x,fUV.t*uValue0.y));\n"
"    Color.a = pow(Color.a,1/2.2);\n"
"    gl_FragColor = Color;"
"}";
char TNS_FRAGMENT_SHADER_SRC_TEXTURE_MULTIPLY_130[] = "#version 130\n"
"\n"
"uniform sampler2D TexSample0;\n"
"in vec4 fColor;\n"
"in vec2 fUV;\n"
"\n"
"void main(){"
"    vec4 Color = texture2D(TexSample0,fUV.st);\n"
"    Color*=fColor;\n"
//"    Color.a = pow(Color.a,1/2.2);\n"
"    gl_FragColor = Color;"
"}";
char TNS_FRAGMENT_SHADER_SRC_SINGLE_TEXTURE_MS_130[] = "#version 150\n"
"\n"
"uniform sampler2DMS TexSample0;\n"
"in vec4 fColor;\n"
"in vec2 fUV;\n"
"\n"
"void main(){"
"    vec4 Color= vec4(0.0);"
"    ivec2 texSize = textureSize(TexSample0);"
"    for(int i=0;i<4;i++)"
"        Color+=texelFetch(TexSample0, ivec2(fUV * texSize),i);\n"
"    gl_FragColor = Color/4;\n"
"}";
char TNS_FRAGMENT_SHADER_SRC_SINGLE_TEXTURE_MS_ALPHA_CTRL_130[] = "#version 150\n"
"\n"
"uniform sampler2DMS TexSample0;\n"
"in vec4 fColor;\n"
"in vec2 fUV;\n"
"\n"
"void main(){"
"    vec4 Color= vec4(0.0);"
"    ivec2 texSize = textureSize(TexSample0);"
"    for(int i=0;i<4;i++)"
"        Color+=texelFetch(TexSample0, ivec2(fUV * texSize),i);\n"
"    Color.a *= fColor.a;\n"
"    gl_FragColor = Color/4;\n"
"}";
char TNS_FRAGMENT_SHADER_SRC_SINGLE_TEXTURE_SS_130[] = "#version 150\n"
"\n"
"uniform sampler2DMS TexSample0;\n"
"in vec4 fColor;\n"
"in vec2 fUV;\n"
"\n"
"void main(){"
"    vec4 Color= vec4(0.0);"
"    ivec2 texSize = textureSize(TexSample0);"
"    for(int i=0;i<16;i++)"
"        Color+=texelFetch(TexSample0, ivec2(fUV * texSize),i);\n"
"    gl_FragColor = Color/16;\n"
"}";
char TNS_VERTEX_SHADER_SRC_SIMPLE_MATCAP_130[] = "#version 130\n"
"uniform mat4 mProjection;\n"
"uniform mat4 mModel;\n"
"uniform mat4 mView;\n"
"\n"
"in vec4 vVertex;\n"
"in vec3 vNormal;\n"
"smooth out vec3 fNormal;\n"
"\n"
"void main(){\n"
"    gl_Position = mProjection * mView  * mModel * vVertex;\n"
"    vec3 N = ( mView * mModel * vec4(vNormal,0)).xyz;\n"
"    fNormal = normalize(N);\n"
"}";
char TNS_FRAGMENT_SHADER_SRC_SIMPLE_MATCAP_130[] = "#version 130\n"
"\n"
"smooth in vec3 fNormal;"
"\n"
"float Interpolate(float between1,float between2,float value1,float value2,float key){\n"
"    float i = (key-between1)/(between2-between1);\n"
"    return value1*(1-i)+value2*i;"
"}\n"
"void main(){\n"
"    float value = dot(vec3(0,0,1),fNormal);\n"
"    if(value<0.65) value=0.15;\n"
"    else if(value>=0.65 && value<0.85) value=Interpolate(0.65,0.85,0.15,0.75,value);\n"
"    else if(value>=0.85 && value<0.95) value=0.75;\n"
"    else if(value>=0.95) value=0.9;\n"
"    gl_FragColor = vec4(vec3(0.84, 0.41, 0.16)*value,1);\n"
"}";


char TNS_VERTEX_SHADER_SRC_GRID_130[] = "#version 130\n"
"\n"
"uniform mat4 mProjection;\n"
"uniform mat4 mModel;\n"
"uniform mat4 mView;\n"
"\n"
"in vec4 vVertex;\n"
"in vec4 vColor;\n"
"in vec2 vUV;\n"
"out vec4 fColor;\n"
"out vec2 uv;\n"
"\n"
"void main(){\n"
"    vec4 pos = mProjection * mView * mModel * vVertex;\n"
"    gl_Position = pos;\n"
"    fColor = vColor;\n"
"    uv = vUV;\n"
"}";

char TNS_FRAGMENT_SHADER_SRC_TRANSPARNT_GRID_130[] = "#version 130\n"
"\n"
"in vec4 fColor;\n"
"in vec2 uv;\n"
"\n"
"void main(){\n"
"    vec4 c = fColor;\n"
"    c.a = sin(uv.x)*sin(uv.y)>0?c.a:0;\n"
"    gl_FragColor = c;\n"
"}";




char TNS_VERTEX_SHADER_SRC_EXTRA_BUFFERS_130[] = "#version 130\n"
"uniform mat4 mProjection;\n"
"uniform mat4 mModel;\n"
"uniform mat4 mView;\n"
"\n"
"in vec4 vVertex;\n"
"in vec3 vNormal;\n"
"in vec4 vColor;\n"
"smooth out vec3 fNormal;\n"
"smooth out vec4 fColor;\n"
"\n"
"void main(){\n"
"    gl_Position = mProjection * mView  * mModel * vVertex;\n"
"    vec3 N = ( mView * mModel * vec4(vNormal,0)).xyz;\n"
"    fNormal = normalize(N);\n"
"    fColor = vColor;\n"
"}";

char TNS_FRAGMENT_SHADER_SRC_EXTRA_BUFFERS_330[] = "#version 130\n"
"\n"
"smooth in vec3 fNormal;"
"smooth in vec4 fColor;"
"\n"
"float Interpolate(float between1,float between2,float value1,float value2,float key){\n"
"    float i = (key-between1)/(between2-between1);\n"
"    return value1*(1-i)+value2*i;"
"}\n"
"void main(){\n"
"    float value = dot(vec3(0,0,1),fNormal);\n"
//"    if(value<0.65) value=0.15;\n"
//"    else if(value>=0.65 && value<0.85) value=Interpolate(0.65,0.85,0.15,0.75,value);\n"
//"    else if(value>=0.85 && value<0.95) value=0.75;\n"
//"    else if(value>=0.95) value=0.9;\n"
"    gl_FragData[0] = vec4(vec3(fColor)*value,1);\n"
"    gl_FragData[1] = vec4((fNormal+vec3(1))*0.5,1);\n"
"}";

const real _ICON_SYMBOL_SIZE[14] = {
    12,
	2,
	2,
	4,
	2,
	2,
	2,
	4,
	2,
	2,
	2,
	2,
	2,
	1.5
};

int tKnlAttatchShader(tnsShader* ts);
tnsShader* tKnlFindShader1i(int CustomIndex);
void tKnlPushMatrix();
void tKnlPopMatrix();

nListHandle* tKnlGetTextureList();

void tMatLoadIdentity44d(tnsMatrix44d m);
void tMatMakeOrthographicMatrix44d(tnsMatrix44d mProjection, real xMin, real xMax, real yMin, real yMax, real zMin, real zMax);
void tMatMakePerspectiveMatrix44d(tnsMatrix44d mProjection, real fFov_rad, real fAspect, real zMin, real zMax);
void tMatMakeTranslationMatrix44d(tnsMatrix44d mTrans, real x, real y, real z);
void tMatMakeRotationMatrix44d(tnsMatrix44d m, real angle_rad, real x, real y, real z);
void tMatMakeScaleMatrix44d(tnsMatrix44d m, real x, real y, real z);
void tMatMultiply44d(tnsMatrix44d result, tnsMatrix44d l, tnsMatrix44d r);



void tMatInitFirstLevel(tnsMatrixStack* tms);
tnsMatrixStackItem* tKnlGetCurrentMatStackItem();
tnsShader* tKnlGetActiveShader();



void tnsAttach2DOffscreenBuffer(tnsOffscreen* target, GLuint attatchment, tnsTexture* use);
void tnsDetach2DOffscreenBuffer(tnsOffscreen* target, GLuint which_attach_point);

//=========================================================================================



tnsMain* T;
real DefaultZRange[2] = { 0.1,1000 };


//=======================================[SYS]
extern NUL MAIN;

void InitGLRenderEnviornment(){
	glEnable(GL_SCISSOR_TEST);
};
void tnsSwitchToCurrentWindowContext(nWindow* wnd){
	HGLRC current = wglGetCurrentContext(), hglrc = wnd->SystemGLRC;
	HDC hdc = wnd->SystemDC;

	//if (hglrc != current)
	int a = GetLastError();

		int s = wglMakeCurrent(hdc, hglrc);
};
void tnsViewportWithScissor(int x, int y, int w, int h){
	tnsShader* current_shader=0;
	glEnable(GL_SCISSOR_TEST);
	glViewport(x, y, w, h);
	glScissor(x, y, w, h);
	if ((current_shader = tKnlGetActiveShader())) {
		tnsResetViewMatrix();
		tSdrApplyView(current_shader, tnsGetViewMatrix());
		tnsResetModelMatrix();
		tSdrApplyModel(current_shader, tnsGetModelMatrix());
	}
}
void tnsViewport(int x, int y, int w, int h){
	glViewport(x, y, w, h);
}
//Must At Origion!

//void UseWindowMatrix_WithScissor(nWndPlacement* wp){
//	tnsViewportWithScissor(0, 0, VAL2_CLIENT_W_H(*wp));
//	SetMatrix2Di(0, VAL2_CLIENT_W_H(*wp), 0);
//};
//void UseBlockMatrix_WithScissor(nBlockPlacement* bp){
//	tnsViewportWithScissor(bp->UpperLeftX,
//		bp->ClientHeight - bp->LowerRightY,
//		bp->LowerRightX - bp->UpperLeftX,
//		bp->LowerRightY - bp->UpperLeftY);
//	SetMatrix2Di(0, bp->LowerRightX - bp->UpperLeftX,
//		0, bp->UpperLeftY - bp->LowerRightY);
//};
void ApplyEulerRotation(real x, real y, real z){
	glRotatef(z, 0, 0, 1);
	glRotatef(y, 0, 1, 0);
	glRotatef(x, 1, 0, 0);
}
void UseColor3d(real r, real g, real b){
	glColor3d(r, g, b);
};
void ClearGLScreen(int d, real r, real g, real b){
	glClearColor(r, g, b, 1);
	glClear(d);
};
void DrawWireRect2dp(real x, real y, real x2, real y2){
	real Verts[8];
	tnsMakeQuad2d(Verts, x, y, x2, y, x2, y2, x, y2);
	tnsVertexArray2d(Verts, 4);
	tnsPackAs(GL_LINE_LOOP);
};
void DrawWireRect4ds(real x, real y, real w, real h){
	real Verts[8];
	tnsMakeQuad2d(Verts, x, y, x+w, y, x+w, y+h, x, y+h);
	tnsVertexArray2d(Verts, 4);
	tnsPackAs(GL_LINE_LOOP);
};
void DrawWireRect4drs(real x, real y, real w, real h){
	real Verts[8];
	tnsMakeQuad2d(Verts, x, y, x+w, y, x+w, y-h, x, y-h);
	tnsVertexArray2d(Verts, 4);
	tnsPackAs(GL_LINE_LOOP);
};
void DrawSoildRect4drs(real x, real y, real w, real h){
	real Verts[8];
	tnsMakeQuad2d(Verts, x, y, x+w, y, x+w, y-h, x, y-h);
	tnsVertexArray2d(Verts, 4);
	tnsPackAs(GL_QUADS);
}
void DrawWireCross4ds(real x, real y, real w, real h){
	real Verts[8];
	tnsMakeQuad2d(Verts, x, y, x + w, y+h, x + w, y, x, y+h);
	tnsVertexArray2d(Verts, 4);
	tnsPackAs(GL_LINES);
};
void DrawWireCross4drs(real x, real y, real w, real h){
	real Verts[8];
	tnsMakeQuad2d(Verts, x, y, x + w, y - h, x + w, y, x, y - h);
	tnsVertexArray2d(Verts, 4);
	tnsPackAs(GL_LINES);
};
void DrawGrid6f(real L, real R, real U, real B, real stepX, real stepY){
	//do nothing;
}
void DrawSolidArrow1i3d(int direction, real xoff, real yoff, real size){
	real size_2 = size / 2.0f;
	real Verts[6];

	switch (direction){
	case TNS_ARROW_L:
		tnsMakeTriangle(Verts,
			xoff + size_2, yoff,
			xoff, yoff - size_2,
			xoff + size_2, yoff - size);
		break;
	case TNS_ARROW_R:
		tnsMakeTriangle(Verts,
			xoff, yoff,
			xoff + size_2, yoff - size_2,
			xoff, yoff - size);
		break;
	}

	tnsPackAs(GL_TRIANGLES);
}


//=======================================[Shader]

int tnsNextPowOf2(int i) {
	int result = 2;
	while (result < i) {
		result *= 2;
	}
	return result;
}

void tSdrMakeIndex(tnsShader* ts){
	int program;

	if (!ts)
		return;

	program = ts->glProgramID;

	ts->modelIndex = glGetUniformLocation(program, "mModel");
	ts->projectionIndex = glGetUniformLocation(program, "mProjection");
	ts->projectionInverseIndex = glGetUniformLocation(program, "mProjectionInverse");
	ts->viewIndex = glGetUniformLocation(program, "mView");
	ts->normalIndex = glGetUniformLocation(program, "mNormalScaler");

	ts->vertexIndex = glGetAttribLocation(program, "vVertex");
	ts->normalIndex = glGetAttribLocation(program, "vNormal");
	ts->colorIndex = glGetAttribLocation(program, "vColor");
	ts->uvIndex = glGetAttribLocation(program, "vUV");

	ts->texture0Index = glGetUniformLocation(program, "TexSample0");
	ts->texture1Index = glGetUniformLocation(program, "TexSample1");
	ts->texture2Index = glGetUniformLocation(program, "TexSample2");
	ts->texture3Index = glGetUniformLocation(program, "TexSample3");
	ts->texture4Index = glGetUniformLocation(program, "TexSample4");

	ts->uniform0Index = glGetUniformLocation(program, "uValue0");
	ts->uniform1Index = glGetUniformLocation(program, "uValue1");
	ts->uniform2Index = glGetUniformLocation(program, "uValue2");
	ts->uniform3Index = glGetUniformLocation(program, "uValue3");
	ts->uniform4Index = glGetUniformLocation(program, "uValue4");
	ts->uniform5Index = glGetUniformLocation(program, "uValue5");
	ts->uniform6Index = glGetUniformLocation(program, "uValue6");
	ts->uniform7Index = glGetUniformLocation(program, "uValue7");
}
void tSdrApplyProjection(tnsShader* ts, tnsMatrix44d m){
	tnsMatrix44f mf;
	if (T->BindedShader != ts) {
		glUseProgram(ts->glProgramID);
		T->BindedShader = ts;
	}
	tMatConvert44df(m, mf);
	glUniformMatrix4fv(ts->projectionIndex, 1, 0, mf);
}
void tMatInverse44d(tnsMatrix44d inverse, tnsMatrix44d mat);
void tSdrApplyProjectionInverse(tnsShader* ts, tnsMatrix44d m) {
	tnsMatrix44f mf;
	tnsMatrix44d i;
	if (T->BindedShader != ts) {
		glUseProgram(ts->glProgramID);
		T->BindedShader = ts;
	}
	tMatInverse44d(i, m);
	tMatConvert44df(i, mf);
	glUniformMatrix4fv(ts->projectionInverseIndex, 1, 0, mf);
}
void tSdrApplyModel(tnsShader* ts, tnsMatrix44d m){
	tnsMatrix44f mf;
	if (T->BindedShader != ts) {
		glUseProgram(ts->glProgramID);
		T->BindedShader = ts;
	}
	tMatConvert44df(m, mf);
	glUniformMatrix4fv(ts->modelIndex, 1, 0, mf);
}
void tSdrApplyView(tnsShader* ts, tnsMatrix44d m){
	tnsMatrix44f mf;
	if (T->BindedShader != ts) {
		glUseProgram(ts->glProgramID);
		T->BindedShader = ts;
	}
	tMatConvert44df(m, mf);
	glUniformMatrix4fv(ts->viewIndex, 1, 0, mf);
}
void tSdrApplyNormalScaler(tnsShader* ts, tnsMatrix44d m){
	tnsMatrix44f mf;
	if (ts->normalIndex == -1) return;
	if (T->BindedShader != ts) {
		glUseProgram(ts->glProgramID);
		T->BindedShader = ts;
	}
	tMatConvert44df(m, mf);
	glUniformMatrix4fv(ts->normalIndex, 1, 0, mf);
}

//--------------------Export
int tnsNewVertexShaderStringBased(char* Content){
	int status = 0;
	char error[1024];
	GLuint VertexShaderObject;
	tnsShader* s = 0;

	if (!Content)
		return -1;

	VertexShaderObject = glCreateShader(GL_VERTEX_SHADER);

	glShaderSource(VertexShaderObject, 1, &Content, 0);
	glCompileShader(VertexShaderObject);
	glGetShaderiv(VertexShaderObject, GL_COMPILE_STATUS, &status);
	if (status == GL_FALSE){
		glGetShaderInfoLog(VertexShaderObject, sizeof(error), 0, error);
		printf("Vertex shader error:\n%s", error);
		glDeleteShader(VertexShaderObject);
		return -1;
	}
	else{
		printf("Vertex OK!\n");
	}

	return VertexShaderObject;
}
int tnsNewFragmentShaderStringBased(char* Content){
	int status = 0;
	char error[1024];
	GLuint FragmentShaderObject;
	tnsShader* s = 0;

	if (!Content)
		return -1;

	FragmentShaderObject = glCreateShader(GL_FRAGMENT_SHADER);

	glShaderSource(FragmentShaderObject, 1, &Content, 0);
	glCompileShader(FragmentShaderObject);
	glGetShaderiv(FragmentShaderObject, GL_COMPILE_STATUS, &status);
	if (status == GL_FALSE){
		glGetShaderInfoLog(FragmentShaderObject, sizeof(error), 0, error);
		printf("Fragment shader error:\n%s", error);
		glDeleteShader(FragmentShaderObject);
		return -1;
	}
	else{
		printf("Fragment OK!\n");
	}

	return FragmentShaderObject;
}
int tnsNewGeometryShaderStringBased(char* Content) {
	int status = 0;
	char error[1024];
	GLuint GeometryShaderObject;
	tnsShader* s = 0;

	if (!Content)
		return -1;

	GeometryShaderObject = glCreateShader(GL_GEOMETRY_SHADER);

	glShaderSource(GeometryShaderObject, 1, &Content, 0);
	glCompileShader(GeometryShaderObject);
	glGetShaderiv(GeometryShaderObject, GL_COMPILE_STATUS, &status);
	if (status == GL_FALSE) {
		glGetShaderInfoLog(GeometryShaderObject, sizeof(error), 0, error);
		printf("Geometry shader error:\n%s", error);
		glDeleteShader(GeometryShaderObject);
		return -1;
	}
	else {
		printf("Geometry OK!\n");
	}

	return GeometryShaderObject;
}

tnsShader* tnsNewShaderProgram(int VertexShaderID, int FragmentShaderID,int GeometryShaderID){
	int vso = VertexShaderID;
	int fso = FragmentShaderID;
	int gso = GeometryShaderID;
	tnsShader* ts = 0;
	int status = 0;
	char error[1024];

	if (!vso || !fso)
		return 0;

	ts = CreateNew(tnsShader);
	ts->vtShaderID = vso;
	ts->fgShaderID = fso;
	ts->glProgramID = glCreateProgram();
	glAttachShader(ts->glProgramID, vso);
	glAttachShader(ts->glProgramID, fso);
	if(GeometryShaderID>-1)glAttachShader(ts->glProgramID, gso);
	glLinkProgram(ts->glProgramID); int a = glGetError();
	glGetProgramiv(ts->glProgramID, GL_LINK_STATUS, &status);
	if (status == GL_FALSE){
		glGetProgramInfoLog(ts->glProgramID, sizeof(error), 0, error);
		printf("Shader Linking error:\n%s", error);
		glDetachShader(ts->glProgramID, vso);
		glDetachShader(ts->glProgramID, fso);
		glDeleteShader(vso);
		glDeleteShader(fso);
		glDeleteProgram(ts->glProgramID);
		FreeMem(ts);
		return 0;
	}
	glUseProgram(ts->glProgramID);

	tSdrMakeIndex(ts);

	return ts;
}
int tnsEnableShader(int index){
	tnsMatrixStackItem* tmsi;
	tnsShader* ts = tKnlFindShader1i(index);
	if (!ts){
		glUseProgram(0);
		T->CurrentShader = 0;
		T->BindedShader = 0;
		return 0;
	}
	glUseProgram(ts->glProgramID);
	T->CurrentShader = ts;
	T->BindedShader = ts;

	tmsi = tKnlGetCurrentMatStackItem();
	tSdrApplyProjection(ts, tmsi->projection);
	tSdrApplyView(ts, tmsi->view);
	tSdrApplyModel(ts, tmsi->model);
	return 1;
}
int tnsEnableShaderv(tnsShader* shader){
	tnsMatrixStackItem* tmsi;
	tnsShader* ts = shader;
	if (!ts){
		glUseProgram(0);
		T->CurrentShader = 0;
		T->BindedShader = 0;
		return 0;
	}
	glUseProgram(ts->glProgramID);
	T->CurrentShader = ts;
	T->BindedShader = ts;

	tmsi = tKnlGetCurrentMatStackItem();
	tSdrApplyProjection(ts, tmsi->projection);
	tSdrApplyProjectionInverse(ts, tmsi->projection);
	tSdrApplyView(ts, tmsi->view);
	tSdrApplyModel(ts, tmsi->model);

	if (ts->vertexIndex != -1) glEnableVertexAttribArray(ts->vertexIndex);
	if (ts->colorIndex != -1) glEnableVertexAttribArray(ts->colorIndex);
	if (ts->normalIndex != -1)glEnableVertexAttribArray(ts->normalIndex);
	if (ts->uvIndex != -1) glEnableVertexAttribArray(ts->uvIndex);

	return 1;
}
int tnsUseShader(tnsShader* shader){
	tnsCommand* c = T->UsingCommand;
	c->ReplaceShader = shader;
}
void tnsUseUiShader(){
	T->UsingCommand->ReplaceShader = T->uiShader;
}
void tnsUseStringShader(){
	//tnsEnableShaderv(T->stringShader);
	T->UsingCommand->ReplaceShader = T->stringShader;
}
void tnsUseTextureShader(){
	//tnsEnableShaderv(T->TextureShader);
	T->UsingCommand->ReplaceShader = T->TextureShader;
}
void tnsUseMultiplyTextureShader() {
	T->UsingCommand->ReplaceShader = T->TextureMultiplyShader;
}
void tnsUseTextureMultisampleShader() {
	//tnsEnableShaderv(T->TextureShader);
	T->UsingCommand->ReplaceShader = T->MSTextureShader;
}
void tnsUseTextureMultisampleAlphaShader() {
	T->UsingCommand->ReplaceShader = T->MSATextureShader;
}
void tnsUseTransparentGridShader() {
	T->UsingCommand->ReplaceShader = T->TransparentGridShader;
}
void tnsUseTextureSobelColorMSShader() {
	T->UsingCommand->ReplaceShader = T->SobelColorShader;
}
void tnsUseTextureSobelMSShader() {
	T->UsingCommand->ReplaceShader = T->SobelShader;
}
void tnsUseExtraBufferShader() {
	T->UsingCommand->ReplaceShader = T->ExtraBuffersShader;
}
void tnsUseFixedShader(){
	tnsMatrixStackItem* tmsi = tKnlGetCurrentMatStackItem();
	glUseProgram(0);
	T->BindedShader = 0;
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glMultMatrixf(tmsi->model);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glMultMatrixf(tmsi->projection);
}
void tnsUseImagePeelShader() {
	//tnsEnableShaderv(T->TextureShader);
	T->UsingCommand->ReplaceShader = T->ImagePeelShader;
}
void tnsBindActiveTexture(){
	if (T && T->CurrentShader)
		glUniform1i(T->CurrentShader->texture0Index, 0);
};

//=======================================[MAT]


void tMatLoadIdentity44d(tnsMatrix44d m) {
	memset(m, 0, sizeof(tnsMatrix44d));
	m[0] = 1.0f;
	m[5] = 1.0f;
	m[10] = 1.0f;
	m[15] = 1.0f;
};

real tMatDistIdv2(real x1, real y1, real x2, real y2) {
	real x = x2 - x1, y = y2 - y1;
	return sqrt((x*x + y*y));
}
real tMatDist3dv(tnsVector3d l, tnsVector3d r) {
	real x = r[0] - l[0];
	real y = r[1] - l[1];
	real z = r[2] - l[2];
	return sqrt(x*x+y*y+z*z);
}
real tMatDist2dv(tnsVector2d l, tnsVector2d r) {
	real x = r[0] - l[0];
	real y = r[1] - l[1];
	return sqrt(x*x + y*y);
}

real tMatLength3d(tnsVector3d l) {
	return (sqrt(l[0] * l[0] + l[1] * l[1] + l[2] * l[2]));
}
real tMatLength2d(tnsVector3d l) {
	return (sqrt(l[0] * l[0] + l[1] * l[1]));
}
void tMatNormalize3d(tnsVector3d result, tnsVector3d l) {
	real r = sqrt(l[0] * l[0] + l[1] * l[1] + l[2] * l[2]);
	result[0] = l[0] / r;
	result[1] = l[1] / r;
	result[2] = l[2] / r;
}
void tMatNormalize2d(tnsVector2d result, tnsVector2d l) {
	real r = sqrt(l[0] * l[0] + l[1] * l[1]);
	result[0] = l[0] / r;
	result[1] = l[1] / r;
}
void tMatNormalizeSelf3d(tnsVector3d result) {
	real r = sqrt(result[0] * result[0] + result[1] * result[1] + result[2] * result[2]);
	result[0] /= r;
	result[1] /= r;
	result[2] /= r;
}
real tMatDot3d(tnsVector3d l, tnsVector3d r, int normalize) {
	tnsVector3d ln, rn;
	if (normalize) {
		tMatNormalize3d(ln, l); tMatNormalize3d(rn, r);
		return (ln[0] * rn[0] + ln[1] * rn[1] + ln[2] * rn[2]);
	}
	return (l[0] * r[0] + l[1] * r[1] + l[2] * r[2]);
}
real tMatDot2d(tnsVector2d l, tnsVector2d r, int normalize) {
	tnsVector3d ln, rn;
	if (normalize) {
		tMatNormalize2d(ln, l); tMatNormalize2d(rn, r);
		return (ln[0] * rn[0] + ln[1] * rn[1]);
	}
	return (l[0] * r[0] + l[1] * r[1]);
}
real tMatVectorCross3d(tnsVector3d result, tnsVector3d l, tnsVector3d r) {
	result[0] = l[1] * r[2] - l[2] * r[1];
	result[1] = l[2] * r[0] - l[0] * r[2];
	result[2] = l[0] * r[1] - l[1] * r[0];
	return tMatLength3d(result);
}
void tMatVectorCrossOnly3d(tnsVector3d result, tnsVector3d l, tnsVector3d r) {
	result[0] = l[1] * r[2] - l[2] * r[1];
	result[1] = l[2] * r[0] - l[0] * r[2];
	result[2] = l[0] * r[1] - l[1] * r[0];
}
real tMatAngleRad3d(tnsVector3d from, tnsVector3d to, tnsVector3d PositiveReference) {
	if (PositiveReference) {
		tnsVector3d res;
		tMatVectorCross3d(res, from, to);
		if (tMatDot3d(res, PositiveReference, 1)>0)
			return acosf(tMatDot3d(from, to, 1));
		else
			return -acosf(tMatDot3d(from, to, 1));
	}
	return acosf(tMatDot3d(from, to, 1));
}
void tMatApplyRotation33d(tnsVector3d result, tnsMatrix44d mat, tnsVector3d v) {
	result[0] = mat[0] * v[0] + mat[1] * v[1] + mat[2] * v[2];
	result[1] = mat[3] * v[0] + mat[4] * v[1] + mat[5] * v[2];
	result[2] = mat[6] * v[0] + mat[7] * v[1] + mat[8] * v[2];
}
void tMatApplyRotation43d(tnsVector3d result, tnsMatrix44d mat, tnsVector3d v) {
	result[0] = mat[0] * v[0] + mat[1] * v[1] + mat[2] * v[2];
	result[1] = mat[4] * v[0] + mat[5] * v[1] + mat[6] * v[2];
	result[2] = mat[8] * v[0] + mat[9] * v[1] + mat[10] * v[2];
}
void tMatApplyTransform43d(tnsVector3d result, tnsMatrix44d mat, tnsVector3d v) {
	real w;
	result[0] = mat[0] * v[0] + mat[4] * v[1] + mat[8] * v[2] + mat[12] * 1;
	result[1] = mat[1] * v[0] + mat[5] * v[1] + mat[9] * v[2] + mat[13] * 1;
	result[2] = mat[2] * v[0] + mat[6] * v[1] + mat[10] * v[2] + mat[14] * 1;
	w = mat[3] * v[0] + mat[7] * v[1] + mat[11] * v[2] + mat[15] * 1;
	result[0] /= w;
	result[1] /= w;
	result[2] /= w;
}
void tMatApplyNormalTransform43d(tnsVector3d result, tnsMatrix44d mat, tnsVector3d v) {
	real w;
	result[0] = mat[0] * v[0] + mat[4] * v[1] + mat[8] * v[2] + mat[12] * 1;
	result[1] = mat[1] * v[0] + mat[5] * v[1] + mat[9] * v[2] + mat[13] * 1;
	result[2] = mat[2] * v[0] + mat[6] * v[1] + mat[10] * v[2] + mat[14] * 1;
}
void tMatApplyTransform44d(tnsVector4d result, tnsMatrix44d mat, tnsVector4d v) {
	result[0] = mat[0] * v[0] + mat[4] * v[1] + mat[8] * v[2] + mat[12] * 1;
	result[1] = mat[1] * v[0] + mat[5] * v[1] + mat[9] * v[2] + mat[13] * 1;
	result[2] = mat[2] * v[0] + mat[6] * v[1] + mat[10] * v[2] + mat[14] * 1;
	result[3] = mat[3] * v[0] + mat[7] * v[1] + mat[11] * v[2] + mat[15] * 1;
}
void tMatApplyTransform44dTrue(tnsVector4d result, tnsMatrix44d mat, tnsVector4d v) {
	result[0] = mat[0] * v[0] + mat[4] * v[1] + mat[8] * v[2] + mat[12] * v[3];
	result[1] = mat[1] * v[0] + mat[5] * v[1] + mat[9] * v[2] + mat[13] * v[3];
	result[2] = mat[2] * v[0] + mat[6] * v[1] + mat[10] * v[2] + mat[14] * v[3];
	result[3] = mat[3] * v[0] + mat[7] * v[1] + mat[11] * v[2] + mat[15] * v[3];
}

void tMatRemoveTranslation44d(tnsMatrix44d result, tnsMatrix44d mat) {
	tMatLoadIdentity44d(result);
	result[0] = mat[0];
	result[1] = mat[1];
	result[2] = mat[2];
	result[4] = mat[4];
	result[5] = mat[5];
	result[6] = mat[6];
	result[8] = mat[8];
	result[9] = mat[9];
	result[10] = mat[10];
}
void tMatClearTranslation44d(tnsMatrix44d mat) {
	mat[3] = 0;
	mat[7] = 0;
	mat[11] = 0;
}


void tMatExtractXYZEuler44d(tnsMatrix44d mat, real* x_result, real* y_result, real* z_result) {
	real xRot, yRot, zRot;

	if (mat[2] < 1) {
		if (mat[2] > -1) {
			yRot = asin(mat[2]);
			xRot = atan2(-mat[6], mat[10]);
			zRot = atan2(-mat[1], mat[0]);
		}
		else {
			yRot = -TNS_PI / 2;
			xRot = -atan2(-mat[4], mat[5]);
			zRot = 0;
		}
	}
	else {
		yRot = TNS_PI / 2;
		xRot = atan2(-mat[4], mat[5]);
		zRot = 0;
	}

	(*x_result) = -xRot;
	(*y_result) = -yRot;
	(*z_result) = -zRot;
}
void tMatExtractLocation44d(tnsMatrix44d mat, real* x_result, real* y_result, real* z_result) {
	*x_result = mat[12];
	*y_result = mat[13];
	*z_result = mat[14];
}
void tMatExtractUniformScale44d(tnsMatrix44d mat, real* result) {
	tnsVector3d v = { mat[0],mat[1],mat[2] };
	*result = tMatLength3d(v);
}



#define L(row,col)  l[(col<<2)+row]
#define R(row,col)  r[(col<<2)+row]
#define P(row,col)  result[(col<<2)+row]

void tMatPrintMatrix44d(tnsMatrix44d l) {
	int i, j;
	for (i = 0; i < 4; i++) {
		for (j = 0; j < 4; j++) {
			printf("%.5f ", L(i, j));
		}
		printf("\n");
	}
}

void tMatCopyMatrix44d(tnsMatrix44d from, tnsMatrix44d to) {
	memcpy(to, from, sizeof(tnsMatrix44d));
}
void tMatMultiply44d(tnsMatrix44d result, tnsMatrix44d l, tnsMatrix44d r) {
	int i;
	for (i = 0; i < 4; i++) {
		real ai0 = L(i, 0), ai1 = L(i, 1), ai2 = L(i, 2), ai3 = L(i, 3);
		P(i, 0) = ai0 * R(0, 0) + ai1 * R(1, 0) + ai2 * R(2, 0) + ai3 * R(3, 0);
		P(i, 1) = ai0 * R(0, 1) + ai1 * R(1, 1) + ai2 * R(2, 1) + ai3 * R(3, 1);
		P(i, 2) = ai0 * R(0, 2) + ai1 * R(1, 2) + ai2 * R(2, 2) + ai3 * R(3, 2);
		P(i, 3) = ai0 * R(0, 3) + ai1 * R(1, 3) + ai2 * R(2, 3) + ai3 * R(3, 3);
	}
};
void tMatInverse44d(tnsMatrix44d inverse, tnsMatrix44d mat){
	int i, j, k;
	double temp;
	tnsMatrix44d tempmat;
	real max;
	int maxj;

	tMatLoadIdentity44d(inverse);

	tMatCopyMatrix44d(mat, tempmat);

	for (i = 0; i < 4; i++) {
		/* Look for row with max pivot */
		max = fabsf(tempmat[i*5]);
		maxj = i;
		for (j = i + 1; j < 4; j++) {
			if (fabsf(tempmat[j*4 + i]) > max) {
				max = fabsf(tempmat[j * 4 + i]);
				maxj = j;
			}
		}

		/* Swap rows if necessary */
		if (maxj != i) {
			for (k = 0; k < 4; k++) {
				real t;
				t = tempmat[i * 4 + k];
				tempmat[i * 4 + k] = tempmat[maxj * 4 + k];
				tempmat[maxj * 4 + k] = t;

				t = inverse[i * 4 + k];
				inverse[i * 4 + k] = inverse[maxj * 4 + k];
				inverse[maxj * 4 + k] = t;
			}
		}

		//if (UNLIKELY(tempmat[i][i] == 0.0f)) {
		//	return false;  /* No non-zero pivot */
		//}

		temp = (double)tempmat[i*5];
		for (k = 0; k < 4; k++) {
			tempmat[i * 4 + k] = (real)((double)tempmat[i * 4 + k] / temp);
			inverse[i * 4 + k] = (real)((double)inverse[i * 4 + k] / temp);
		}
		for (j = 0; j < 4; j++) {
			if (j != i) {
				temp = tempmat[j*4 + i];
				for (k = 0; k < 4; k++) {
					tempmat[j * 4 + k] -= (real)((double)tempmat[i * 4 + k] * temp);
					inverse[j * 4 + k] -= (real)((double)inverse[i * 4 + k] * temp);
				}
			}
		}
	}
}
void tMatMakeTranslationMatrix44d(tnsMatrix44d mTrans, real x, real y, real z) {
	tMatLoadIdentity44d(mTrans);
	mTrans[12] = x;
	mTrans[13] = y;
	mTrans[14] = z;
}
void tMatMakePerspectiveMatrix44d(tnsMatrix44d mProjection, real fFov_rad, real fAspect, real zMin, real zMax) {
	real yMax = zMin * tanf(fFov_rad * 0.5f);
	real yMin = -yMax;
	real xMin = yMin * fAspect;
	real xMax = -xMin;

	tMatLoadIdentity44d(mProjection);

	mProjection[0] = (2.0f * zMin) / (xMax - xMin);
	mProjection[5] = (2.0f * zMin) / (yMax - yMin);
	mProjection[8] = (xMax + xMin) / (xMax - xMin);
	mProjection[9] = (yMax + yMin) / (yMax - yMin);
	mProjection[10] = -((zMax + zMin) / (zMax - zMin));
	mProjection[11] = -1.0f;
	mProjection[14] = -((2.0f * (zMax*zMin)) / (zMax - zMin));
	mProjection[15] = 0.0f;
}
void tMatMakeZTrackingMatrix44d(tnsMatrix44d mat, tnsVector3d this,tnsVector3d that, tnsVector3d up) {
	tnsVector4d fwd, l, t,rt;
	fwd[3] = l[3] = t[3] = rt[3] = 1;
	t[0] = up[0];
	t[1] = up[1];
	t[2] = up[2];
	fwd[0] = that[0] - this[0];
	fwd[1] = that[1] - this[1];
	fwd[2] = that[2] - this[2];

	tMatLoadIdentity44d(mat);
	
	tMatVectorCross3d(l, t, fwd);
	tMatVectorCross3d(rt, fwd, l);

	tMatNormalizeSelf3d(l);
	tMatNormalizeSelf3d(rt);
	tMatNormalizeSelf3d(fwd);

	mat[0] = l[0];
	mat[1] = l[1];
	mat[2] = l[2];

	mat[4] = rt[0];
	mat[5] = rt[1];
	mat[6] = rt[2];

	mat[8] = fwd[0];
	mat[9] = fwd[1];
	mat[10] = fwd[2];
}
void tMatMakeZTrackingMatrixDelta44d(tnsMatrix44d mat, tnsVector3d delta, tnsVector3d up) {
	tnsVector4d fwd, l, t, rt;
	fwd[3] = l[3] = t[3] = rt[3] = 1;
	t[0] = up[0];
	t[1] = up[1];
	t[2] = up[2];
	fwd[0] = delta[0];
	fwd[1] = delta[1];
	fwd[2] = delta[2];

	tMatLoadIdentity44d(mat);

	tMatVectorCross3d(l, t, fwd);
	tMatVectorCross3d(rt, fwd, l);

	tMatNormalizeSelf3d(l);
	tMatNormalizeSelf3d(rt);
	tMatNormalizeSelf3d(fwd);

	mat[0] = l[0];
	mat[1] = l[1];
	mat[2] = l[2];

	mat[4] = rt[0];
	mat[5] = rt[1];
	mat[6] = rt[2];

	mat[8] = fwd[0];
	mat[9] = fwd[1];
	mat[10] = fwd[2];
}
void tMatMakeOrthographicMatrix44d(tnsMatrix44d mProjection, real xMin, real xMax, real yMin, real yMax, real zMin, real zMax) {
	tMatLoadIdentity44d(mProjection);

	mProjection[0] = 2.0f / (xMax - xMin);
	mProjection[5] = 2.0f / (yMax - yMin);
	mProjection[10] = -2.0f / (zMax - zMin);
	mProjection[12] = -((xMax + xMin) / (xMax - xMin));
	mProjection[13] = -((yMax + yMin) / (yMax - yMin));
	mProjection[14] = -((zMax + zMin) / (zMax - zMin));
	mProjection[15] = 1.0f;
}
void tMatMakeRotationMatrix44d(tnsMatrix44d m, real angle_rad, real x, real y, real z)
{
	real mag, s, c;
	real xx, yy, zz, xy, yz, zx, xs, ys, zs, one_c;

	s = (real)sin(angle_rad);
	c = (real)cos(angle_rad);

	mag = (real)sqrt(x*x + y*y + z*z);

	// Identity matrix
	if (mag == 0.0f) {
		tMatLoadIdentity44d(m);
		return;
	}

	// Rotation matrix is normalized
	x /= mag;
	y /= mag;
	z /= mag;

#define M(row,col)  m[col*4+row]

	xx = x * x;
	yy = y * y;
	zz = z * z;
	xy = x * y;
	yz = y * z;
	zx = z * x;
	xs = x * s;
	ys = y * s;
	zs = z * s;
	one_c = 1.0f - c;

	M(0, 0) = (one_c * xx) + c;
	M(0, 1) = (one_c * xy) + zs;
	M(0, 2) = (one_c * zx) + ys;
	M(0, 3) = 0.0f;

	M(1, 0) = (one_c * xy) - zs;
	M(1, 1) = (one_c * yy) + c;
	M(1, 2) = (one_c * yz) + xs;
	M(1, 3) = 0.0f;

	M(2, 0) = (one_c * zx) + ys;
	M(2, 1) = (one_c * yz) - xs;
	M(2, 2) = (one_c * zz) + c;
	M(2, 3) = 0.0f;

	M(3, 0) = 0.0f;
	M(3, 1) = 0.0f;
	M(3, 2) = 0.0f;
	M(3, 3) = 1.0f;

#undef M
}
void tMatMakeRotationXMatrix44d(tnsMatrix44d m, real angle_rad) {
	tMatLoadIdentity44d(m);
	m[5] = cos(angle_rad);
	m[6] = sin(angle_rad);
	m[9] = -sin(angle_rad);
	m[10] = cos(angle_rad);
}
void tMatMakeRotationYMatrix44d(tnsMatrix44d m, real angle_rad) {
	tMatLoadIdentity44d(m);
	m[0] = cos(angle_rad);
	m[2] = -sin(angle_rad);
	m[8] = sin(angle_rad);
	m[10] = cos(angle_rad);
}
void tMatMakeRotationZMatrix44d(tnsMatrix44d m, real angle_rad) {
	tMatLoadIdentity44d(m);
	m[0] = cos(angle_rad);
	m[1] = sin(angle_rad);
	m[4] = -sin(angle_rad);
	m[5] = cos(angle_rad);
}
void tMatMakeScaleMatrix44d(tnsMatrix44d m, real x, real y, real z) {
	tMatLoadIdentity44d(m);
	m[0] = x;
	m[5] = y;
	m[10] = z;
}
void tMatMakeViewportMatrix44d(tnsMatrix44d m, real w, real h,real Far,real Near) {
	tMatLoadIdentity44d(m);
	m[0] = w / 2;
	m[5] = h / 2;
	m[10] = (Far-Near)/2;
	m[12] = w/2;
	m[13] = h/2;
	m[14] = (Far+Near)/2;
	m[15] = 1;
	//m[0] = 2/w;
	//m[5] = 2/h;
	//m[10] = 1;
	//m[12] = 2/w;
	//m[13] = 2/h;
	//m[14] = 1;
	//m[15] = 1;
}


void tMatInitFirstLevel(tnsMatrixStack* tms) {
	tMatLoadIdentity44d(tms->level[0].model);
	tMatLoadIdentity44d(tms->level[0].projection);
	tMatLoadIdentity44d(tms->level[0].view);
}

//-------------------Export

void tnsOrtho(real xMin, real xMax, real yMin, real yMax, real zMin, real zMax){
	tnsShader* current_shader = 0;
	real* mat = tnsGetProjectionMatrix();

	tMatMakeOrthographicMatrix44d(mat, xMin, xMax, yMin, yMax, zMin, zMax);

	if (current_shader = tKnlGetActiveShader()){
		tSdrApplyProjection(current_shader, mat);
	}
}
void tnsPerspective(real fFov_rad, real fAspect, real zMin, real zMax){
	tnsShader* current_shader = 0;
	real* mat = tnsGetProjectionMatrix();

	tMatMakePerspectiveMatrix44d(mat, fFov_rad, fAspect, zMin, zMax);

	if (current_shader = tKnlGetActiveShader()){
		tSdrApplyProjection(current_shader, mat);
	}
}

void tnsTranslate3d(real x, real y, real z){
	tnsShader* current_shader = 0;
	real* mat = tnsGetModelMatrix();
	tnsMatrix44d transmat, result;


	tMatMakeTranslationMatrix44d(transmat, x, y, z);
	tMatMultiply44d(result, mat, transmat);
	memcpy(mat, result, sizeof(tnsMatrix44d));

	if (current_shader = tKnlGetActiveShader()){
		tSdrApplyModel(current_shader, result);
	}
}
void tnsPreTranslate3d(real x, real y, real z){
	tnsShader* current_shader = 0;
	real* mat = tnsGetModelMatrix();
	tnsMatrix44d transmat, result;

	//glTranslatef(x, y, z);

	tMatMakeTranslationMatrix44d(transmat, x, y, z);
	tMatMultiply44d(result, mat, transmat);
	memcpy(mat, result, sizeof(tnsMatrix44d));
}
void tnsRotate4d(real degrees, real x, real y, real z){
	tnsShader* current_shader = 0;
	real* mat = tnsGetModelMatrix();
	tnsMatrix44d rotmat, result;

	glRotatef(degrees, x, y, z);

	tMatMakeRotationMatrix44d(rotmat, rad(degrees), x, y, z);
	tMatMultiply44d(result, mat, rotmat);
	memcpy(mat, result, sizeof(tnsMatrix44d));

	if (current_shader = tKnlGetActiveShader()){
		tSdrApplyModel(current_shader, result);
	}
}
void tnsPreRotate4d(real degrees, real x, real y, real z){
	tnsShader* current_shader = 0;
	real* mat = tnsGetModelMatrix();
	tnsMatrix44d rotmat, result;

	glRotatef(degrees, x, y, z);

	tMatMakeRotationMatrix44d(rotmat, rad(degrees), x, y, z);
	tMatMultiply44d(result, mat, rotmat);
	memcpy(mat, result, sizeof(tnsMatrix44d));
}
void tnsScale3d(real x, real y, real z){
	tnsShader* current_shader = 0;
	real* mat  = tnsGetModelMatrix();
	tnsMatrix44d scalemat, normal_scaler, result;

	glScalef(x, y, z);

	tMatMakeScaleMatrix44d(scalemat, x, y, z);

	tMatMultiply44d(result, mat, scalemat);
	memcpy(mat, result, sizeof(tnsMatrix44d));

	if (current_shader = tKnlGetActiveShader()){
		tSdrApplyModel(current_shader, result);
	}
}
void tnsPreScale3d(real x, real y, real z){
	tnsShader* current_shader = 0;
	real* mat = tnsGetModelMatrix();
	tnsMatrix44d scalemat, normal_scaler, result;

	glScalef(x, y, z);

	tMatMakeScaleMatrix44d(scalemat, x, y, z);

	tMatMultiply44d(result, mat, scalemat);
	memcpy(mat, result, sizeof(tnsMatrix44d));
}

void tKnlPopMatrix();
void tKnlPushMatrix();
void tnsPushMatrix(){
	tKnlPushMatrix();
};
void tnsPopMatrix(){
	tKnlPopMatrix();
}

//========================================[KNL]

void tnsInitRenderKernel(int matrixStackLevel){
	tnsCommand* c;
	T = HyperNew(tnsMain);
	GLuint m_nQuadVAO;

	nutMakeHyperData(&T->World);

	T->GLVersionStr = glGetString(GL_VERSION);
	T->GLVendorStr = glGetString(GL_VENDOR);
	T->GLRendererStr = glGetString(GL_RENDERER);
	T->GLSLVersionStr = glGetString(GL_SHADING_LANGUAGE_VERSION);

	T->stack.max_level = matrixStackLevel;
	T->stack.level = CreateNewBuffer(tnsMatrixStackItem, matrixStackLevel);
	T->NextShaderIndex = 1;
	T->BindedShader = -1;
	tMatInitFirstLevel(&T->stack);

	//Immediate-style capable:

	c = T->UsingCommand = T->DrawingCommand;
	c->ColorBegin = c->ColorEnd = T->NextColor = 0;
	c->NormalBegin = c->NormalEnd = T->NextNormal = 0;
	c->TexCoordBegin = c->TexCoordEnd = T->NextTexCoord = 0;
	c->VertBegin = c->VertEnd = T->NextVert = 0;
	c->IndexBegin = c->IndexEnd = T->NextIndex = 0;

	tnsGetGaussianKernel(&T->EdgeGaussFilter, TNS_SNAKE_FILTER_SIZE, 0.7);

	/*glGenVertexArrays(1,&m_nQuadVAO);
    glBindVertexArray(m_nQuadVAO);*/


	glGenBuffers(1,&T->VertBufObject);
	glBindBuffer(GL_ARRAY_BUFFER, T->VertBufObject);
	glBufferData(GL_ARRAY_BUFFER, 8192 * sizeof(GLfloat), 0, GL_DYNAMIC_DRAW);

	glGenBuffers(1, &T->ColorBufObject);
	glBindBuffer(GL_ARRAY_BUFFER, T->ColorBufObject);
	glBufferData(GL_ARRAY_BUFFER, 8192 * sizeof(GLfloat), 0, GL_DYNAMIC_DRAW);

	glGenBuffers(1, &T->NormalBufObject);
	glBindBuffer(GL_ARRAY_BUFFER, T->NormalBufObject);
	glBufferData(GL_ARRAY_BUFFER, 8192 * sizeof(GLfloat), 0, GL_DYNAMIC_DRAW);

	glGenBuffers(1, &T->TexCoordBufObject);
	glBindBuffer(GL_ARRAY_BUFFER, T->TexCoordBufObject);
	glBufferData(GL_ARRAY_BUFFER, 8192 * sizeof(GLfloat), 0, GL_DYNAMIC_DRAW);

	glGenBuffers(1, &T->IndexBufObject);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, T->IndexBufObject);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, 1024 * sizeof(GLushort), 0, GL_DYNAMIC_DRAW);
}
void tnsInitBuiltinShaders(){
	T->uiShader = tnsNewShaderProgram(
		tnsNewVertexShaderStringBased(TNS_VERTEX_SHADER_SRC_UI_DEFALUT_130),
		tnsNewFragmentShaderStringBased(TNS_FRAGMENT_SHADER_SRC_UI_DEFALUT_130), -1);

	T->stringShader = tnsNewShaderProgram(
		tnsNewVertexShaderStringBased(TNS_VERTEX_SHADER_SRC_UI_TEXT_130),
		tnsNewFragmentShaderStringBased(TNS_FRAGMENT_SHADER_SRC_UI_TEXT_130), -1);

	T->TextureShader = tnsNewShaderProgram(
		tnsNewVertexShaderStringBased(TNS_VERTEX_SHADER_SRC_UI_TEXT_130),
		tnsNewFragmentShaderStringBased(TNS_FRAGMENT_SHADER_SRC_SINGLE_TEXTURE_130), -1);

	T->TextureMultiplyShader = tnsNewShaderProgram(
		tnsNewVertexShaderStringBased(TNS_VERTEX_SHADER_SRC_UI_TEXT_130),
		tnsNewFragmentShaderStringBased(TNS_FRAGMENT_SHADER_SRC_TEXTURE_MULTIPLY_130), -1);

	T->RectangleTextureShader = tnsNewShaderProgram(
		tnsNewVertexShaderStringBased(TNS_VERTEX_SHADER_SRC_UI_TEXT_130),
		tnsNewFragmentShaderStringBased(TNS_FRAGMENT_SHADER_SRC_RECTANGLE_TEXTURE_130), -1);

	T->MSTextureShader = tnsNewShaderProgram(
		tnsNewVertexShaderStringBased(TNS_VERTEX_SHADER_SRC_UI_TEXT_130),
		tnsNewFragmentShaderStringBased(TNS_FRAGMENT_SHADER_SRC_SINGLE_TEXTURE_MS_130), - 1);

	T->MSATextureShader = tnsNewShaderProgram(
		tnsNewVertexShaderStringBased(TNS_VERTEX_SHADER_SRC_UI_TEXT_130),
		tnsNewFragmentShaderStringBased(TNS_FRAGMENT_SHADER_SRC_SINGLE_TEXTURE_MS_ALPHA_CTRL_130), -1);

	T->TEST_MatcapShader = tnsNewShaderProgram(
		tnsNewVertexShaderStringBased(TNS_VERTEX_SHADER_SRC_SIMPLE_MATCAP_130),
		tnsNewFragmentShaderStringBased(TNS_FRAGMENT_SHADER_SRC_SIMPLE_MATCAP_130), -1);

	T->TransparentGridShader = tnsNewShaderProgram(
		tnsNewVertexShaderStringBased(TNS_VERTEX_SHADER_SRC_GRID_130),
		tnsNewFragmentShaderStringBased(TNS_FRAGMENT_SHADER_SRC_TRANSPARNT_GRID_130), -1);

	T->SobelColorShader = tnsNewShaderProgram(
		tnsNewVertexShaderStringBased(TNS_VERTEX_SHADER_SRC_UI_TEXT_130),
		tnsNewFragmentShaderStringBased(txtReadFileAsString("../NUL4/TNS_SnakeEdge.fragment")), -1);

	T->ExtraBuffersShader = tnsNewShaderProgram(
		tnsNewVertexShaderStringBased(TNS_VERTEX_SHADER_SRC_EXTRA_BUFFERS_130),
		tnsNewFragmentShaderStringBased(TNS_FRAGMENT_SHADER_SRC_EXTRA_BUFFERS_330), -1);

	T->AtlasTransformShader = tnsNewShaderProgram(
		//tnsNewVertexShaderStringBased(TNS_VERTEX_SHADER_SRC_UI_DEFALUT_130),
		tnsNewVertexShaderStringBased(txtReadFileAsString("../NUL4/TNS_AtlasProjectPassThrough.vertex")),
		tnsNewFragmentShaderStringBased(txtReadFileAsString("../NUL4/TNS_AtlasProjectClip.fragment")), -1);

	T->AtlasPreviewShader = tnsNewShaderProgram(
		//tnsNewVertexShaderStringBased(TNS_VERTEX_SHADER_SRC_UI_DEFALUT_130),
		tnsNewVertexShaderStringBased(TNS_VERTEX_SHADER_SRC_ATLAS_PREVIEW_130),
		tnsNewFragmentShaderStringBased(TNS_FRAGMENT_SHADER_SRC_UI_DEFALUT_130),
		//-1);
		tnsNewGeometryShaderStringBased(txtReadFileAsString("../NUL4/TNS_AtlasPreview.geometry")));

	T->ImagePeelShader = tnsNewShaderProgram(
		tnsNewVertexShaderStringBased(TNS_VERTEX_SHADER_SRC_UI_TEXT_130),
		tnsNewFragmentShaderStringBased(txtReadFileAsString("../NUL4/TNS_ImagePeel.fragment")), -1);

	T->LineConnectionShader = tnsNewShaderProgram(
		tnsNewVertexShaderStringBased(TNS_VERTEX_SHADER_SRC_LINE_CONNECTION_130),
		tnsNewFragmentShaderStringBased(TNS_FRAGMENT_SHADER_SRC_UI_DEFALUT_130),
		//-1);
		tnsNewGeometryShaderStringBased(txtReadFileAsString("../NUL4/TNS_LineConnection.geometry")));



	tnsEnableShaderv(T->uiShader);
}
void tnsInitWindowDefaultRenderConfig(){
	tnsInitBuiltinShaders();
};
void tnsQuit(){
	FreeMem(T->stack.level);
	FreeMem(T);
}

tnsMatrixStackItem* tKnlGetCurrentMatStackItem(){
	return &T->stack.level[T->stack.current_level];
}
tnsShader* tKnlGetActiveShader(){
	return T->CurrentShader;
}
real* tnsGetModelMatrix() {
	return tKnlGetCurrentMatStackItem()->model;
}
real* tnsGetViewMatrix() {
	return tKnlGetCurrentMatStackItem()->view;
}
real* tnsGetProjectionMatrix() {
	return tKnlGetCurrentMatStackItem()->projection;
}
void tnsResetModelMatrix() {
	tMatLoadIdentity44d(tnsGetModelMatrix());
}
void tnsResetViewMatrix() {
	tMatLoadIdentity44d(tnsGetViewMatrix());
}
void tnsResetProjectionMatrix() {
	tMatLoadIdentity44d(tnsGetProjectionMatrix());
}
int tKnlAttatchShader(tnsShader* ts){
	lstAppendItem(&T->Shaders, ts);
	ts->CustomID = T->NextShaderIndex;
	T->NextShaderIndex++;
	return ts->CustomID;
};
int CMP_NUM_IsThisShader(tnsShader* ts, int* index){
	return (ts->CustomID == *index);
}
tnsShader* tKnlFindShader1i(int CustomIndex){
	return lstFindItem(&CustomIndex, CMP_NUM_IsThisShader, &T->Shaders);
}
void tKnlPushMatrix(){
	tnsMatrixStack* tms = &T->stack;
	tnsShader* current_shader = tKnlGetActiveShader();

	if (tms->current_level == tms->max_level || !current_shader)
		return;

	memcpy(&tms->level[tms->current_level + 1], &tms->level[tms->current_level], sizeof(tnsMatrixStackItem));

	tms->current_level++;

	//tSdrApplyModel(current_shader,tms->level[tms->current_level].model);
	//tSdrApplyView(current_shader,tms->level[tms->current_level].view);
	//tSdrApplyProjection(current_shader,tms->level[tms->current_level].projection);
};
void tKnlPopMatrix(){
	tnsMatrixStack* tms = &T->stack;
	tnsShader* current_shader = tKnlGetActiveShader();

	if (tms->current_level == 0 || !current_shader)
		return;

	tms->current_level--;

	tSdrApplyModel(current_shader, tms->level[tms->current_level].model);
	tSdrApplyView(current_shader, tms->level[tms->current_level].view);
	tSdrApplyProjection(current_shader, tms->level[tms->current_level].projection);

}
nListHandle* tKnlGetTextureList(){
	return &T->Textures;
}

tnsBatch* tnsCreateBatch(u32bit NumVert, int Dimention, float* Data, float* Normal) {
	tnsBatch* b = CreateNew(tnsBatch);
	b->Dimention = Dimention;
	b->NumVert = NumVert;
	
	glGenBuffers(1, &b->VBO);
	glBindBuffer(GL_ARRAY_BUFFER, b->VBO);
	glBufferData(GL_ARRAY_BUFFER, NumVert * Dimention * sizeof(GLfloat), Data, GL_DYNAMIC_DRAW);

	if (Normal) {
		glGenBuffers(1, &b->NBO);
		glBindBuffer(GL_ARRAY_BUFFER, b->NBO);
		glBufferData(GL_ARRAY_BUFFER, NumVert * Dimention * sizeof(GLfloat), Normal, GL_DYNAMIC_DRAW);
	}

	return b;
}
tnsBatch* tnsCreateBatchi(u32bit NumVert, int Dimention, int* Data) {
	tnsBatch* b = CreateNew(tnsBatch);
	b->Dimention = Dimention;
	b->NumVert = NumVert;

	glGenBuffers(1, &b->VBO);
	glBindBuffer(GL_ARRAY_BUFFER, b->VBO);
	glBufferData(GL_ARRAY_BUFFER, NumVert * Dimention * sizeof(GLint), Data, GL_DYNAMIC_DRAW);

	return b;
}

tnsBranchedCommand* tnsCreateCommand(tnsBatch* b, u32bit ElementCount,int Dimention, GLenum DrawAs,int InternalMode, u32bit* Elements) {
	tnsBranchedCommand* bc = CreateNew(tnsBranchedCommand);
	bc->ElementCount = ElementCount;
	bc->DrawAs = DrawAs;
	bc->InternalMode = InternalMode;
	bc->Dimention = Dimention;

	glGenBuffers(1, &bc->EBO);
	//glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, bc->EBO);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, ElementCount * Dimention * sizeof(GLuint), Elements, GL_DYNAMIC_DRAW);
	//int a = GetLastError();
	lstAppendItem(&b->Branches, bc);

	return bc;
}
void tnsDeleteBatch(tnsBatch* b) {
	tnsBranchedCommand *bc, *NextBc;

	glDeleteBuffers(1, &b->VBO);
	if(b->NBO > -1) glDeleteBuffers(1, &b->NBO);

	for (bc = b->Branches.pFirst; bc; bc = NextBc) {
		NextBc = bc->Item.pNext;
		lstRemoveItem(&b->Branches, bc);
		if(bc->EBO > -1) glDeleteBuffers(1, &bc->EBO);
		FreeMem(bc);
	}

	FreeMem(b);
}


//==========================*=======================[Texture]

int tnsGetTextureMemoryComponetCount(tnsTexture* t) {
	int Comp = 0;
	int CompSize;
	switch (t->GLTexBitsType){
	case GL_RGB:
		Comp = 3;
		break;
	case GL_RGBA:
		Comp = 4;
		break;
	}
	return t->Width*t->Height*Comp;
}

int tnsInit2DRectTextureAdvanced(tnsTexture* t, GLuint InternalFormat,GLuint DataFormat,GLuint DataType, int w, int h, void* bits) {
	if (!t) return 0;

	glGenTextures(1, &t->GLTexHandle);

	glBindTexture(GL_TEXTURE_RECTANGLE, t->GLTexHandle);

	t->Width = w;
	t->Height = h;
	t->GLTexBitsType = InternalFormat;
	t->DataFormat = DataFormat;
	t->DataType = DataType;
	t->GLTexType = GL_TEXTURE_RECTANGLE;

	glTexImage2D(GL_TEXTURE_RECTANGLE,
		0,
		InternalFormat,
		w, h,
		0,
		DataFormat,
		DataType,
		bits);

	if (InternalFormat != GL_DEPTH_COMPONENT) {
		glTexParameteri(GL_TEXTURE_RECTANGLE, GL_TEXTURE_WRAP_S, GL_CLAMP);
		glTexParameteri(GL_TEXTURE_RECTANGLE, GL_TEXTURE_WRAP_T, GL_CLAMP);
		glTexParameteri(GL_TEXTURE_RECTANGLE, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_RECTANGLE, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		//glTexEnvi(GL_TEXTURE_2D, GL_TEXTURE_ENV_MODE, GL_REPLACE);
	}

	glBindTexture(GL_TEXTURE_RECTANGLE, 0);

	return 1;
}
int tnsInit2DTexture(tnsTexture* t, GLuint BitsFormat, int w, int h, void* bits){
	if (!t) return 0;

	glGenTextures(1, &t->GLTexHandle);

	glBindTexture(GL_TEXTURE_2D, t->GLTexHandle);

	t->Width = w;
	t->Height = h;
	t->GLTexBitsType = t->DataFormat = BitsFormat;
	t->DataType = GL_UNSIGNED_BYTE;
	t->GLTexType = GL_TEXTURE_2D;

	glTexImage2D(GL_TEXTURE_2D,
		0,
		BitsFormat,
		w, h,
		0,
		BitsFormat,
		GL_UNSIGNED_BYTE,
	    bits);

	if (BitsFormat != GL_DEPTH_COMPONENT) {
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		//glTexEnvi(GL_TEXTURE_2D, GL_TEXTURE_ENV_MODE, GL_REPLACE);
	}

	glBindTexture(GL_TEXTURE_2D, 0);

	return 1;
}
int tnsInit2DTextureMultisample(tnsTexture* t, GLuint BitsFormat, int w, int h, void* bits) {
	if (!t) return 0;

	glGenTextures(1, &t->GLTexHandle);

	glBindTexture(GL_TEXTURE_2D_MULTISAMPLE, t->GLTexHandle);

	t->Width = w;
	t->Height = h;
	t->GLTexBitsType = t->DataFormat = BitsFormat;
	t->DataType = GL_UNSIGNED_BYTE;
	t->GLTexType = GL_TEXTURE_2D_MULTISAMPLE;

	glTexImage2DMultisample(GL_TEXTURE_2D_MULTISAMPLE, 4, BitsFormat, w, h, GL_TRUE);

	glBindTexture(GL_TEXTURE_2D_MULTISAMPLE, 0);

	return 1;
}
int tnsInit2DTextureSupersample(tnsTexture* t, GLuint BitsFormat, int w, int h, void* bits) {
	if (!t) return 0;

	glGenTextures(1, &t->GLTexHandle);

	glBindTexture(GL_TEXTURE_2D_MULTISAMPLE, t->GLTexHandle);

	t->Width = w;
	t->Height = h;
	t->GLTexBitsType = t->DataFormat = BitsFormat;
	t->DataType = GL_UNSIGNED_BYTE;
	t->GLTexType = GL_TEXTURE_2D_MULTISAMPLE;

	glTexImage2DMultisample(GL_TEXTURE_2D_MULTISAMPLE, 16, BitsFormat, w, h, GL_TRUE);

	glBindTexture(GL_TEXTURE_2D_MULTISAMPLE, 0);

	return 1;
}
int tnsInit2DDepthBufferMultisample(tnsTexture* t, int w, int h) {
	if (!t) return 0;

	glGenRenderbuffers(1, &t->GLTexHandle);

	glBindRenderbufferEXT(GL_RENDERBUFFER, t->GLTexHandle);

	t->Width = w;
	t->Height = h;
	t->GLTexBitsType = t->DataFormat = GL_DEPTH_COMPONENT32F;
	t->DataType = GL_FLOAT;
	t->GLTexType = GL_RENDERBUFFER;

	glRenderbufferStorageMultisample(GL_RENDERBUFFER, 4, GL_DEPTH_COMPONENT32F, w, h);

	glBindRenderbufferEXT(GL_RENDERBUFFER, 0);

	return 1;
}
int tnsInit2DDepthBufferSupersample(tnsTexture* t, int w, int h) {
	if (!t) return 0;

	glGenRenderbuffers(1, &t->GLTexHandle);

	glBindRenderbufferEXT(GL_RENDERBUFFER, t->GLTexHandle);

	t->Width = w;
	t->Height = h;
	t->GLTexBitsType = t->DataFormat = GL_DEPTH_COMPONENT32F;
	t->DataType = GL_FLOAT;
	t->GLTexType = GL_RENDERBUFFER;

	glRenderbufferStorageMultisample(GL_RENDERBUFFER, 16, GL_DEPTH_COMPONENT32F, w, h);

	glBindRenderbufferEXT(GL_RENDERBUFFER, 0);

	return 1;
}



void tnsFeed2DRectTextureData(tnsTexture* t, void* data) {
	glBindTexture(GL_TEXTURE_RECTANGLE, t->GLTexHandle);

	glTexImage2D(GL_TEXTURE_RECTANGLE,
		0,
		t->GLTexBitsType,
		t->Width, t->Height,
		0,
		t->DataFormat?t->DataFormat:t->GLTexBitsType,
		t->DataType?t->DataType:GL_UNSIGNED_BYTE,
		data);

	glBindTexture(GL_TEXTURE_RECTANGLE, 0);
}
void tnsFeed2DTextureData(tnsTexture* t, void* data) {
	glBindTexture(GL_TEXTURE_2D, t->GLTexHandle);

	glTexImage2D(GL_TEXTURE_2D,
		0,
		t->GLTexBitsType,
		t->Width, t->Height,
		0,
		t->DataFormat ? t->DataFormat : t->GLTexBitsType,
		t->DataType ? t->DataType : GL_UNSIGNED_BYTE,
		data);

	glBindTexture(GL_TEXTURE_2D, 0);
}
void tnsConfigure2DTexture(tnsTexture* t, GLuint BitsFormat, int w, int h, void* bits){
	glBindTexture(GL_TEXTURE_2D, t->GLTexHandle);

	t->Width = w;
	t->Height = h;
	t->GLTexBitsType = BitsFormat;
	t->GLTexType = GL_TEXTURE_2D;

	glTexImage2D(GL_TEXTURE_2D,
		0,
		BitsFormat,
		w, h,
		0,
		BitsFormat,
		GL_UNSIGNED_BYTE,
		bits);
	glBindTexture(GL_TEXTURE_2D, 0);
}
tnsTexture* tnsCreate2DRectTextureAdvanced(GLuint InternalFormat, GLuint DataFormat, GLuint DataType, int w, int h, void* bits) {
	tnsTexture* tex = CreateNew(tnsTexture);
	tnsInit2DRectTextureAdvanced(tex, InternalFormat, DataFormat, DataType, w, h, bits);
	nulNotifyUsers("tns.texture_list");
	lstAppendItem(tKnlGetTextureList(), tex);
	return tex;
};
tnsTexture* tnsCreate2DTexture(GLuint BitsFormat, int w, int h, void* bits){
	tnsTexture* tex = CreateNew(tnsTexture);
	tnsInit2DTexture(tex, BitsFormat, w, h, bits);
	nulNotifyUsers("tns.texture_list");
	lstAppendItem(tKnlGetTextureList(), tex);
	return tex;
};
tnsTexture* tnsCreate2DTextureMultisample(GLuint BitsFormat, int w, int h, void* bits) {
	tnsTexture* tex = CreateNew(tnsTexture);
	if (BitsFormat == GL_DEPTH_COMPONENT)
		tnsInit2DDepthBufferMultisample(tex, w, h);
	else tnsInit2DTextureMultisample(tex, BitsFormat, w, h, bits);
	lstAppendItem(tKnlGetTextureList(), tex);
	return tex;
};
tnsTexture* tnsCreate2DDepthBufferMultisample(GLuint BitsFormat, int w, int h, void* bits) {
	tnsTexture* tex = CreateNew(tnsTexture);
	tnsInit2DDepthBufferMultisample(tex, w, h);
	lstAppendItem(tKnlGetTextureList(), tex);
	return tex;
};
tnsTexture* tnsCreate2DTextureSupersample(GLuint BitsFormat, int w, int h, void* bits) {
	tnsTexture* tex = CreateNew(tnsTexture);
	if (BitsFormat == GL_DEPTH_COMPONENT)
		tnsInit2DDepthBufferMultisample(tex, w, h);
	else tnsInit2DTextureSupersample(tex, BitsFormat, w, h, bits);
	lstAppendItem(tKnlGetTextureList(), tex);
	return tex;
};
tnsTexture* tnsCreate2DDepthBufferSupersample(GLuint BitsFormat, int w, int h, void* bits) {
	tnsTexture* tex = CreateNew(tnsTexture);
	tnsInit2DDepthBufferSupersample(tex, w, h);
	lstAppendItem(tKnlGetTextureList(), tex);
	return tex;
};
void tnsCopyScreenTo2DTexture(tnsTexture* target, int x_lower_left, int y_lower_left, int w, int h){
	glBindTexture(GL_TEXTURE_2D, target->GLTexHandle);
	glReadBuffer(GL_BACK);
	glCopyTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, x_lower_left, y_lower_left, w, h);
	glBindTexture(GL_TEXTURE_2D, 0);
}

void tnsActiveTexture(GLenum tex) {
	if(T->GlTextureSets!=tex)glActiveTexture(tex);
	T->GlTextureSets = tex;
}
int a;
void tnsBindTexture0(tnsTexture* t){
	if (!t) return;
	glBindTexture(t->GLTexType, t->GLTexHandle); a = glGetError();
	if (T->CurrentShader->texture0Index)glUniform1i(T->CurrentShader->texture0Index, 0);
	T->Texture0Enabled = t;
}
void tnsUseTexture0(tnsTexture* t){
	T->StateTexture0 = t;
	T->UsingCommand->Texture0 = t;
}
void tnsBindTexture1(tnsTexture* t) {
	if (!t) return;
	glBindTexture(t->GLTexType, t->GLTexHandle); a = glGetError();
	if(T->CurrentShader->texture1Index)glUniform1i(T->CurrentShader->texture1Index, 1);
	T->Texture1Enabled = t;
}
void tnsUseTexture1(tnsTexture* t) {
	T->StateTexture1 = t;
	T->UsingCommand->Texture1 = t;
}
void tnsBindTexture2(tnsTexture* t) {
	if (!t) return;
	glBindTexture(t->GLTexType, t->GLTexHandle); a = glGetError();
	if (T->CurrentShader->texture2Index)glUniform1i(T->CurrentShader->texture2Index, 2);
	T->Texture2Enabled = t;
}
void tnsUseTexture2(tnsTexture* t) {
	T->StateTexture2 = t;
	T->UsingCommand->Texture2 = t;
}
void tnsBindTexture3(tnsTexture* t) {
	if (!t) return;
	glBindTexture(t->GLTexType, t->GLTexHandle); a = glGetError();
	if (T->CurrentShader->texture3Index)glUniform1i(T->CurrentShader->texture3Index, 3);
	T->Texture3Enabled = t;
}
void tnsBindTexture4(tnsTexture* t) {
	if (!t) return;
	glBindTexture(t->GLTexType, t->GLTexHandle); a = glGetError();
	if (T->CurrentShader->texture4Index)glUniform1i(T->CurrentShader->texture4Index, 4);
	T->Texture4Enabled = t;
}

void tnsUnbindTexture0(tnsTexture* t){
	if (!T->Texture0Enabled) return;
	glBindTexture(T->Texture0Enabled->GLTexType, 0);
	T->Texture0Enabled = 0;
}
void tnsUnbindTexture1(tnsTexture* t) {
	if (!T->Texture1Enabled) return;
	glBindTexture(T->Texture1Enabled->GLTexType, 0);
	T->Texture1Enabled = 0;
}
void tnsUnbindTexture2(tnsTexture* t) {
	if (!T->Texture2Enabled) return;
	glBindTexture(T->Texture2Enabled->GLTexType, 0);
	T->Texture2Enabled = 0;
}
void tnsUnbindTexture3(tnsTexture* t) {
	if (!T->Texture3Enabled) return;
	glBindTexture(T->Texture3Enabled->GLTexType, 0);
	T->Texture3Enabled = 0;
}
void tnsUnbindTexture4(tnsTexture* t) {
	if (!T->Texture4Enabled) return;
	glBindTexture(T->Texture4Enabled->GLTexType, 0);
	T->Texture4Enabled = 0;
}

void tnsDraw2DTextureDirectly(tnsTexture* t, int x_upper_right, int y_upper_right, int w, int h){
	real Verts[8];
	real UV[8] = {
		0.0f, 1.0f,
		1.0f, 1.0f,
		1.0f, 0.0f,
		0.0f, 0.0f
	};
	tnsMakeQuad2d(Verts,
		x_upper_right, y_upper_right,
		x_upper_right + w, y_upper_right,
		x_upper_right + w, y_upper_right + h,
		x_upper_right, y_upper_right + h);

	tnsUseTexture0(t);
	tnsVertexArray2d(Verts, 4);
	tnsTexCoordArray2d(UV, 4);
	tnsPackAs(GL_TRIANGLE_FAN);
}
void tnsDraw2DTextureArg(tnsTexture* t,
	real x_upper_right, real y_upper_right, int w, int h,
	real* MultiplyColor,
	real LPadding, real RPadding, real TPadding, real BPadding) {
	real Verts[8];
	real UV[8] = {
		0.0f+ LPadding, 1.0f- TPadding,
		1.0f- RPadding, 1.0f- TPadding,
		1.0f- RPadding, 0.0f+ BPadding,
		0.0f+ LPadding, 0.0f+ BPadding
	};
	tnsMakeQuad2d(Verts,
		x_upper_right, y_upper_right,
		x_upper_right + w, y_upper_right,
		x_upper_right + w, y_upper_right + h,
		x_upper_right, y_upper_right + h);

	if (MultiplyColor) tnsColor4dv(MultiplyColor);

	tnsUseTexture0(t);
	tnsVertexArray2d(Verts, 4);
	tnsTexCoordArray2d(UV, 4);
	tnsPackAs(GL_TRIANGLE_FAN);
}
void tnsDeleteTexture(tnsTexture* t){
	nListHandle* lst = tKnlGetTextureList();

	if (!t)
		return;

	if (t->GLTexType == GL_RENDERBUFFER) {
		glBindRenderbufferEXT(GL_RENDERBUFFER, t->GLTexHandle);
		glDeleteRenderbuffersEXT(1, &t->GLTexHandle);
	}else {
		glBindTexture(t->GLTexType, 0);
		glDeleteTextures(1, &t->GLTexHandle);
	}
	lstRemoveItem(lst, t);
	if (t->DrawData) FreeMem(t->DrawData);
	FreeMem(t);
}


int tnsTextureMemorySize(tnsTexture* t, int Mem) {
	int ElemSize;

	if (Mem) return t->Width*t->Height * sizeof(void*);

	switch (t->GLTexBitsType) {
	default:
	case GL_RED: ElemSize = sizeof(char) * 1; break;
	case GL_RG: ElemSize = sizeof(char) * 2; break;
	case GL_RGB: ElemSize = sizeof(char) * 3; break;
	case GL_RGBA: ElemSize = sizeof(char) * 4; break;
	case GL_RGBA32F: ElemSize = sizeof(float) * 4; break;
	case GL_DEPTH_COMPONENT32F: ElemSize = sizeof(float); break;
	}

	t->ElemSize = ElemSize;

	return t->Width*t->Height*ElemSize;
}

void tnsCreateTextureReadbackBuffer(tnsTexture* t){
	if (t->SamplePtr) {
		tnsTextureSample* ts;
		memset(t->SamplePtr, 0, tnsTextureMemorySize(t, 1));
		while (ts = lstPopItem(&t->ErasedSamples)) FreeMem(ts);
		while (ts = lstPopItem(&t->PendingSamples)) FreeMem(ts);
	}else {
		t->TextureReadBack = calloc(1, tnsTextureMemorySize(t, 0));
		t->SamplePtr = calloc(1, tnsTextureMemorySize(t, 1));
	}
	t->RBWidth = t->Width;
	t->RBHeight = t->Height;
}
void tnsDeleteTextureReadbackBuffer(tnsTexture* t) {
	FreeMem(t->TextureReadBack);
	FreeMem(t->SamplePtr);
	tnsTextureSample* ts;
	while (ts = lstPopItem(&t->ErasedSamples)) FreeMem(ts);
	while (ts = lstPopItem(&t->PendingSamples)) FreeMem(ts);
	t->RBWidth = t->RBHeight = 0;
}
void tnsReadbackTexture(tnsTexture* t) {
	if (t->RBWidth != t->Width || t->RBHeight != t->Height)
		tnsDeleteTextureReadbackBuffer(t);
	
	tnsCreateTextureReadbackBuffer(t);

	glActiveTexture(GL_TEXTURE0);
	glBindTexture(t->GLTexType, t->GLTexHandle);
	glGetTexImage(t->GLTexType, 0, t->DataFormat, t->DataType, t->TextureReadBack);
	int a = glGetError();
	glBindTexture(t->GLTexType, 0);

	int w, h;
	tnsTextureSample* ts;
	u8bit sample;
	for (h = 0; h < t->Height; h++) {
		for (w = 0; w < t->Width; w++) {
			int index = h*t->Width + w;
			if ((sample = *tnsTextureSampleU8(t, w, h)) > TNS_SNAKE_STRENGTH_LIMIT) {
				ts = CreateNew(tnsTextureSample);
				lstAppendItem(&t->PendingSamples, ts);
				t->SamplePtr[h*t->Width + w] = ts;
				ts->X = w;
				ts->Y = h;
				ts->Sample = sample;
			}
		}
	}
	
}


//====================================================[NEW RENDER KERNEL]
//=================[Immediate-style api]

void tnsColor4d(real r, real g, real b, real a){
	tnsCommand* c = T->UsingCommand;
	T->StateColor[0] = r;
	T->StateColor[1] = g;
	T->StateColor[2] = b;
	T->StateColor[3] = a;
	c->UniformColor[0] = r;
	c->UniformColor[1] = g;
	c->UniformColor[2] = b;
	c->UniformColor[3] = a;
}
void tnsColor4dv(real* rgba){
	tnsCommand* c = T->UsingCommand;
	T->StateColor[0] = rgba[0];
	T->StateColor[1] = rgba[1];
	T->StateColor[2] = rgba[2];
	T->StateColor[3] = rgba[3];
	c->UniformColor[0] = rgba[0];
	c->UniformColor[1] = rgba[1];
	c->UniformColor[2] = rgba[2];
	c->UniformColor[3] = rgba[3];
}

void tnsPolygonMode(GLenum PolyMode){
	glPolygonMode(GL_FRONT_AND_BACK, PolyMode);
	glPolygonMode(GL_FRONT_AND_BACK, PolyMode);
	T->StatePolyMode = PolyMode;
}
void tnsShadeMode(GLenum ShadeMode){
	glShadeModel(ShadeMode);
	T->StateShadeMode = ShadeMode;
}
void tnsVertex3d(real x, real y, real z){
	tnsCommand* c = T->UsingCommand;
	short vend = c->VertEnd;
	GLfloat* varr = T->Vert;

	c->UseVert = 1;

	if (!c->Dimentions) c->Dimentions = 3;

	if (c->Dimentions == 3){
		varr[vend]  = x;
		varr[vend+1] = y;
		varr[vend+2] = z;
		c->NumVert++;
		T->NextVert += 3;
		c->VertEnd += 3;
	}else{
		varr[vend] = x;
		varr[vend + 1] = y;
		c->NumVert++;
		T->NextVert += 2;
		c->VertEnd += 2;
	}
}
void tnsVertex2d(real x, real y){
	tnsCommand* c = T->UsingCommand;
	short vend = c->VertEnd;
	GLfloat* varr = T->Vert;

	c->UseVert = 1;

	if (!c->Dimentions) c->Dimentions = 2;

	if (c->Dimentions == 2){
		varr[vend] = x;
		varr[vend + 1] = y;
		c->NumVert++;
		T->NextVert += 2;
		c->VertEnd+=2;
	}
	else{
		tnsVertex3d(x, y, 0.0f);
	}
}
void tnsVertexArray2d(real* verts, int amount){
	tnsCommand* c = T->UsingCommand;
	int trans = 2 * amount;
	short vend = c->VertEnd;
	GLfloat* varr = T->Vert;

	c->UseVert = 1;

	if (!c->Dimentions) c->Dimentions = 2;

	if (c->Dimentions == 2){
		int i;
		for (i = 0; i < trans; i++) {
			varr[vend] = verts[i];
			vend++;
		}
		//memcpy(&varr[vend], verts, trans*sizeof(real));
		c->VertEnd  += trans;
		c->NumVert += amount;
		T->NextVert += trans;
	}
}
void tnsVertexArray3d(real* verts, int amount){
	tnsCommand* c = T->UsingCommand;
	int trans = 3 * amount;
	short vend = c->VertEnd;
	GLfloat* varr = T->Vert;

	c->UseVert = 1;

	if (!c->Dimentions) c->Dimentions = 3;

	if (c->Dimentions == 3){
		int i;
		for (i = 0; i < trans; i++) {
			varr[vend] = verts[i];
			vend++;
		}
		//memcpy(&varr[vend], verts, trans*sizeof(real));
		c->VertEnd  += trans;
		c->NumVert += amount;
		T->NextVert += trans;
	}
}
void tnsColorArray4d(real* colors, int amount){
	tnsCommand* c = T->UsingCommand;
	int trans = 4 * amount;
	GLfloat* carr = T->Color;
	int ofst = c->ColorEnd;

	c->UseColor = 1;

	int i;
	for (i = 0; i < trans; i++) {
		carr[ofst] = colors[i];
		ofst++;
	}
	//memcpy(&T->Color[c->ColorEnd], colors, trans*sizeof(real));
	c->ColorEnd += trans;
	T->NextColor += trans;
}
void tnsNormalArray3d(real* normals, int amount){
	tnsCommand* c = T->UsingCommand;
	int trans = 3 * amount;
	GLfloat* narr = T->Normal;
	int ofst = c->NormalEnd;

	c->UseNormal = 1;

	int i;
	for (i = 0; i < trans; i++) {
		narr[ofst] = normals[i];
		ofst++;
	}
	//memcpy(&T->Normal[c->NormalEnd], normals, trans*sizeof(real));
	c->NormalEnd += trans;
	T->NextNormal += trans;
}
void tnsTexCoordArray2d(real* coords, int amount){
	tnsCommand* c = T->UsingCommand;
	int trans = 2 * amount;
	GLfloat* carr = T->TexCoord;
	int ofst = c->TexCoordEnd;

	c->UseTexCoord = 1;

	int i;
	for (i = 0; i < trans; i++) {
		carr[ofst] = coords[i];
		ofst++;
	}
	//memcpy(&T->TexCoord[c->TexCoordEnd], coords, trans*sizeof(real));
	c->TexCoordEnd += trans;
	T->NextTexCoord += trans;
}
void tnsIndexArray(GLushort* index,short amount){
	tnsCommand* c = T->UsingCommand;

	if (c->UseIndex) return;

	c->UseIndex = 1;

	memcpy(&T->Index[c->IndexEnd], index, amount*sizeof(GLushort));
	c->IndexEnd += amount;
	c->NumIndex += amount;
	T->NextIndex += amount;
}

void tnsPackAs(GLenum Mode){
	tnsCommand* c = T->UsingCommand;
	tnsCommand* nc;

	c->Mode = Mode;
	T->UsingCommand++;
	if (!T->UsingCommand)
		return;
	nc = T->UsingCommand;

	memset(nc, 0, sizeof(tnsCommand));

	nc->VertBegin = nc->VertEnd = c->VertEnd;
	nc->NormalBegin = nc->NormalEnd = c->NormalEnd;
	nc->ColorBegin = nc->ColorEnd = c->ColorEnd;
	nc->TexCoordBegin = nc->TexCoordEnd = c->TexCoordEnd;
	nc->IndexBegin = nc->IndexEnd = c->IndexEnd;

	nc->PolyMode = c->PolyMode;
	nc->Shade = c->Shade;
	nc->Texture0 = c->Texture0;

	nc->UniformColor[0] = c->UniformColor[0];
	nc->UniformColor[1] = c->UniformColor[1];
	nc->UniformColor[2] = c->UniformColor[2];
	nc->UniformColor[3] = c->UniformColor[3];

	//memcpy(nc->UniformColor, c->UniformColor, sizeof(real) * 4);

	if (T->UsingCommand == &T->DrawingCommand[127]){
		tnsFlush();
	}
}
void tnsFlush(){
	tnsShader* cs = T->CurrentShader;
	tnsCommand* tc = &T->DrawingCommand[0];
	tnsCommand* c = tc;
	int PrevDimentions = 0;
	int LastVertBegin = 0;

	if (!c || !cs) 
		return;

	if (T->NextVert){
		glBindBuffer(GL_ARRAY_BUFFER, T->VertBufObject);
		//glEnableVertexAttribArray(cs->vertexIndex);
		glBufferSubData(GL_ARRAY_BUFFER, 0, T->NextVert*sizeof(GLfloat), T->Vert);
	}
	if (T->NextColor){
		glBindBuffer(GL_ARRAY_BUFFER, T->ColorBufObject);
		//glEnableVertexAttribArray(cs->colorIndex);
		glBufferSubData(GL_ARRAY_BUFFER, 0, T->NextColor*sizeof(GLfloat), T->Color);
	}
	if (T->NextNormal){
		glBindBuffer(GL_ARRAY_BUFFER, T->NormalBufObject);
		//glEnableVertexAttribArray(cs->colorIndex);
		glBufferSubData(GL_ARRAY_BUFFER, 0, T->NextNormal*sizeof(GLfloat), T->Normal);
	}
	if (T->NextTexCoord){
		glBindBuffer(GL_ARRAY_BUFFER, T->TexCoordBufObject);
		//glEnableVertexAttribArray(cs->colorIndex);
		glBufferSubData(GL_ARRAY_BUFFER, 0, T->NextTexCoord*sizeof(GLfloat), T->TexCoord);
	}


	for (c; c != T->UsingCommand; c++){

		if (c->PolyMode != T->StatePolyMode) tnsPolygonMode(c->PolyMode);
		if (c->Shade != T->StateShadeMode) tnsShadeMode(c->Shade);
		if (c->ReplaceShader && c->ReplaceShader != T->CurrentShader){
			tnsEnableShaderv(c->ReplaceShader);
			cs = c->ReplaceShader;
			if (!cs) return 0;
		}
		
		glBindBuffer(GL_ARRAY_BUFFER, T->VertBufObject);
		if (c->UseVert){
			glEnableVertexAttribArray(cs->vertexIndex);
			glVertexAttribPointer(cs->vertexIndex, (c->Dimentions ? c->Dimentions : PrevDimentions),
				GL_FLOAT, 0, 0, c->VertBegin*sizeof(GLfloat));
			LastVertBegin = c->VertBegin;
		}else{
			glEnableVertexAttribArray(cs->vertexIndex);
			glVertexAttribPointer(cs->vertexIndex, (c->Dimentions ? c->Dimentions : PrevDimentions),
				GL_FLOAT, 0, 0, LastVertBegin*sizeof(GLfloat));
		}

		PrevDimentions = (c->Dimentions?c->Dimentions:PrevDimentions);
		
		if (c->UseIndex){
			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, T->IndexBufObject);
			glBufferSubData(GL_ELEMENT_ARRAY_BUFFER, 0, c->NumIndex*sizeof(GLushort), &T->Index[c->IndexBegin]);
		}

		if (cs->normalIndex != -1){
			glBindBuffer(GL_ARRAY_BUFFER, T->NormalBufObject);
			if (c->UseNormal){
				glEnableVertexAttribArray(cs->normalIndex);
			}
			glVertexAttribPointer(cs->normalIndex, 3, GL_FLOAT, 0, 0, c->NormalBegin*sizeof(GLfloat));
			if (!c->UseNormal && c->Dimentions == 2){
				glDisableVertexAttribArray(cs->normalIndex);
				glVertexAttrib3f(cs->normalIndex, 0, 0, -1);
			}
		}

		if (cs->colorIndex != -1){
			glBindBuffer(GL_ARRAY_BUFFER, T->ColorBufObject);
			if (c->UseColor){
				glEnableVertexAttribArray(cs->colorIndex);
				glVertexAttribPointer(cs->colorIndex, 4, GL_FLOAT, 0, 0, c->ColorBegin*sizeof(GLfloat));
			}
			if(!c->UseColor){
				glDisableVertexAttribArray(cs->colorIndex);
				glVertexAttrib4fv(cs->colorIndex, c->UniformColor);
			}
		}
		
		if (c->Texture0 && cs->texture0Index != -1){
			tnsActiveTexture(GL_TEXTURE0);
			if (c->Texture0 != T->StateTexture0 || c->Texture0 != T->Texture0Enabled) 
				tnsBindTexture0(c->Texture0);
		}else{
			tnsUnbindTexture0(T->StateTexture0);
		}
		
		if (c->Texture1 && cs->texture1Index != -1) {
			tnsActiveTexture(GL_TEXTURE1);
			if (c->Texture1 != T->StateTexture1 || c->Texture1 != T->Texture1Enabled)
				tnsBindTexture1(c->Texture1);
		}
		else {
			tnsUnbindTexture1(T->StateTexture1);
		}

		if (c->Texture2 && cs->texture2Index != -1) {
			tnsActiveTexture(GL_TEXTURE2);
			if (c->Texture2 != T->StateTexture2 || c->Texture2 != T->Texture2Enabled)
				tnsBindTexture2(c->Texture2);
		}
		else {
			tnsUnbindTexture2(T->StateTexture2);
		}

		if (c->UseTexCoord) {
			glBindBuffer(GL_ARRAY_BUFFER, T->TexCoordBufObject);
			if (c->UseTexCoord) {
				glEnableVertexAttribArray(cs->uvIndex);
			}
			glVertexAttribPointer(cs->uvIndex, 2, GL_FLOAT, 0, 0, c->TexCoordBegin*sizeof(GLfloat));
			if (!c->UseTexCoord) {
				glDisableVertexAttribArray(cs->uvIndex);
			}
		}
		
		if (c->UseIndex){
			glDrawElements(c->Mode, c->NumIndex, GL_UNSIGNED_SHORT,0);
		}else{
		    if (c->Mode == GL_QUADS) c->Mode = GL_TRIANGLE_STRIP;
			glDrawArrays(c->Mode, 0, c->NumVert);
		}

		tnsActiveTexture(GL_TEXTURE0);
	}

	T->UsingCommand = T->DrawingCommand;
	c = T->UsingCommand;
	if (!T->UsingCommand)
		return 0;

	memset(c, 0, sizeof(tnsCommand));

	c->ColorBegin = c->ColorEnd = T->NextColor = 0;
	c->NormalBegin = c->NormalEnd = T->NextNormal = 0;
	c->TexCoordBegin = c->TexCoordEnd = T->NextTexCoord = 0;
	c->VertBegin = c->VertEnd = T->NextVert = 0;
	c->IndexBegin = c->IndexEnd = T->NextIndex = 0;

	//must --why?
	//T->BindedShader = 0;
};




//============================================================================================[offscr]


const GLuint TNS_ATTACHMENT_ARRAY_NONE[] = { GL_NONE };
const GLuint TNS_ATTACHMENT_ARRAY[] = { GL_COLOR_ATTACHMENT0 };
const GLuint TNS_ATTACHMENT_ARRAY_1[] = { GL_COLOR_ATTACHMENT1 };
const GLuint TNS_ATTACHMENT_ARRAY_2[] = { GL_COLOR_ATTACHMENT2 };
const GLuint TNS_ATTACHMENT_ARRAY_1_2[] = { GL_COLOR_ATTACHMENT1, GL_COLOR_ATTACHMENT2 };
const GLuint TNS_ATTACHMENT_ARRAY_0_1_2[] = { GL_COLOR_ATTACHMENT0, GL_COLOR_ATTACHMENT1, GL_COLOR_ATTACHMENT2 };
const GLenum TNS_WINDOW_DRAWBUFFER_ARRAY[] = { GL_BACK };


tnsOffscreen* tnsCreateOffscreenHandle(){
	tnsOffscreen* toff = CreateNew(tnsOffscreen);
	glGenFramebuffers(1, &toff->FboHandle);
	return toff;
}
void tnsAttach2DOffscreenBuffer(tnsOffscreen* target, GLuint attatchment, tnsTexture* use){
	if (!target || !use || target->FboHandle == -1 || use->GLTexHandle == -1) return;

	if (attatchment >= GL_COLOR_ATTACHMENT0 && attatchment <= GL_COLOR_ATTACHMENT15){
		if (target->pColorTextures[attatchment - GL_COLOR_ATTACHMENT0]) return;

		glBindFramebuffer(GL_FRAMEBUFFER, target->FboHandle);
		glBindTexture(use->GLTexType, use->GLTexHandle);
		glFramebufferTexture2D(GL_FRAMEBUFFER, attatchment, use->GLTexType, use->GLTexHandle, 0);

		target->pColorTextures[attatchment - GL_COLOR_ATTACHMENT0] = use;

		glBindTexture(use->GLTexType, 0);
		glBindFramebuffer(GL_FRAMEBUFFER, 0);

	}elif(attatchment == GL_DEPTH_ATTACHMENT){
		if (target->pDepthTexture) return;

		glBindFramebuffer(GL_FRAMEBUFFER, target->FboHandle);
		glBindTexture(use->GLTexType, use->GLTexHandle);
		glFramebufferTexture2D(GL_FRAMEBUFFER, attatchment, use->GLTexType, use->GLTexHandle, 0);


		//glBindRenderbufferEXT(GL_RENDERBUFFER, use->GLTexHandle);
		//glBindFramebuffer(GL_FRAMEBUFFER, target->FboHandle);
		//glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, use->GLTexHandle);

		target->pDepthTexture = use;

		glBindTexture(use->GLTexType, 0);
		glBindFramebuffer(GL_FRAMEBUFFER, 0);


		//glBindRenderbufferEXT(GL_RENDERBUFFER, 0);
		//glBindFramebuffer(GL_FRAMEBUFFER, 0);
	}
	GLenum result = glCheckFramebufferStatus(GL_FRAMEBUFFER);
	//if (result == GL_FRAMEBUFFER_COMPLETE) {
	//	printf("Framebuffer complete!\n");
	//}
	//else {
	//	printf("Framebuffer incomplete!\n");
	//}
	
}
void tnsDetach2DOffscreenBuffer(tnsOffscreen* target, GLuint which_attach_point){
	if (which_attach_point >= GL_COLOR_ATTACHMENT0 && which_attach_point <= GL_COLOR_ATTACHMENT15){
		if (target->pColorTextures[which_attach_point - GL_COLOR_ATTACHMENT0] == 0) return;

		glBindFramebuffer(GL_FRAMEBUFFER, target->FboHandle);
		glFramebufferTexture2D(GL_FRAMEBUFFER, which_attach_point, GL_TEXTURE_2D, 0, 0);

		tnsDeleteTexture(target->pColorTextures[which_attach_point - GL_COLOR_ATTACHMENT0]);

		target->pColorTextures[which_attach_point - GL_COLOR_ATTACHMENT0] = 0;
	}elif(which_attach_point == GL_DEPTH_ATTACHMENT){
		if (target->pDepthTexture) return;

		glBindFramebuffer(GL_FRAMEBUFFER, target->FboHandle);
		glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, 0, 0);

		tnsDeleteTexture(target->pDepthTexture);

		target->pDepthTexture = 0;
	}
	glBindFramebuffer(GL_FRAMEBUFFER, 0);
	glBindTexture(GL_TEXTURE_2D, 0);
}

tnsOffscreen* tnsCreate2DOffscreenBasic(int ColorBitsFormat, int w, int h, void* ColorData) {
	tnsOffscreen* toff = tnsCreateOffscreenHandle();
	tnsTexture* color = tnsCreate2DTexture(ColorBitsFormat, w, h, ColorData);

	tnsAttach2DOffscreenBuffer(toff, GL_COLOR_ATTACHMENT0, color);

	return toff;
}
tnsOffscreen* tnsCreate2DOffscreenWithDepthMultisample(int ColorBitsFormat, int w, int h, void* ColorData, void* DepthData){
	tnsOffscreen* toff = tnsCreateOffscreenHandle();
	tnsTexture* color = tnsCreate2DTextureMultisample(ColorBitsFormat, w, h, ColorData);
	tnsTexture* depth = tnsCreate2DTextureMultisample(GL_DEPTH_COMPONENT32F, w, h, DepthData);
	
	tnsAttach2DOffscreenBuffer(toff, GL_COLOR_ATTACHMENT0, color);
	tnsAttach2DOffscreenBuffer(toff, GL_DEPTH_ATTACHMENT, depth);

	return toff;
}
tnsOffscreen* tnsCreate2DOffscreenMultisample(int ColorBitsFormat, int w, int h, void* ColorData) {
	tnsOffscreen* toff = tnsCreateOffscreenHandle();
	tnsTexture* color = tnsCreate2DTextureMultisample(ColorBitsFormat, w, h, ColorData);

	tnsAttach2DOffscreenBuffer(toff, GL_COLOR_ATTACHMENT0, color);

	return toff;
}
void tnsCreate2DOffscreenMSAttachmentExtraColor(tnsOffscreen* off) {
	tnsOffscreen* toff = off;
	tnsTexture* color = tnsCreate2DTextureMultisample(
		off->pColorTextures[0]->GLTexBitsType, off->pColorTextures[0]->Width, off->pColorTextures[0]->Height, 0);

	tnsAttach2DOffscreenBuffer(toff, GL_COLOR_ATTACHMENT1, color);

	return toff;
}
void tnsCreate2DOffscreenMSAttachmentExtraNormal(tnsOffscreen* off) {
	tnsOffscreen* toff = off;
	tnsTexture* color = tnsCreate2DTextureMultisample(
		off->pColorTextures[0]->GLTexBitsType, off->pColorTextures[0]->Width, off->pColorTextures[0]->Height, 0);

	tnsAttach2DOffscreenBuffer(toff, GL_COLOR_ATTACHMENT2, color);

	return toff;
}
tnsOffscreen* tnsCreate2DOffscreenSupersample(int ColorBitsFormat, int w, int h, void* ColorData) {
	tnsOffscreen* toff = tnsCreateOffscreenHandle();
	tnsTexture* color = tnsCreate2DTextureSupersample(ColorBitsFormat, w, h, ColorData);

	tnsAttach2DOffscreenBuffer(toff, GL_COLOR_ATTACHMENT0, color);

	return toff;
}
tnsOffscreen* tnsCreate2DOffscreenWithDepthSupersample(int ColorBitsFormat, int w, int h, void* ColorData, void* DepthData) {
	tnsOffscreen* toff = tnsCreateOffscreenHandle();
	tnsTexture* color = tnsCreate2DTextureSupersample(ColorBitsFormat, w, h, ColorData);
	tnsTexture* depth = tnsCreate2DDepthBufferSupersample(GL_DEPTH_COMPONENT, w, h, DepthData);

	tnsAttach2DOffscreenBuffer(toff, GL_COLOR_ATTACHMENT0, color);
	tnsAttach2DOffscreenBuffer(toff, GL_DEPTH_ATTACHMENT, depth);

	return toff;
}

void tnsDelete2DOffscreen(tnsOffscreen* o){
	if (!o) return;

	glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);

	if (o->pColorTextures[1]) { tnsDetach2DOffscreenBuffer(o, GL_COLOR_ATTACHMENT1); tnsDeleteTexture(o->pColorTextures[1]); }
	if (o->pColorTextures[2]) { tnsDetach2DOffscreenBuffer(o, GL_COLOR_ATTACHMENT2); tnsDeleteTexture(o->pColorTextures[2]); }
	if (o->pColorTextures[3]) { tnsDetach2DOffscreenBuffer(o, GL_COLOR_ATTACHMENT3); tnsDeleteTexture(o->pColorTextures[3]); }

	tnsDetach2DOffscreenBuffer(o, GL_COLOR_ATTACHMENT0);
	tnsDetach2DOffscreenBuffer(o, GL_DEPTH_ATTACHMENT);

	tnsDeleteTexture(o->pColorTextures[0]);
	tnsDeleteTexture(o->pDepthTexture);

	glDeleteFramebuffers(1, &o->FboHandle);
	//tnsDeleteTexture(o->pStencilTexture);

	FreeMem(o);
}

void tnsDrawToOffscreen(tnsOffscreen* toff, int HowMany, GLuint* AttachmentArray){
	if (!toff) return;
	glBindFramebuffer(GL_DRAW_FRAMEBUFFER, toff->FboHandle);
	glDrawBuffers((HowMany ? HowMany : 1), (AttachmentArray ? AttachmentArray : TNS_ATTACHMENT_ARRAY));
	T->IsOffscreen = 1;
	T->BindedShader = 0;
}
void tnsDrawToExtraColorAttachment(tnsOffscreen* toff) {
	if (!toff) return;
	if (!toff->pColorTextures[1]) {
		tnsCreate2DOffscreenMSAttachmentExtraColor(toff);
	}
	glBindFramebuffer(GL_DRAW_FRAMEBUFFER, toff->FboHandle);
	glDrawBuffers(1 , TNS_ATTACHMENT_ARRAY_1);
	T->IsOffscreen = 1;
	T->BindedShader = 0;
}
void tnsDrawToExtraNormalAttachment(tnsOffscreen* toff) {
	if (!toff) return;
	if (!toff->pColorTextures[2]) {
		tnsCreate2DOffscreenMSAttachmentExtraNormal(toff);
	}
	glBindFramebuffer(GL_DRAW_FRAMEBUFFER, toff->FboHandle);
	glDrawBuffers(1, TNS_ATTACHMENT_ARRAY_2);
	T->IsOffscreen = 1;
	T->BindedShader = 0;
}
void tnsDrawToAllExtraAttachments(tnsOffscreen* toff) {
	if (!toff) return;
	if (!toff->pColorTextures[1]) {
		tnsCreate2DOffscreenMSAttachmentExtraColor(toff);
	}
	if (!toff->pColorTextures[2]) {
		tnsCreate2DOffscreenMSAttachmentExtraNormal(toff);
	}
	glBindFramebuffer(GL_DRAW_FRAMEBUFFER, toff->FboHandle);
	glDrawBuffers(2, TNS_ATTACHMENT_ARRAY_1_2);
	T->IsOffscreen = 1;
	T->BindedShader = 0;
}
void tnsDrawToOffscreenOnlyBind(tnsOffscreen* toff, int HowMany, GLuint* AttachmentArray) {
	if (!toff) return;
	glBindFramebuffer(GL_DRAW_FRAMEBUFFER, toff->FboHandle);
}
void tnsReadFromOffscreen(tnsOffscreen* toff) {
	if (!toff) return;
	glBindFramebuffer(GL_READ_FRAMEBUFFER, toff->FboHandle);
}
void tnsPassColorBetweenOffscreens(tnsOffscreen* from,tnsOffscreen* to,
	GLint srcX0, GLint srcY0, GLint srcX1, GLint srcY1, GLint dstX0, GLint dstY0, GLint dstX1, GLint dstY1, GLenum FilterMode) {
	if (!from || !to) return;
	glBlitFramebufferEXT(
		srcX0, srcY0, srcX1, srcY1,
		dstX0, dstY0, dstX1, dstY1,
		GL_COLOR_BUFFER_BIT, FilterMode);
}
void tnsDrawToScreen(){
	glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
	glDrawBuffer(GL_BACK);
	T->IsOffscreen = 0;
	T->BindedShader = 0;
}






//===========================================================================[FONT]


tnsFontManager* FM;

void tnsSetuptnsFontManager(){
	FM = CreateNew(tnsFontManager);
};	

void tnsAddToFontManager(HDC hdc, unsigned dliststart, char* name){
	tnsFont* f = CreateNew(tnsFont);
	if (!f)
		SEND_PANIC_ERROR("Can't Create tnsFont!");

	//f->dListBegin=dliststart;
	f->fontName = name;
	//f->hdc=hdc;

	lstAppendItem(&FM->Fonts, f);
};

int next_p2(int a){
	int rval = 1;
	while (rval<a) rval <<= 1;
	return rval;
}

int tnsLoadSystemFont(const char* name,const char* IconName, unsigned int size){
	unsigned char i;

	tnsFont* f = CreateNew(tnsFont);

	f->height = (int)((real)size*1.4);
	f->fontName = name;
	f->IconFontName = IconName;

	if (FT_Init_FreeType(&f->ftlib))
		SEND_PANIC_ERROR("Can't Load Main Font:Freetype Init Failed!");
	if (FT_New_Face(f->ftlib, name, 0, &f->ftface))
		SEND_PANIC_ERROR("Can't Load Main Font:Freetype Can't Init Face");

	//if (FT_Init_FreeType(&f->ftlib))
	//	SEND_PANIC_ERROR("Can't Load Icon Font:Freetype Init Failed!");
	if (FT_New_Face(f->ftlib, IconName, 0, &f->iconftface))
		SEND_PANIC_ERROR("Can't Load Icon Font:Freetype Can't Init Face");


	FT_Select_Charmap(f->ftface, FT_ENCODING_UNICODE);
	FT_Set_Char_Size(f->ftface, size << 6, size << 6, 96, 96);

	FT_Select_Charmap(f->iconftface, FT_ENCODING_UNICODE);
	FT_Set_Char_Size(f->iconftface, size << 6, size << 6, 96, 96);


	tnsInit2DTexture(&f->TexBuffer, GL_ALPHA, TNS_FONT_BUFFER_W, TNS_FONT_BUFFER_H, 0);

	lstAppendItem(&FM->Fonts, f);

	tnsUseFont(name);
};
int tnsLoadVectorGraphPackage(const char* name, unsigned int size) {
	unsigned char i;

	tnsFont* f = CreateNew(tnsFont);

	f->height = (int)((real)size*1.4);
	f->fontName = name;
	
	if (FT_Init_FreeType(&f->ftlib))
		SEND_PANIC_ERROR("Can't Load Main Font:Freetype Init Failed!");
	if (FT_New_Face(f->ftlib, name, 0, &f->ftface))
		SEND_PANIC_ERROR("Can't Load Main Font:Freetype Can't Init Face");

	FT_Select_Charmap(f->ftface, FT_ENCODING_UNICODE);
	FT_Set_Char_Size(f->ftface, size << 6, size << 6, 96, 96);

	//tnsInit2DTexture(&f->TexBuffer, GL_ALPHA, TNS_FONT_BUFFER_W, TNS_FONT_BUFFER_H, 0);

	//lstAppendItem(&FM->Fonts, f);

	//tnsUseFont(name);

	FM->VectorsGrapghs = f;
};


int tfntBufferWidthEnough(int total_width, int current_width, int this_width){
	return (current_width + this_width < total_width);
}

void tfntApplyCaracterBufferOffset(tnsFont* f, tnsFontSingleCharacter* fsc){
	if (!tfntBufferWidthEnough(TNS_FONT_BUFFER_W, f->CurrentX, fsc->width)){
		f->CurrentY += f->height*2;
		f->CurrentX = 0;
	}
	fsc->bufferx = f->CurrentX;
	fsc->buffery = f->CurrentY;
	f->CurrentX += fsc->width;
}

tnsFontSingleCharacter* tfntFetchVectorGraphTextureIDW(int ID) {
	GLuint revel = 0;
	tnsFont* f = FM->VectorsGrapghs;
	tnsFontSingleCharacter* fsc = 0;
	FT_Glyph glyph = 0;
	FT_BitmapGlyph bitmap_glyph;
	FT_Bitmap bm;
	FT_Face face;
	int w, h, i, j;
	GLubyte* buf = 0;
	int a;
	int Size = ((int)(_ICON_SYMBOL_SIZE[ID - ICON_SYMBOL_START] * MAIN.UiRowHeight)) << 6;
	int ret;

	if (!f)
		return 0;

	if (revel = f->icons[ID].Generated)
		return &f->icons[ID];
		
	FT_Set_Char_Size(f->ftface, Size, Size, 96, 96);

	if (ret = FT_Load_Char(f->ftface, ID, FT_LOAD_TARGET_NORMAL | FT_LOAD_NO_HINTING | FT_LOAD_NO_BITMAP))
		return 0;

	fsc = &f->icons[ID];

	face = f->ftface;


	if (FT_Get_Glyph(face->glyph, &glyph))
		return 0;

	FT_Render_Glyph(face->glyph, FT_RENDER_MODE_NORMAL);

	FT_Glyph_To_Bitmap(&glyph, FT_RENDER_MODE_NORMAL, 0, 1);
	bitmap_glyph = glyph;
	bm = bitmap_glyph->bitmap;

	w = bm.width;
	h = bm.rows;
	fsc->charID = ID;
	fsc->width = w;
	fsc->height = h;
	fsc->advx = face->glyph->advance.x / 64.0f;
	fsc->advy = face->size->metrics.y_ppem;
	fsc->deltax = bitmap_glyph->left;
	fsc->deltay = bitmap_glyph->top - h;

	//tfntApplyCaracterBufferOffset(f, fsc);

	//glBindTexture(GL_TEXTURE_2D, f->TexBuffer.GLTexHandle);

	buf = CreateNewBuffer(GLubyte, w*h);

	for (j = 0; j<h; j++) {
		for (i = 0; i<w; i++) {
			unsigned char _vl = bm.buffer[i + w*j];
			buf[i + (h - j - 1)*w] = _vl;
		}
	}

	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

	fsc->Tex = tnsCreate2DTexture(GL_ALPHA, w, h, buf);

	glBindTexture(GL_TEXTURE_2D, 0);

	//glTexImage2D(GL_TEXTURE_2D, 0, GL_ALPHA, w, h, 0, GL_ALPHA, GL_UNSIGNED_BYTE, buf);
	//glBindTexture(GL_TEXTURE_2D, 0);

	FreeMem(buf);

	fsc->Generated = 1;
	return fsc;
}
tnsFontSingleCharacter* tfntFetchCharTextureIDW(wchar_t ch,char IsIcon){
	GLuint revel = 0;
	tnsFont* f = FM->UsingFont;
	tnsFontSingleCharacter* fsc = 0;
	FT_Glyph glyph=0;
	FT_BitmapGlyph bitmap_glyph;
	FT_Bitmap bm;
	FT_Face face;
	int w, h, i, j;
	GLubyte* buf = 0;

	if (!f)
		return 0;

	if (IsIcon) {
		if (revel = f->icons[ch].Generated)
			return &f->icons[ch];

		if (FT_Load_Char(f->iconftface, ch, FT_LOAD_TARGET_NORMAL | FT_LOAD_NO_HINTING | FT_LOAD_NO_BITMAP))
			return 0;

		fsc = &f->icons[ch];

		face = f->iconftface;

	}else {
		if (revel = f->characters[ch].Generated)
			return &f->characters[ch];

		if (FT_Load_Char(f->ftface, ch, FT_LOAD_TARGET_NORMAL | FT_LOAD_NO_HINTING | FT_LOAD_NO_BITMAP))
			return 0;

		fsc = &f->characters[ch];

		face = f->ftface;
	}

	if (FT_Get_Glyph(face->glyph, &glyph))
		return 0;

	FT_Render_Glyph(face->glyph, FT_RENDER_MODE_NORMAL);

	FT_Glyph_To_Bitmap(&glyph, FT_RENDER_MODE_NORMAL, 0, 1);
	bitmap_glyph = glyph;
	bm = bitmap_glyph->bitmap;

	w = bm.width;
	h = bm.rows;
	fsc->charID = ch;
	fsc->width = w;
	fsc->height = h;
	fsc->advx = face->glyph->advance.x / 64.0f;
	fsc->advy = face->size->metrics.y_ppem;
	fsc->deltax = bitmap_glyph->left;
	fsc->deltay = bitmap_glyph->top - h;

	tfntApplyCaracterBufferOffset(f, fsc);

	glBindTexture(GL_TEXTURE_2D, f->TexBuffer.GLTexHandle);

	buf = CreateNewBuffer(GLubyte, w*h);

	for (j = 0; j<h; j++){
		for (i = 0; i<w; i++){
			unsigned char _vl = bm.buffer[i + w*j];
			buf[i + (h - j - 1)*w] = _vl;
		}
	}

	//glPixelStorei(GL_UNPACK_LSB_FIRST, GL_FALSE);
	//glPixelStorei(GL_UNPACK_ROW_LENGTH, 0);
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	glTexSubImage2D(GL_TEXTURE_2D, 0, fsc->bufferx, fsc->buffery, w, h, GL_ALPHA, GL_UNSIGNED_BYTE, buf);

	FreeMem(buf);

	fsc->Generated = 1;
	return fsc;
}
tnsFontSingleCharacter* tfntFetchCharacterW(wchar_t ch,char IsIcon){
	return tfntFetchCharTextureIDW(ch, IsIcon);
}

int CMP_NAME_IsThisFont(tnsFont* enumed, char* name){
	return (!strcmp(enumed->fontName, name));
};

int tnsUseFont(char* name){
	tnsFont* f = lstFindItem(name, CMP_NAME_IsThisFont, &FM->Fonts);

	if (!f)
		return 0;

	FM->UsingFont = f;

	return 1;
};

int tnsStringGetWidthU(wchar_t* content, nThemeState* ts,int Count){
	int sx = 0, sy = FM->UsingFont->height;
	int i;
	int C=0;
	int len = wcslen(content);
	char hb;
	int RestoreI;

	for (i = 0; i<len; i++){
		tnsFontSingleCharacter* fsc;
		if (content[i] == '\n'){
			sx = 0;
			sy += FM->UsingFont->height;
			continue;
		}
		else{
			if (content[i] != '$') {
				fsc = tfntFetchCharTextureIDW(content[i], 0);
			}
			else {
				int id=0;
				RestoreI = i;
				i++;
				while (content[i] != '$') {
					if (!(content[i] >= '0' && content[i] <= '9')) {
						i = RestoreI; fsc = tfntFetchCharTextureIDW(content[i], 0); break;
					}
					id = id * 10 + content[i] - '0';
					i++;
				}				
				if (i != RestoreI) fsc = tfntFetchCharTextureIDW(id, 1);
			}

			sx += fsc->advx + ts->TextSpacing;
			if ((content[i]) > 127)C += 2; else C += 1;
			if (Count && C == Count)
				return sx;
		}
	}

	return sx;
}
int tnsStringGetWidthM(char* content, nThemeState* ts, int Count){
	wchar_t tempdata[512] = { 0 };
	if (!content) return 0;
	MultiByteToWideChar(CP_ACP, 0, content, strlen(content), tempdata, 512);
	return tnsStringGetWidthU(tempdata,ts, Count);
}
//int tnsCalcSingleLineWidth(wchar_t* content){
//	int revel = 0;
//	wchar_t* p = content;
//	while (*p != L'\n' && *p != L'\0'){
//		tnsFontSingleCharacter* fsc;
//		if (*p != '$') {
//			fsc = tfntFetchCharTextureIDW(*p, 0);
//		}else {
//			int id;
//			vswscanf(content, "%d", &id);
//			fsc = tfntFetchCharTextureIDW(id, 1);
//			p++;
//			while (*p != '$') p++;
//		}
//		revel += fsc->advx;
//		p++;
//	}
//	return revel;
//}
void tnsDrawVectorGraphPackage(int ID, nThemeState* ts, int L, int R, int T, int B, int Align, int RevY) {
	int sx = L, sy;// = (FM->VectorsGrapghs->height + T);
	int i;
	//int line_width = tnsCalcSingleLineWidth(content);
	real xo, yo, xl, yl;
	real TexCoord[8];
	real VertArr[8];
	int total_width = 0;
	int Index[4] = { 0, 1, 2, 3 };
	tnsFont* f = FM->UsingFont;
	int FscHeight;

	tnsFontSingleCharacter* fsc = tfntFetchVectorGraphTextureIDW(ID);

	int Size = fsc->height;
	sy = Size + T;

	int cx = sx + fsc->deltax;
	int cy = sy;// -fsc->deltay;
	
	tnsMakeQuad2d(TexCoord, 0, 1, 0, 0, 1, 1, 1, 0);

	//if (RevY)
	//	tnsMakeQuad2d(VertArr,
	//		cx, 2 * T - cy + fsc->height,
	//		cx, 2 * T - cy,
	//		cx + fsc->width, 2 * T - cy + fsc->height,
	//		cx + fsc->width, 2 * T - cy);
	//else


	tnsMakeQuad2d(VertArr,
		cx, +cy - fsc->height,
		cx, +cy,
		cx + fsc->width, +cy - fsc->height,
		cx + fsc->width, +cy);

	//SWITCH SHADER
	tnsUseTexture0(fsc->Tex);
	tnsColor4dv(ts->TextColor->RGBA);
	tnsVertexArray2d(VertArr, 4);
	tnsTexCoordArray2d(TexCoord, 4);
	tnsPackAs(GL_TRIANGLE_STRIP);
}
void tnsDrawStringDirectU(wchar_t* content,nThemeState* ts,int L,int R,int T,int RevY){
	int sx = L, sy = (FM->UsingFont->height + T);
	int i;
	int len = wcslen(content);
	//int line_width = tnsCalcSingleLineWidth(content);
	real xo, yo, xl, yl;
	real TexCoord[8];
	real VertArr[8];
	int total_width = 0;
	tnsFont* f = FM->UsingFont;
	int FscHeight;
	int RestoreI;


	for (i = 0; i<len; i++){
		tnsFontSingleCharacter* fsc;
		int cx, cy;
		
		if (content[i] == '\n'){
			continue;
		}

		if (content[i] != '$') {
			fsc = tfntFetchCharTextureIDW(content[i], 0);
		}else {
			int id=0;
			RestoreI = i;
			i++;
			while (content[i] != '$') {
				if (!(content[i] >= '0' && content[i] <= '9')) {
					i = RestoreI; fsc = tfntFetchCharTextureIDW(content[i], 0); break;
				}
				id = id * 10 + content[i]-'0';
				i++;
			}
			if(i!=RestoreI) fsc = tfntFetchCharTextureIDW(id, 1);
		}

		if (sx + fsc->advx /*+ ts->TextSpacing */> R) break;

		cx = sx + fsc->deltax;
		cy = sy - fsc->deltay;
		xo = (real)fsc->bufferx / (real)TNS_FONT_BUFFER_W;
		yo = (real)fsc->buffery / (real)TNS_FONT_BUFFER_H;
		xl = (real)(fsc->bufferx + fsc->width) / (real)TNS_FONT_BUFFER_W;
		yl = (real)(fsc->buffery + fsc->height) / (real)TNS_FONT_BUFFER_H;

		tnsMakeQuad2d(TexCoord, xo, yl, xo, yo, xl, yl, xl, yo);
		if (RevY)
			tnsMakeQuad2d(VertArr,
				cx, 2 * T - cy + fsc->height,
				cx, 2 * T - cy,
				cx + fsc->width, 2 * T - cy + fsc->height,
				cx + fsc->width, 2 * T - cy);
		else
			tnsMakeQuad2d(VertArr,
				cx, +cy - fsc->height,
				cx, +cy,
				cx + fsc->width, +cy - fsc->height,
				cx + fsc->width, +cy);

		//SWITCH SHADER
		tnsUseTexture0(&f->TexBuffer);
		tnsColor4dv(ts->TextColor->RGBA);
		tnsVertexArray2d(VertArr, 4);
		tnsTexCoordArray2d(TexCoord, 4);
		tnsPackAs(GL_TRIANGLE_STRIP);

		sx += fsc->advx + ts->TextSpacing;
		if (sx > R) break;
	}
}
void tnsDrawStringAutoM(char* content, nThemeState* ts, int L, int R, int T, int RevY, nStringSplitor* Instructions) {
	char buf[512] = { 0 };
	int LL;
	int al = Instructions ? (strArgumentMatch(Instructions, "text_align", "left") ? NUL_TEXT_ALIGN_LEFT :
							(strArgumentMatch(Instructions, "text_align", "right") ? NUL_TEXT_ALIGN_RIGHT :
							(strArgumentMatch(Instructions, "text_align", "center") ? NUL_TEXT_ALIGN_CENTER :
							ts->TextAlign))) : ts->TextAlign;
	if (!content) return;
	switch (al) {
	case NUL_TEXT_ALIGN_AUTO:
	case NUL_TEXT_ALIGN_CENTER:
		LL = L + (R - L - tnsStringGetWidthM(content, ts, 0)) / 2;
		break;
	case NUL_TEXT_ALIGN_LEFT:
		LL = L;
		break;
	case NUL_TEXT_ALIGN_RIGHT:
		LL = R - (tnsStringGetWidthM(content, ts, 0));
		break;
	}
	strToWideChar(buf, content);
	tnsDrawStringDirectU(buf, ts, LL, R, T, RevY);
}
void tnsDrawStringAutoU(wchar_t* content, nThemeState* ts, int L, int R, int T, int RevY, nStringSplitor* Instructions) {
	int LL;
	int al = Instructions ? (strArgumentMatch(Instructions, "text_align", "left") ? NUL_TEXT_ALIGN_LEFT :
							(strArgumentMatch(Instructions, "text_align", "right") ? NUL_TEXT_ALIGN_RIGHT :
							(strArgumentMatch(Instructions, "text_align", "center") ? NUL_TEXT_ALIGN_CENTER :
							ts->TextAlign))) : ts->TextAlign;
	if (!content) return;
	switch (al) {
	case NUL_TEXT_ALIGN_AUTO:
	case NUL_TEXT_ALIGN_CENTER:
		LL = L+(R-L-tnsStringGetWidthU(content,ts, 0)) / 2;
		break;
	case NUL_TEXT_ALIGN_LEFT:
		LL = L;
		break;
	case NUL_TEXT_ALIGN_RIGHT:
		LL = R - (tnsStringGetWidthU(content,ts, 0));
		break;
	}
	tnsDrawStringDirectU(content, ts, LL, R, T, RevY);
}
void tnsDrawStringWithPriority(wchar_t* Label, wchar_t* MajorContent,nThemeState* ts, int L, int R, int T, int RevY) {
	int W = R - L;
	int Str1W = tnsStringGetWidthU(Label, ts, 0);
	int Str2W = tnsStringGetWidthU(MajorContent, ts, 0);
	if (Str1W + Str2W > W) {
		if (Str2W < W) {
			tnsDrawStringDirectU(MajorContent, ts, R - Str2W, R, T, RevY);
			tnsDrawStringDirectU(Label, ts, L, R - Str2W, T, RevY);
		}
		else
			tnsDrawStringDirectU(MajorContent, ts, L, R, T, RevY);
	}else {
		int LL = L, ML;
		switch (ts->TextAlign) {
		case NUL_TEXT_ALIGN_CENTER:
			ML = L+Str1W+(W - (Str1W + Str2W)) / 2;
			break;
		case NUL_TEXT_ALIGN_LEFT:
			ML = L + Str1W;
			break;
		case NUL_TEXT_ALIGN_AUTO:
		case NUL_TEXT_ALIGN_RIGHT:
			ML = R - (Str2W);
			break;
		}

		tnsDrawStringDirectU(Label, ts, LL, R, T, RevY);
		tnsDrawStringDirectU(MajorContent, ts, ML, R, T, RevY);
	}
}

void tnsDrawIcon(short ID, nThemeState* ts, int L, int T, int RevY) {
	int sx = L, sy = (FM->UsingFont->height + T);
	int i;
	//int line_width = tnsCalcSingleLineWidth(content);
	real xo, yo, xl, yl;
	real TexCoord[8];
	real VertArr[8];
	int total_width = 0;
	tnsFont* f = FM->UsingFont;
	int FscHeight;


	tnsFontSingleCharacter* fsc = tfntFetchCharacterW(ID, 1);
	int cx, cy;

	cx = sx + fsc->deltax;
	cy = sy - fsc->deltay;
	xo = (real)fsc->bufferx / (real)TNS_FONT_BUFFER_W;
	yo = (real)fsc->buffery / (real)TNS_FONT_BUFFER_H;
	xl = (real)(fsc->bufferx + fsc->width) / (real)TNS_FONT_BUFFER_W;
	yl = (real)(fsc->buffery + fsc->height) / (real)TNS_FONT_BUFFER_H;

	tnsMakeQuad2d(TexCoord, xo, yl, xo, yo, xl, yl, xl, yo);
	if (RevY)
		tnsMakeQuad2d(VertArr,
			cx, 2 * T - cy + fsc->height,
			cx, 2 * T - cy,
			cx + fsc->width, 2 * T - cy + fsc->height,
			cx + fsc->width, 2 * T - cy);
	else
		tnsMakeQuad2d(VertArr,
			cx, +cy - fsc->height,
			cx, +cy,
			cx + fsc->width, +cy - fsc->height,
			cx + fsc->width, +cy);

	//SWITCH SHADER
	tnsUseTexture0(&f->TexBuffer);
	tnsColor4dv(ts->TextColor->RGBA);
	tnsVertexArray2d(VertArr, 4);
	tnsTexCoordArray2d(TexCoord, 4);
	tnsPackAs(GL_TRIANGLE_STRIP);

}




//=====================================================================================[Object]

tns3DObject* tnsFindRootObject(tnsScene* s,char* Name) {
	tns3DObject* o;
	for (o = s->AllObjects.pFirst; o; o = o->Item.pNext) {
		if (strIsTheSame(o->Name->Ptr, Name)) return o;
	}
	return 0;
}
tns3DObject* tnsFindObjectRecursive(tns3DObject* o, char* Name, tns3DObject** From) {
	tns3DObject* io;
	tns3DObject* ro;
	for (io = o->ChildObjects.pFirst; io; io = io->Item.pNext) {
		if (strIsTheSame(io->Name->Ptr, Name)) {
			*From = &io;
			return o;
		}
		if (ro = tnsFindObjectRecursive(io, Name, From))return ro;
	}
	return 0;
}
tns3DObject* tnsFindObject(tnsScene* s, char* Name, tns3DObject** From) {
	tns3DObject* o;
	tns3DObject* ro;
	for (o = s->Objects.pFirst; o; o = o->Item.pNext) {
		if (strIsTheSame(o->Name->Ptr, Name)) {
			if(From)*From = 0;
			return o;
		}
		if (ro = tnsFindObjectRecursive(o, Name, From))return ro;
	}
	return 0;
}

void tObjExtractSelfEulerRotation(tns3DObject* o) {
	switch (o->RotationMode) {
	case TNS_ROTATION_XYZ_EULER:
		tMatExtractXYZEuler44d(o->SelfTransform,&o->Rotation[0],&o->Rotation[1],&o->Rotation[2]);
		break;
	}
}
void tObjExtractSelfLocation(tns3DObject* o) {
	tMatExtractLocation44d(o->SelfTransform, &o->Location[0], &o->Location[1], &o->Location[2]);
}
void tObjExtractSelfScaling(tns3DObject* o) {
	tMatExtractUniformScale44d(o->SelfTransform, &o->Scale);
}
void tObjApplySelfTransformMatrix(tns3DObject* o, real* inverse_result);
void tObjApplyGlobalTransformMatrix(tns3DObject* o) {
	if (!o->ParentObject) {
		memcpy(o->GlobalTransform, o->SelfTransform, sizeof(tnsMatrix44d));
	}
	else {
		tMatMultiply44d(o->GlobalTransform, o->ParentObject->GlobalTransform, o->SelfTransform);
	}
	tMatExtractLocation44d(o->GlobalTransform, &o->GLocation[0], &o->GLocation[1], &o->GLocation[2]);
	tMatExtractXYZEuler44d(o->GlobalTransform, &o->GRotation[0], &o->GRotation[1], &o->GRotation[2]);
	tMatExtractUniformScale44d(o->GlobalTransform, &o->GScale);
}
void tObjApplyGlobalAndSelfTransformMatrix(tns3DObject* o) {
	if (!o->ParentObject) {
		memcpy(o->Location, o->GLocation, sizeof(tnsVector3d));
		memcpy(o->Rotation, o->GRotation, sizeof(tnsVector3d));
		o->Scale = o->GScale;
		tObjApplySelfTransformMatrix(o, 0);
	}
	else {
		tMatMultiply44d(o->GlobalTransform, o->ParentObject->GlobalTransform, o->SelfTransform);
	}
	tMatExtractLocation44d(o->GlobalTransform, &o->GLocation[0], &o->GLocation[1], &o->GLocation[2]);
	tMatExtractXYZEuler44d(o->GlobalTransform, &o->GRotation[0], &o->GRotation[1], &o->GRotation[2]);
	tMatExtractUniformScale44d(o->GlobalTransform, &o->GScale);
}
void tObjApplyGlobalTransformMatrixRecursive(tns3DObject* o) {
	tns3DObject* oc;
	if (!o->ParentObject) {
		memcpy(o->GlobalTransform, o->SelfTransform, sizeof(tnsMatrix44d));
	}
	else {
		tMatMultiply44d(o->GlobalTransform, o->ParentObject->GlobalTransform, o->SelfTransform);
	}
	if (o->ChildObjects.pFirst) {
		for (oc = o->ChildObjects.pFirst; oc; oc = oc->Item.pNext) {
			tObjApplyGlobalTransformMatrixRecursive(oc);
		}
	}
	tMatExtractLocation44d(o->GlobalTransform, &o->GLocation[0], &o->GLocation[1], &o->GLocation[2]);
	tMatExtractXYZEuler44d(o->GlobalTransform, &o->GRotation[0], &o->GRotation[1], &o->GRotation[2]);
	tMatExtractUniformScale44d(o->GlobalTransform, &o->GScale);
}
void tObjApplyGlobalTransformMatrixReverted(tns3DObject* o) {
	if (!o->ParentObject) {
		memcpy(o->GlobalTransform, o->SelfTransform, sizeof(tnsMatrix44d));
	}
	else {
		tObjApplyGlobalTransformMatrixReverted(o->ParentObject);
		tMatMultiply44d(o->GlobalTransform, o->ParentObject->GlobalTransform, o->SelfTransform);
	}
	tMatExtractLocation44d(o->GlobalTransform, &o->GLocation[0], &o->GLocation[1], &o->GLocation[2]);
	tMatExtractXYZEuler44d(o->GlobalTransform, &o->GRotation[0], &o->GRotation[1], &o->GRotation[2]);
	tMatExtractUniformScale44d(o->GlobalTransform, &o->GScale);
}
void tObjApplySelfTransformMatrix(tns3DObject* o, real* inverse_result) {
	tnsMatrix44d Trans, Rot1, Rot2, Rot3, Scale, Res1, Res2;
	real LX, LY, LZ, RX, RY, RZ, SC;
	if (!o) return;
	else {

		LX = o->Location[0];
		LY = o->Location[1];
		LZ = o->Location[2];
		RX = o->Rotation[0];
		RY = o->Rotation[1];
		RZ = o->Rotation[2];
		SC = o->Scale;

		tMatLoadIdentity44d(o->SelfTransform);
		tMatMakeTranslationMatrix44d(Trans, LX, LY, LZ);
		tMatMakeScaleMatrix44d(Scale, SC, SC, SC);

		tMatMakeRotationXMatrix44d(Rot1, RX);
		tMatMakeRotationYMatrix44d(Rot2, RY);
		tMatMakeRotationZMatrix44d(Rot3, RZ);
		switch (o->RotationMode) {
		case TNS_ROTATION_ZYX_EULER:
			tMatMultiply44d(Res1, Trans, Rot1);
			tMatMultiply44d(Res2, Res1, Rot2);
			tMatMultiply44d(Res1, Res2, Rot3);
			tMatMultiply44d(o->SelfTransform, Res1, Scale);
			break;
		case TNS_ROTATION_XZY_EULER:
			tMatMultiply44d(Res1, Trans, Rot1);
			tMatMultiply44d(Res2, Res1, Rot3);
			tMatMultiply44d(Res1, Res2, Rot2);
			tMatMultiply44d(o->SelfTransform, Res1, Scale);
			break;
		case TNS_ROTATION_YXZ_EULER:
			tMatMultiply44d(Res1, Trans, Rot2);
			tMatMultiply44d(Res2, Res1, Rot1);
			tMatMultiply44d(Res1, Res2, Rot3);
			tMatMultiply44d(o->SelfTransform, Res1, Scale);
			break;
		case TNS_ROTATION_YZX_EULER:
			tMatMultiply44d(Res1, Trans, Rot2);
			tMatMultiply44d(Res2, Res1, Rot3);
			tMatMultiply44d(Res1, Res2, Rot1);
			tMatMultiply44d(o->SelfTransform, Res1, Scale);
			break;
		case TNS_ROTATION_ZXY_EULER:
			tMatMultiply44d(Res1, Trans, Rot3);
			tMatMultiply44d(Res2, Res1, Rot1);
			tMatMultiply44d(Res1, Res2, Rot2);
			tMatMultiply44d(o->SelfTransform, Res1, Scale);
			break;
		case TNS_ROTATION_XYZ_EULER:
			tMatMultiply44d(Res1, Trans, Rot3);
			tMatMultiply44d(Res2, Res1, Rot2);
			tMatMultiply44d(Res1, Res2, Rot1);
			tMatMultiply44d(o->SelfTransform, Res1, Scale);
			break;
		}
	}

	tObjApplyGlobalTransformMatrix(o);

	if (o->ChildObjects.pFirst) {
		tns3DObject* io;
		for (io = o->ChildObjects.pFirst; io; io = io->Item.pNext) {
			tObjApplySelfTransformMatrix(io,0);
		}
	}
}
void tObjApplySelfTransformValue(tns3DObject* o) {
	if (!o) return;
	
	tMatExtractLocation44d(o->SelfTransform, &o->Location[0], &o->Location[1], &o->Location[2]);
	tMatExtractXYZEuler44d(o->SelfTransform, &o->Rotation[0], &o->Rotation[1], &o->Rotation[2]);
	tMatExtractUniformScale44d(o->SelfTransform, &o->Scale);
	
}
void tObjApplyGlobalTransformValue(tns3DObject* o) {
	if (!o) return;

	tMatExtractLocation44d(o->GlobalTransform, &o->GLocation[0], &o->GLocation[1], &o->GLocation[2]);
	tMatExtractXYZEuler44d(o->GlobalTransform, &o->GRotation[0], &o->GRotation[1], &o->GRotation[2]);
	tMatExtractUniformScale44d(o->GlobalTransform, &o->GScale);

}
void tObjApplySameGlobalTransform(tns3DObject* to, tns3DObject* from) {
	memcpy(to->GlobalTransform, from->GlobalTransform, sizeof(tnsMatrix44d));
	memcpy(to->GLocation, from->GLocation, sizeof(tnsVector3d));
	memcpy(to->GRotation, from->GRotation, sizeof(tnsVector3d));
	to->GScale = from->GScale;

	tObjApplyGlobalTransformValue(to);
	tObjGetSelfTransformFromGlobal(to);
	tObjApplySelfTransformValue(to);
}
void tObjGetSelfTransformFromGlobal(tns3DObject* o) {
	tnsMatrix44d LastGInv;
	if (o->ParentObject) {
		tMatInverse44d(LastGInv, o->ParentObject->GlobalTransform);
		tMatMultiply44d(o->SelfTransform, LastGInv, o->GlobalTransform);
		tObjApplySelfTransformValue(o);
		tObjApplyGlobalTransformValue(o);
	}else{
		tMatCopyMatrix44d(o->GlobalTransform, o->SelfTransform);
		tMatVectorCopy4d(o->GRotation, o->Rotation);
		tMatVectorCopy3d(o->GLocation, o->Location);
		o->Scale = o->GScale;
	}
}
void tObjInitObjectBase(tns3DObject* o, tnsScene* Scene, char* Name, int Type,
	real AtX, real AtY, real AtZ,
	real RotX, real RotY, real RotZ, real RotW, u8bit RotationMode,
	real Scale) {
	if (!o) return;
	strSafeSet(&o->Name, Name);
	o->Type = Type;
	o->InScene = Scene;
	o->Location[0] = AtX;
	o->Location[1] = AtY;
	o->Location[2] = AtZ;
	o->Rotation[0] = RotX;
	o->Rotation[1] = RotY;
	o->Rotation[2] = RotZ;
	o->Rotation[3] = RotW;
	o->RotationMode = RotationMode;
	o->Scale = Scale;
	o->Show = 1;
	o->DrawMode = GL_LINE_LOOP;
	tObjApplySelfTransformMatrix(o, 0);
	if (Scene) {
		lstAppendItem2(&Scene->AllObjects, o);
		lstAppendItem(&Scene->Objects, o);
	}
}
void tnsParentObject(tns3DObject* child, tns3DObject* parent) {
	if (!child || !parent || lstFindItem(child, nutSameAddress, &parent->ChildObjects)) return;

	if (parent->Type == TNS_OBJECT_CAMERA) return;

	if (child->ParentObject) lstRemoveItem(&child->ParentObject->ChildObjects, child);
	else switch (child->Type) {
	case TNS_OBJECT_PLACEHOLDER:
	case TNS_OBJECT_LAMP:
	case TNS_OBJECT_MESH:
	case TNS_OBJECT_CAMERA:
		lstRemoveItem(&child->InScene->Objects, child);
		break;
	}

	lstAppendItem(&parent->ChildObjects, child);
	child->ParentObject = parent;
	tObjApplyGlobalTransformMatrix(child);
}
void tnsUnparentObject(tns3DObject* o) {
	if (!o || !o->ParentObject) return;
	lstRemoveItem(&o->ParentObject->ChildObjects, o);
	switch (o->Type) {
	case TNS_OBJECT_PLACEHOLDER:
	case TNS_OBJECT_MESH:
	case TNS_OBJECT_CAMERA:
	case TNS_OBJECT_LAMP:
		lstAppendItem(&o->InScene->Objects, o);
		break;
	}
	tObjApplyGlobalTransformMatrix(o);
}
void tnsRotateObjectSelf(tns3DObject* o, real x, real y, real z) {
	o->Rotation[0] += x;
	o->Rotation[1] += y;
	o->Rotation[2] += z;
	tObjApplySelfTransformMatrix(o, 0);
}
void tnsMoveObjectPrior(tns3DObject* o, real x, real y, real z) {
	o->Location[0] += x;
	o->Location[1] += y;
	o->Location[2] += z;
	tObjApplySelfTransformMatrix(o, 0);
}
void tnsRotateObjectGlobalXYZ(tns3DObject* o, real x, real y, real z) {
	tnsMatrix44d mat, res1, res2;
	real xs, ys, zs;

	//save location data;

	xs = o->GlobalTransform[12];
	ys = o->GlobalTransform[13];
	zs = o->GlobalTransform[15];
	
	tMatLoadIdentity44d(res1);

	tMatMakeRotationZMatrix44d(mat, z);
	tMatMultiply44d(res2, res1, mat);
	tMatMakeRotationYMatrix44d(mat, y);
	tMatMultiply44d(res1, res2, mat);
	tMatMakeRotationXMatrix44d(mat, x);
	tMatMultiply44d(res2, res1, mat);

	tMatMultiply44d(res1, res2,o->GlobalTransform);
	memcpy(o->GlobalTransform, res1, sizeof(tnsMatrix44d));

	o->GlobalTransform[12] = xs;
	o->GlobalTransform[13] = ys;
	o->GlobalTransform[14] = zs;

	//tObjApplyGlobalTransformValue(o);
	tObjGetSelfTransformFromGlobal(o);
	tObjApplySelfTransformValue(o);
	tObjApplyGlobalTransformMatrix(o);
}
void tnsRotateObjectGlobalBaseXYZ(tns3DObject* o, tnsVector3d BaseLocation, real x, real y, real z) {
	tnsMatrix44d mat, res1, res2,t1,t2;
	tnsVector4d orig = { 0,0,0,1 }, trans = {0,0,0,1}, oo = { 0,0,0,1 },trans2;
	real xs, ys, zs;

	//save location data;

	

	tMatVectorMinus3d(trans, o->GLocation, BaseLocation);

	xs = o->GlobalTransform[12];
	ys = o->GlobalTransform[13];
	zs = o->GlobalTransform[14];

	tMatLoadIdentity44d(res1);
	tMatMakeRotationZMatrix44d(mat, z);
	tMatMultiply44d(res2, res1, mat);
	tMatMakeRotationYMatrix44d(mat, y);
	tMatMultiply44d(res1, res2, mat);
	tMatMakeRotationXMatrix44d(mat, x);
	tMatMultiply44d(res2, res1, mat);
	tMatMultiply44d(res1, res2, o->GlobalTransform);

	tMatApplyTransform44d(trans2, res2, trans);

	memcpy(o->GlobalTransform, res1, sizeof(tnsMatrix44d));

	o->GlobalTransform[12] = trans2[0] + BaseLocation[0];
	o->GlobalTransform[13] = trans2[1] + BaseLocation[1];
	o->GlobalTransform[14] = trans2[2] + BaseLocation[2];
	
	//tObjApplyGlobalTransformValue(o);
	tObjGetSelfTransformFromGlobal(o);
	tObjApplySelfTransformValue(o);
	//tObjApplyGlobalTransformMatrix(o);//?necessary?
}
void tnsTranslateObjectGlobal(tns3DObject* o, real x, real y, real z) {
	tnsMatrix44d mat,res1;


	tMatLoadIdentity44d(res1);

	tMatMakeTranslationMatrix44d(mat, x, y, z);

	tMatMultiply44d(res1, mat, o->GlobalTransform);
	memcpy(o->GlobalTransform, res1, sizeof(tnsMatrix44d));


	//tObjApplyGlobalTransformValue(o);
	tObjGetSelfTransformFromGlobal(o);
	tObjApplySelfTransformValue(o);
	tObjApplyGlobalTransformMatrix(o);
}
void tnsTranslateObjectLocal(tns3DObject* o, real x, real y, real z) {
	tnsMatrix44d mat, res1;


	tMatLoadIdentity44d(res1);

	tMatMakeTranslationMatrix44d(mat, x, y, z);

	tMatMultiply44d(res1, o->GlobalTransform, mat);
	memcpy(o->GlobalTransform, res1, sizeof(tnsMatrix44d));


	//tObjApplyGlobalTransformValue(o);
	tObjGetSelfTransformFromGlobal(o);
	tObjApplySelfTransformValue(o);
	tObjApplyGlobalTransformMatrix(o);
}

void tnsMoveObjectGlobal(tns3DObject* o, real x, real y, real z) {
	tnsMatrix44d mat, res1, res2, res3;

	tMatMakeTranslationMatrix44d(res2, x, y, z);

	tMatMultiply44d(res1, res2, o->SelfTransform);//if reverse then local
	memcpy(o->SelfTransform, res1, sizeof(tnsMatrix44d));
	tObjApplySelfTransformValue(o);
	tObjApplyGlobalTransformMatrix(o);
}
void tnsScaleObject(tns3DObject* o, real fac) {
	o->Scale *= fac;
	tObjApplySelfTransformMatrix(o, 0);
	tObjApplyGlobalTransformMatrix(o);
}
void tnsSetObjectLocation(tns3DObject* o, real x, real y, real z) {
	o->Location[0] = x;
	o->Location[1] = y;
	o->Location[2] = z;
	tObjApplySelfTransformMatrix(o, 0);
	tObjApplyGlobalTransformMatrix(o);
}

void tnsZoomViewingCamera(tnsCamera* c, real Ratio) {
	if (c->FocusDistance < 0.1)return;
	tnsTranslateObjectLocal(c, 0, 0, -c->FocusDistance*Ratio);
	c->FocusDistance *= (1-Ratio);
}
void tnsRotateViewingCamera(tnsCamera* c, real x,real z) {
	tnsTranslateObjectLocal(c, 0, 0, -c->FocusDistance);
	tnsRotateObjectSelf(c, x, 0, z);
	tnsTranslateObjectLocal(c, 0, 0, c->FocusDistance);
}
void tnsTranslateViewingCamera(tnsCamera* c, int ViewportW,int ViewportH,real x, real y) {
	tnsVector4d p = { 0 }, tp = { 0 };
	tnsMatrix44d combine,projection,projinv,vpinv,vp;

	tMatLoadIdentity44d(projection);
	tMatMakePerspectiveMatrix44d(projection, c->FOV, (real)ViewportW / (real)ViewportH, c->ZMin, c->ZMax);
	tMatMakeViewportMatrix44d(vp, ViewportW, ViewportH,c->ZMax,c->ZMin);
	tMatInverse44d(projinv, projection);
	tMatInverse44d(vpinv, vp);
	//tMatInverse44d(inv, projection);

	//p[0] = x;
	//p[1] = y;
	p[2] = -c->FocusDistance;
	p[3] = 1;
	tMatApplyTransform44d(tp, projection, p);
	tMatVectorMultiSelf3d(tp, 1 / tp[3]);
	tMatApplyTransform43d(p, vp, tp);
	p[0] += x;
	p[1] += y;

	tMatApplyTransform43d(tp, vpinv, p);
	tMatVectorMultiSelf3d(tp, tp[3]);
	tMatApplyTransform44d(p, projinv, tp);
	//printf("%d,%d\n", tp[0]?x/tp[0]:0, tp[1]?y/tp[1]:0);
	//tnsTranslateObjectLocal(c, tp[0] ? fabs(tp[0])/x : 0, tp[1] ? fabs(tp[1])/y : 0, 0);

	tnsTranslateObjectLocal(c, p[0], p[1], 0);
}

tnsScene* tnsCreateScene(char* name) {
	tnsWorld* w = &T->World;
	tnsScene* nc = memAquireHyper(sizeof(tnsScene));

	strSafeSet(&nc->Name, name);
	lstAppendItem(&w->Scenes, nc);

	w->ActiveScene = nc;

	tns_InitAtlasInput(nc, TNS_ATLAS_DEFAULT_INPUT_WIDTH);

	return nc;
}
void tnsDestroyScene(tnsScene* s) {
	tnsWorld* w = &T->World;
	tnsScene* nc = s;
	tns3DObject* o,*NextO;

	tns_DestroyAtlasInputs(s);
	
	w->ActiveScene = nc->Item.pPrev ? nc->Item.pPrev : nc->Item.pNext ? nc->Item.pNext : 0;

	lstRemoveItem(&w->Scenes, nc);

	while( o = lstPopItem2(&s->AllObjects)) {
		
		tnsDestroyObjectRecursive(o, s);
	}
	strSafeDestroy(&s->Name);
	memFree(s);
}
void tnsActiveScene(tnsScene* s) {
	tnsWorld* w = &T->World;
	if (!s) return;
	w->ActiveScene = s;
}
tnsCamera* tnsCreateXYZCamera(tnsScene* Scene, char* Name, real FOV,
	real AtX, real AtY, real AtZ,
	real RotX, real RotY, real RotZ,
	real FocusDistance) {
	tnsCamera* c;
	tnsWorld* w = &T->World;

	//if (!Scene) return 0;

	c = memAquireHyper(sizeof(tnsCamera));
	tObjInitObjectBase(&c->Base, Scene, Name, TNS_OBJECT_CAMERA,
		AtX, AtY, AtZ, RotX, RotY, RotZ, 1.0f, TNS_ROTATION_XYZ_EULER, 1.0f);
	c->FOV = FOV;
	c->CameraType = TNS_PRESPECTIVE_CAMERA;
	c->ZMin = 0.1f;
	c->ZMax = 1000.0f;
	c->FocusDistance = FocusDistance;
	c->OrthScale = 1.0f;

	if (Scene && !Scene->ActiveCamera)Scene->ActiveCamera = c;
		
	return c;
}
tns3DObject* tnsCreateXYZEmpty(tnsScene* Scene, char* Name, real AtX, real AtY, real AtZ) {
	tns3DObject* o;
	tnsWorld* w = &T->World;

	if (!Scene) return 0;

	o = memAquireHyper(sizeof(tns3DObject));
	tObjInitObjectBase(o, Scene, Name, TNS_OBJECT_PLACEHOLDER,
		AtX, AtY, AtZ, 0, 0, 0, 1.0f, TNS_ROTATION_XYZ_EULER, 1.0f);

	return o;
}
void tnsApplyCameraView(int W,int H, tnsCamera* Camera) {
	tnsMatrixStackItem* tmsi = tKnlGetCurrentMatStackItem();
	tnsShader* current_shader = 0;
	real* mat;
	tnsMatrix44d result,inv;

	if (!Camera) return;


	mat = tnsGetProjectionMatrix();
	tMatMakePerspectiveMatrix44d(mat, Camera->FOV, (real)W / (real)H, Camera->ZMin, Camera->ZMax);
    if (current_shader = tKnlGetActiveShader()) {
		tSdrApplyProjection(current_shader, mat);
		tSdrApplyProjectionInverse(current_shader, mat);
	}

	tnsResetViewMatrix();
	mat = tnsGetViewMatrix();
	tObjApplySelfTransformMatrix(Camera,0);
	tObjApplyGlobalTransformMatrix(Camera);
	tMatInverse44d(inv, Camera->Base.GlobalTransform);
	tMatMultiply44d(result, mat, inv);
	memcpy(mat, result, sizeof(tnsMatrix44d));

	if (current_shader = tKnlGetActiveShader()) {
		tSdrApplyView(current_shader, result);
	}
}

void tnsSetActiveCamera(tnsScene* s, tnsCamera* Camera) {
	if (!s || !Camera) return 0;
	if (!tnsFindObject(s, Camera->Base.Name->Ptr, 0))return 0;
	s->ActiveCamera = Camera;
}


void tnsApplyObjectMatrix(tns3DObject* Object) {
	tnsShader* current_shader = 0;
	real* mat;
	tnsMatrix44d result;

	if (!Object) return;

	mat = tnsGetModelMatrix();

	tMatMultiply44d(result, mat, Object->SelfTransform);
	memcpy(mat, result, sizeof(tnsMatrix44d));

	//Actually This Works Pretty Fast,but we are currently testing its functionality;
	//memcpy(mat, Object->GlobalTransform, sizeof(tnsMatrix44d));
	
	if (current_shader = tKnlGetActiveShader()) {
		tSdrApplyModel(current_shader, mat);
	}
}


void tnsDestroyObjectRecursive(tns3DObject* Object,tnsScene* Scene) {
	tns3DObject* o,*NextO;
	tnsUnparentObject(Object);
	lstRemovePointer(&Scene->AllObjects, Object);
	lstRemoveItem(&Scene->Objects, Object);
	for (o = Object->ChildObjects.pFirst; o; o = NextO) {
		NextO = o->Item.pNext;
		tnsDestroyObjectRecursive(o, Scene);
	}
}

void tnsTrackZTo(tns3DObject* o,tnsVector3d Target,tnsVector3d Up) {
	tnsVector4d target = { 0 };
	tnsMatrix44d mat, res,res2;
	real* d;
	tnsShader* current_shader = 0;

	d = tnsGetModelMatrix();

	tMatMakeZTrackingMatrix44d(mat, o->Location, target, Up);
	tMatMultiply44d(res, o->SelfTransform,mat);
	tMatMultiply44d(res2, res, d);
	memcpy(d, res2, sizeof(tnsMatrix44d));

	//memcpy(o->SelfTransform, res, sizeof(tnsMatrix44d));

	//d = tKnlGetModifyingMatirx(tKnlGetCurrentMatStackItem());
	//memcpy(d, res, sizeof(tnsMatrix44d));

	if (current_shader = tKnlGetActiveShader()) {
		tSdrApplyModel(current_shader, d);
	}
}




void tnsDrawCamera(tns3DObject*o, tnsScene* s) {
	tnsCamera* c = o;
	real fov_2 = c->FOV/2;
	real ex, ey, ez;
	ey = 10 * sin(fov_2);
	ex = ey;
	ez = 10 * cos(fov_2);

	//if (T->CurrentShader != T->uiShader) tnsEnableShaderv(T->uiShader);

	tnsColor4d(1, 1, 1, 1);

	tnsVertex3d(ex, ey, -ez);
	tnsVertex3d(ex, -ey, -ez);
	tnsVertex3d(0.0, 0.0, 0.0);
	tnsVertex3d(ex, ey, -ez);
	tnsVertex3d(-ex, ey, -ez);
	tnsVertex3d(0.0, 0.0, 0.0);
	tnsVertex3d(-ex, -ey, -ez);
	tnsVertex3d(-ex, ey, -ez);

	tnsPackAs(GL_LINE_STRIP);
	tnsFlush();
}
void tnsDrawThisObject(tns3DObject* o,tnsScene* s) {
	if (!o->Show) return;
	switch (o->Type) {
	case TNS_OBJECT_MESH:
		//if (T->CurrentShader != T->TEST_MatcapShader) tnsEnableShaderv(T->TEST_MatcapShader);
		tnsDrawMeshObject(o, 0);
		break;
	case TNS_OBJECT_CAMERA:
		//tnsDrawCamera(o, s);
		break;
	case TNS_OBJECT_BEZIER_CURVE:
		//tnsDrawBezierObject(o, 0);
		break;
	case TNS_OBJECT_PLACEHOLDER:
	default:
		////if (T->CurrentShader != T->uiShader) tnsEnableShaderv(T->uiShader);
		//tnsVertex3d(-1.0, 0.0, 0.0);
		//tnsVertex3d(5.0, 0.0, 0.0);
		//tnsVertex3d(0.0, 5.0, 0.0);
		//tnsVertex3d(0.0, -1.0, 0.0);
		//tnsVertex3d(0.0, 0.0, -1.0);
		//tnsVertex3d(0.0, 0.0, 5.0);
		//tnsColor4d(1, 1, 1, 1);
		//tnsPackAs(GL_LINES);
		//tnsFlush();
		break;
	}
	
}

void tnsDrawObjectTree(tnsScene* Scene, tns3DObject* from) {
	tns3DObject* o;

	if (from) {
		for (o = from->ChildObjects.pFirst; o; o = o->Item.pNext) {
			//if (o != Scene->ActiveCamera) {
				tnsPushMatrix();
				tnsApplyObjectMatrix(o);
				tnsDrawThisObject(o, Scene);
				if (o->ChildObjects.pFirst) {
					tnsDrawObjectTree(Scene, o);
				}
				tnsPopMatrix();
			//}
		}
	}
	else {
		if (!Scene) return;
		for (o = Scene->Objects.pFirst; o; o = o->Item.pNext) {
			//if (o != Scene->ActiveCamera) {
				tnsPushMatrix();
				tnsApplyObjectMatrix(o);
				tnsDrawThisObject(o, Scene);
				if (o->ChildObjects.pFirst) {
					tnsDrawObjectTree(Scene, o);
				}
				tnsPopMatrix();
			//}
		}
	}
}

void tnsDrawScene(int W, int H, tnsScene* Scene) {
	if (!Scene->ActiveCamera) return;

	tnsApplyCameraView(W, H, Scene->ActiveCamera);

	tnsDrawObjectTree(Scene, 0);
}
void tnsDrawWorld(int W, int H) {
	tnsWorld* w = &T->World;

	if (!w->ActiveScene) return;

	tnsDrawScene(W, H, w->ActiveScene);
}




//===================================================================[Material]

tnsMaterial* tnsCreateMaterial(tnsScene* s,char* Name) {
	tnsMaterial* m;

	if (!s || !Name || !Name[0]) return 0;

	m = memAquireHyper(sizeof(tnsMaterial));
	strSafeSet(&m->Name, Name);
	lstAppendItem(&s->Materials, m);
	return m;
}
tnsMaterial* tnsFindMaterialByIndex(tnsScene* s, int index) {
	tnsMaterial* m;
	int i = 0;
	if (!s)return 0;
	for (m = s->Materials.pFirst; m; m = m->Item.pNext) {
		if (i == index) {
			m->ID = i;
			return m;
		}
		i++;
	}
	return 0;
}
tnsMaterial* tnsFindMaterial(tnsScene* s, char* name) {
	tnsMaterial* m;
	int i = 0;
	if (!s)return 0;
	for (m = s->Materials.pFirst; m; m = m->Item.pNext) {
		if (strIsTheSame(m->Name->Ptr,name)) {
			return m;
		}
	}
	return 0;
}


//===================================================================[Group]

tnsGroup* tnsCreateGroup(tnsScene* s, char* Name) {
	tnsGroup* g;

	if (!s || !Name || !Name[0]) return 0;

	g = memAquireHyper(sizeof(tnsGroup));
	strSafeSet(&g->Name, Name);
	lstAppendItem(&s->Groups, g);
	return g;
}
tnsGroup* tnsFindGroup(tnsScene* s, char* Name) {
	tnsGroup* g;

	if (!s || !Name || !Name[0]) return 0;

	for (g = s->Groups.pFirst; g;g=g->Item.pNext) {
		if (strIsTheSame(g->Name->Ptr, Name)) return g;
	}

	return 0;
}
void tnsAddToGroup(tns3DObject* o, tnsGroup* g) {
	if (!g || !o) return;

	lstAppendPointer(&g->ObjectLinks, o);
	lstAppendPointer(&o->InGroups, g);
}
void tnsRemoveFromGroup(tns3DObject* o, tnsGroup* g) {
	if (!g || !o) return;

	lstRemovePointer(&g->ObjectLinks, o);
	lstRemovePointer(&o->InGroups, g);
}
int tnsGroupHaveObject(tnsGroup* g, tns3DObject* o) {
	nListItemPointer* lip;
	for (lip = g->ObjectLinks.pFirst; lip; lip = lip->pNext) {
		if (lip->p == o) return 1;
	}
	return 0;
}
int tnsObjectInGroup(tns3DObject* o, tnsGroup* g) {
	nListItemPointer* lip;
	for (lip = o->InGroups.pFirst; lip; lip = lip->pNext) {
		if (lip->p == g) return 1;
	}
	return 0;
}
void tnsCleanObjectFinishMarksRecursive(nListHandle* Objects) {
	tns3DObject* o;
	for (o = Objects->pFirst; o; o = o->Item.pNext) {
		o->MaterialRenderingDone = 0;
		tnsCleanObjectFinishMarksRecursive(&o->ChildObjects);
	}
}
void tnsCleanObjectFinishMarks(tnsScene* s) {
	tns3DObject* o;
	if (!s) return;
	for (o = s->Objects.pFirst; o; o = o->Item.pNext) {
		o->MaterialRenderingDone = 0;
		o->LineRenderingDone = 0;
		tnsCleanObjectFinishMarksRecursive(&o->ChildObjects);
	}
}


//=====================================================================[ACTUATORS]

int ACTCHK_tnsAlignView(nPropPack* This, nStringSplitor* ss) {
	tnsScene* s;
	n3DViewUiExtra* ex;
	tnsCamera* c, *vc;
	if (This) {
		ex = This->EndInstance;
		s = ex->CurrentScene;
		if (!s) return 0;
	}else return 0;

	c = s->ActiveCamera;
	vc = ex->ViewingCamera;

	if (!c||!vc) return 0;
}
int ACTINV_tnsAlignCameraToView(nActuatorIntern* a, nEvent* e) {
	tnsScene* s;
	n3DViewUiExtra* ex;
	tnsCamera* c, *vc;

	if (a->This) {
		ex = a->This->EndInstance;
		s = ex->CurrentScene;
		if (!s) return 0;
	}
	else return NUL_FINISHED;

	c = s->ActiveCamera;
	vc = ex->ViewingCamera;

	if (!c) return NUL_FINISHED;

	tObjApplySameGlobalTransform(c, vc);
	c->FOV = vc->FOV;
	c->ZMax = vc->ZMax;
	c->ZMin = vc->ZMin;


	//nulNotifySubPropUsers()

	return NUL_FINISHED;
}
int ACTINV_tnsAlignViewToCamera(nActuatorIntern* a, nEvent* e) {
	tnsScene* s;
	n3DViewUiExtra* ex;
	tnsCamera* c, *vc;

	if (a->This) {
		ex = a->This->EndInstance;
		s = ex->CurrentScene;
		if (!s) return 0;
	}
	else return NUL_FINISHED;

	c = s->ActiveCamera;
	vc = ex->ViewingCamera;

	if (!c) return NUL_FINISHED;

	tObjApplySameGlobalTransform(vc, c);
	vc->FOV = c->FOV;
	vc->ZMax = c->ZMax;
	vc->ZMin = c->ZMin;

	//nulNotifySubPropUsers()

	return NUL_FINISHED;
}

int ACTCHK_IsScene(nPropPack* This, nStringSplitor* ss) {
	if (This && This->LastPs->p->SubProp == _TNS_PROP_SCENE) return 1;
	return 0;
}
int ACTINV_ActiveScene(nActuatorIntern* a, nEvent*e) {
	tnsScene* s = a->This ? a->This->LastPs->UseInstance : 0;
	if (s)
		T->World.ActiveScene = s;
	else {
		;//Find Scene Based On Name
	}
	nulNotifyUsers("tns.world.scenes");
	return NUL_FINISHED;
}
int ACTINV_DestroyScene(nActuatorIntern* a, nEvent*e) {
	tnsScene* s = a->This ? a->This->LastPs->UseInstance : 0;

	if (!s) {
		;//Find Scene Based On Name
	}

	tnsDestroyScene(s);
	nulNotifyUsers("tns.world.scenes");

	return NUL_FINISHED;
}
int ACTINV_UnlinkScene(nActuatorIntern* a, nEvent*e) {
	tnsScene* s = a->This ? a->This->LastPs->UseInstance : 0;

	if (!s) {
		;//Find Scene Based On Name
	}

	T->World.ActiveScene = 0;
	lstRemoveItem(&T->World.Scenes, s);
	nulThrowToTrashBin(s, "tns_scene");
	nulNotifyUsers("tns.world.scenes");

	return NUL_FINISHED;
}
int ACTINV_NewScene(nActuatorIntern* a, nEvent*e) {
	tnsScene* s = memAquireHyper(sizeof(tnsScene));

	tnsCreateScene("New Scene");

	nulNotifyUsers("tns.world.scenes");

	return NUL_FINISHED;
}
int ACTCHK_Is3DViewExtra(nPropPack* This, nStringSplitor* ss);

int ACTINV_SetActiveCamera(nActuatorIntern* a, nEvent*e) {
	tnsScene* s = T->World.ActiveScene;
	tnsCamera* c = a->This->LastPs->UseInstance;

	if (!c) return NUL_FINISHED;

	tnsSetActiveCamera(s, c);

	nulNotifyUsers("tns.world.scenes.active_camera");

	return NUL_FINISHED;
}


int ACTINV_SaveRenderBufferConfig(nActuatorIntern* a, nEvent*e) {
	nulInvoke(a, "NUL_file_dialog", e, 0, 0, 0);

	return NUL_RUNNING;
}
int ACTMOD_SaveRenderBufferConfig(nActuatorIntern* a, nEvent*e) {
	if (a->ConfirmData) {
		if (a->ConfirmData->StrData) {
			nUDF* udf = nulPrepareUDF(a->ConfirmData->StrData);
			nulWriteProp(udf, "tns.render_buffer_list");
			nulPackUDF(udf, 1, 0);
			return NUL_FINISHED_PASS;
		}
		return NUL_FINISHED_PASS;
	}
	return NUL_RUNNING;
}
int ACTINV_LoadRenderBufferConfig(nActuatorIntern* a, nEvent*e) {
	nulInvoke(a, "NUL_file_dialog", e, 0, 0, 0);

	return NUL_RUNNING;
}
int ACTMOD_LoadRenderBufferConfig(nActuatorIntern* a, nEvent*e) {
	if (a->ConfirmData) {
		if (a->ConfirmData->StrData) {
			nUDF* udf = nulOpenUDF(a->ConfirmData->StrData);
			nulExtractUDF(udf, NUL_UDF_APPEND, 0);
			nulCloseUDF(udf);
			return NUL_FINISHED_PASS;
		}
		return NUL_FINISHED_PASS;
	}
	return NUL_RUNNING;
}



void tnsRegisterAllActuators() {
	nulCreateActuatorType("TNS_activate_scene", "Active Scene", "Set A Scene As The Active Scene",
		ACTCHK_IsScene, 0, 0, ACTINV_ActiveScene, 0, ICON_HAND_POINTING_RIGHT, 0);
	nulCreateActuatorType("TNS_destroy_scene", "Destroy Scene", "Destroy A Scene",
		ACTCHK_IsScene, 0, 0, ACTINV_DestroyScene, 0, ICON_CLOSE_CIRCLE_OUTLINE, 0);
	nulCreateActuatorType("TNS_unlink_scene", "Unlink Scene", "Unlink A Scene And Throw It To Trash Bin",
		ACTCHK_IsScene, 0, 0, ACTINV_UnlinkScene, 0, ICON_LINK_VARIANT_OFF, 0);
	nulCreateActuatorType("TNS_new_scene", "New Scene", "Create A New Scene In The World",
		0, 0, 0, ACTINV_NewScene, 0, ICON_CREATION, 0);

	nulCreateActuatorType("TNS_set_active_camera", "Set Active Camera", "Set active camera as selected one",
		0, 0, 0, ACTINV_SetActiveCamera, 0, ICON_THUMB_DOWN, 0);

	nulCreateActuatorType("TNS_align_camera_to_view", "Align Camera To View", "Align Camera To 3D Ui View Direction",
		ACTCHK_tnsAlignView, 0, 0, ACTINV_tnsAlignCameraToView, 0,ICON_FORMAT_ALIGN_RIGHT, 0);
	nulCreateActuatorType("TNS_align_view_to_camera", "Align View To Camera", "Align 3D View Direction To Active Camera",
		ACTCHK_tnsAlignView, 0, 0, ACTINV_tnsAlignViewToCamera, 0, ICON_FORMAT_ALIGN_RIGHT, 0);

	nulCreateActuatorType("TNS_save_render_buffer_config", "Save Render Buffer Config", "Save render buffer config under this workspace",
		0, 0, 0, ACTINV_SaveRenderBufferConfig, ACTMOD_SaveRenderBufferConfig, ICON_FLOPPY, 0);
	nulCreateActuatorType("TNS_load_render_buffer_config", "Load Render Buffer Config", "Load render buffer config into this workspace",
		0, 0, 0, ACTINV_LoadRenderBufferConfig, ACTMOD_LoadRenderBufferConfig, ICON_PLUS_BOX, 0);

}



//==================================================================[Util]


#define MAX3(a,b,c)\
a>b?(a>c?a:c):(b>c?b:c)

#define MIN3(a,b,c)\
a<b?(a<c?a:c):(b<c?b:c)

extern NUL MAIN;

void tnsSingleLinearToLog(real* a, real gamma) {
	if (*a < 0) *a = 0;
	*a = powf(*a, gamma);
}
void tnsSingleLogToLinear(real* a, real gamma) {
	if (*a < 0) *a = 0;
	*a = powf(*a, 1.0f / gamma);
}
void tnsLinearToLog(real* rgb,real gamma) {
	if (rgb[0] < 0) rgb[0] = 0;
	if (rgb[1] < 0) rgb[1] = 0;
	if (rgb[2] < 0) rgb[2] = 0;
	rgb[0] = powf(rgb[0], gamma);
	rgb[1] = powf(rgb[1], gamma);
	rgb[2] = powf(rgb[2], gamma);
}
void tnsLogToLinear(real* rgb, real gamma) {
	if (rgb[0] < 0) rgb[0] = 0;
	if (rgb[1] < 0) rgb[1] = 0;
	if (rgb[2] < 0) rgb[2] = 0;
	rgb[0] = powf(rgb[0], 1.0f / gamma);
	rgb[1] = powf(rgb[1], 1.0f / gamma);
	rgb[2] = powf(rgb[2], 1.0f / gamma);
}
void tnsRgbToLuminance(real* rgb) {
	real l;
	if (rgb[0] < 0) rgb[0] = 0;
	if (rgb[1] < 0) rgb[1] = 0;
	if (rgb[2] < 0) rgb[2] = 0;
	l = rgb[0] * 0.299 + rgb[1] * 0.587 + rgb[2] * 0.114;
	rgb[0] = rgb[1] = rgb[2] = l;
}
void tnsHSVtoRGB(real* hsv,real* rgb) {
	int i;
	real f, p, q, t;
	real h = (hsv[0]>1?(hsv[0]-1.0)/6:hsv[0]*6), s = hsv[1], v = hsv[2];

	if (MAIN.ColorAccessCorrectGamma)
		s = hsv[1] = powf(hsv[1], MAIN.Gamma);

	if (hsv[1] == 0) {
		rgb[0] = rgb[1] = rgb[2] = hsv[2];
		return;
	}
	i = (int)floor(h);
	f = h - i;
	p = v*(1 - s);
	q = v*(1 - s*f);
	t = v*(1 - s*(1 - f));
	switch (i) {
	case 0:
		rgb[0] = v;
		rgb[1] = t;
		rgb[2] = p;
		break;
	case 1:
		rgb[0] = q;
		rgb[1] = v;
		rgb[2] = p;
		break;
	case 2:
		rgb[0] = p;
		rgb[1] = v;
		rgb[2] = t;
		break;
	case 3:
		rgb[0] = p;
		rgb[1] = q;
		rgb[2] = v;
		break;
	case 4:
		rgb[0] = t;
		rgb[1] = p;
		rgb[2] = v;
		break;
	default:
		rgb[0] = v;
		rgb[1] = p;
		rgb[2] = q;
		break;
	}
}
void tnsRGBtoHSV(real* rgb, real* hsv) {
	real min, max, delta;
	real r = rgb[0], g = rgb[1], b = rgb[2];
	min = MIN3(r, g, b);
	max = MAX3(r, g, b);
	hsv[2] = max;
	delta = max - min;
	if (max != 0) hsv[1] = delta / max;
	else {
		hsv[1] = .00005;
		hsv[0] = 0;
		return hsv;
	}
	if (delta == 0) {
		hsv[1] = .00005;
		hsv[0] = 0;
		return hsv;
	}
	if (r == max) hsv[0] = (g - b) / delta;
	else if (g == max) hsv[0] = 2 + (b - r) / delta;
	else hsv[0] = 4 + (r - g) / delta;
	hsv[0] /= 6;
	if (hsv[0]<0) hsv[0] +=1;
	if (hsv[0] >= 1) hsv[0] -=1;

	if (MAIN.ColorAccessCorrectGamma)
		hsv[1] = powf(hsv[1], 1.0f/MAIN.Gamma);
}


void tnsClearAll() {
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
}
void tnsClearColorv(real* rgba) {
	glClearColor(rgba[0], rgba[1], rgba[2], rgba[3]);
}
void tnsClearColor(real r,real g,real b,real a) {
	glClearColor(r,g,b,a);
}

//Assuming *result is long enough
void tnsMakeIndexUShort(unsigned short* result, short num, ...){
	short i;
	va_list list;

	va_start(list,num);

	for (i = 0; i < num; i++){
		result[i] = va_arg(list, unsigned short);
	}

	va_end(list);
}
void tnsMakeTriangle(real* arr, real x1, real y1, real x2, real y2, real x3, real y3){
	arr[0] = x1;
	arr[1] = y1;
	arr[2] = x2;
	arr[3] = y2;
	arr[4] = x3;
	arr[5] = y3;
}
void tnsMakeQuad2d(real* arr, real x1, real y1, real x2, real y2, real x3, real y3,real x4,real y4){
	arr[0] = x1;
	arr[1] = y1;
	arr[2] = x2;
	arr[3] = y2;
	arr[4] = x3;
	arr[5] = y3;
	arr[6] = x4;
	arr[7] = y4;
}
void tnsMakeQuad3d(real* arr, real x1, real y1, real z1, real x2, real y2, real z2, real x3, real y3, real z3, real x4, real y4, real z4){
	arr[0] = x1;
	arr[1] = y1;
	arr[2] = z1;
	arr[3] = x2;
	arr[4] = y2;
	arr[5] = z2;
	arr[6] = x3;
	arr[7] = y3;
	arr[8] = z3;
	arr[9] = x4;
	arr[10] = y4;
	arr[11] = z4;
}
void tnsMakeQuad4d(real* arr, real x1, real y1, real z1, real w1, real x2, real y2, real z2, real w2, real x3, real y3, real z3, real w3, real x4, real y4, real z4, real w4){
	arr[0] = x1;
	arr[1] = y1;
	arr[2] = z1;
	arr[3] = w1;
	arr[4] = x2;
	arr[5] = y2;
	arr[6] = z2;
	arr[7] = w2;
	arr[8] = x3;
	arr[9] = y3;
	arr[10] = z3;
	arr[11] = w3;
	arr[12] = x4;
	arr[13] = y4;
	arr[14] = z4;
	arr[15] = w4;
}
void tnsMakeCircle2d(real* arr,int slices, real ctrX, real ctrY, real r) {
	int i;
	real radstep = 2 * TNS_PI / (real)slices;
	real rd = 0;
	real x, y;
	for (i = 0; i < slices;i++) {
		x = ctrX + cos(rd) * r;
		y = ctrY + sin(rd) * r;
		arr[2 * i] = x;
		arr[2 * i + 1] = y;
		rd += radstep;
	}
}
void tnsMakeArc2d(real* arr, int slices, real ctrX, real ctrY, real r,real rad_begin,real rad_end) {
	int i;
	real radstep = (rad_end-rad_begin) / (real)slices;
	real rd = rad_begin;
	real x, y;
	for (i = 0; i <= slices; i++) {
		x = ctrX + cos(rd) * r;
		y = ctrY + sin(rd) * r;
		arr[2 * i] = x;
		arr[2 * i + 1] = y;
		rd += radstep;
	}
}

void tnsMakeBridgedIndex(unsigned short* result, short num,int revert, short begin) {
	short i;

	if (!revert) {
		for (i = 0; i < num; i++) {
			result[i * 2] = begin + i;
			result[i * 2 + 1] = begin + i + num;
		}
	}else{
		for (i = 0; i < num; i++) {
			result[i * 2] = begin + i;
			result[i * 2 + 1] = begin - i + num * 2;
		}
	}
}

void tnsMakeLinerGradient3d(real* arr,int num_points, real r0, real g0, real b0, real r1, real g1, real b1){
	int i = 0;
	real sr = (r1 - r0) / (real)num_points
		, sg = (g1 - g0) / (real)num_points
		, sb = (b1 - b0) / (real)num_points;
	for (i; i < num_points; i++) {
		int a = 3 * i;
		arr[a] = r0 + sr*i;
		arr[a + 1] = g0 + sg*i;
		arr[a + 2] = b0 + sb*i;
	}
}
void tnsMakeLinerGradient4d(real* arr, int num_points, real r0, real g0, real b0, real a0, real r1, real g1, real b1, real a1) {
	int i = 0;
	real sr = (r1 - r0) / (real)num_points
		, sg = (g1 - g0) / (real)num_points
		, sb = (b1 - b0) / (real)num_points
		, sa = (a1 - a0) / (real)num_points;
	for (i; i < num_points; i++) {
		int a = 4 * i;
		arr[a] = r0 + sr*i;
		arr[a + 1] = g0 + sg*i;
		arr[a + 2] = b0 + sb*i;
		arr[a + 3] = a0 + sa*i;
	}
}
void tnsMakeLinerGradient3dv(real* arr, int num_points, real* rgb0,real* rgb1) {
	int i = 0;
	real sr = (rgb1[0] - rgb0[0]) / (real)num_points
		, sg = (rgb1[1] - rgb0[1]) / (real)num_points
		, sb = (rgb1[2] - rgb0[2]) / (real)num_points;
	for (i; i < num_points; i++) {
		int a = 3 * i;
		arr[a] = rgb0[0] + sr*i;
		arr[a + 1] = rgb0[1] + sg*i;
		arr[a + 2] = rgb0[2] + sb*i;
	}
}
void tnsMakeLinerGradient4dv(real* arr, int num_points, real* rgb0, real* rgb1) {
	int i = 0;
	real sr = (rgb1[0] - rgb0[0]) / (real)(num_points - 1)
		, sg = (rgb1[1] - rgb0[1]) / (real)(num_points - 1)
		, sb = (rgb1[2] - rgb0[2]) / (real)(num_points - 1)
		, sa = (rgb1[3] - rgb0[3]) / (real)(num_points - 1);
	for (i; i < num_points; i++) {
		int a = 4 * i;
		arr[a] = rgb0[0] + sr*i;
		arr[a + 1] = rgb0[1] + sg*i;
		arr[a + 2] = rgb0[2] + sb*i;
		arr[a + 3] = rgb0[3] + sa*i;
	}
}

void tnsMakeFoucsSquare(int L, int R, int U, int B,int W) {
	int w = W;
	//int al = Len + w;
	real v[16];

	//v[0] = L - w; v[1] = U - l;
	//v[2] = L - w; v[3] = U + w;
	//v[4] = L + l; v[4] = U + w;
	//v[6] = L + l; v[7] = U - w;
	//v[8] = L + w; v[9] = U - w;
	//v[10] = L + w; v[11] = U - l;

	//v[12] = R - l; v[13] = U + w;
	//v[14] = R + w; v[15] = U + w;
	//v[16] = R + w; v[17] = U - l;
	//v[18] = R - w; v[19] = U - l;
	//v[20] = R - w; v[21] = U - w;
	//v[22] = R - l; v[23] = U - w;

	tnsMakeQuad2d(v, L, U, R, U, R, B, L, B);
	tnsVertexArray2d(v, 4);
	tnsPackAs(GL_LINE_LOOP);
	tnsMakeQuad2d(&v[8], L + W, U - W, R - W, U - W, R - W, B + W, L + W, B + W);
	tnsVertexArray2d(&v[8], 4);
	tnsPackAs(GL_LINE_LOOP);
}

void tnsDrawFloor(int Size, int Span, int* ShowAxis) {
	int i = 0;
	int Lim = Span * 2 + 1;
	int Dist = Size * Span;

	for (i; i < Lim; i++) {
		if (i == Span && ShowAxis[0])continue;
		tnsVertex3d(-Dist, i*Size - Dist, 0);
		tnsVertex3d(Dist, i*Size - Dist, 0);
	}

	for (i = 0; i < Lim; i++) {
		if (i == Span && ShowAxis[1])continue;
		tnsVertex3d(i*Size -Dist, -Dist, 0);
		tnsVertex3d(i*Size - Dist, Dist, 0);
	}

	tnsPackAs(GL_LINES);

	if (ShowAxis[0]) {
		tnsColor4d(1, 0, 0, 1);
		tnsVertex3d(-Dist, 0, 0);
		tnsVertex3d(Dist, 0, 0);
		tnsPackAs(GL_LINES);
	}
	if (ShowAxis[1]) {
		tnsColor4d(0, 1, 0, 1);
		tnsVertex3d(0, -Dist, 0);
		tnsVertex3d(0, Dist, 0);
		tnsPackAs(GL_LINES);
	}
	if (ShowAxis[2]) {
		tnsColor4d(0, 0, 1, 1);
		tnsVertex3d(0, 0, -Dist);
		tnsVertex3d(0, 0, Dist);
		tnsPackAs(GL_LINES);
	}

}
void tnsDrawFloorWithTheme(int Size, int Span, int* ShowAxis,nThemeState* ts) {
	int i = 0;
	int Lim = Span * 2 + 1;
	int Dist = Size * Span;
	if (!ts || !ts->BackColor || !ts->BorderColor || !ts->TextColor || !ts->SecondColor) return;

	tnsColor4dv(ts->TextColor->RGBA);

	for (i; i < Lim; i++) {
		if (i == Span && ShowAxis[0])continue;
		tnsVertex3d(-Dist, i*Size - Dist, 0);
		tnsVertex3d(Dist, i*Size - Dist, 0);
	}

	for (i = 0; i < Lim; i++) {
		if (i == Span && ShowAxis[1])continue;
		tnsVertex3d(i*Size - Dist, -Dist, 0);
		tnsVertex3d(i*Size - Dist, Dist, 0);
	}

	tnsPackAs(GL_LINES);

	if (ShowAxis[0]) {
		tnsColor4dv(ts->BorderColor->RGBA);
		tnsVertex3d(-Dist, 0, 0);
		tnsVertex3d(Dist, 0, 0);
		tnsPackAs(GL_LINES);
	}
	if (ShowAxis[1]) {
		tnsColor4dv(ts->BackColor->RGBA);
		tnsVertex3d(0, -Dist, 0);
		tnsVertex3d(0, Dist, 0);
		tnsPackAs(GL_LINES);
	}
	if (ShowAxis[2]) {
		tnsColor4dv(ts->SecondColor->RGBA);
		tnsVertex3d(0, 0, -Dist);
		tnsVertex3d(0, 0, Dist);
		tnsPackAs(GL_LINES);
	}

}



void tnsGetGaussianKernel(tnsFilterKernel* fk, int size, real sigma){
	const double PI = TNS_PI;
	int center = size / 2;

	if (fk->Kernel) FreeMem(fk->Kernel);
	fk->Kernel = CreateNewBuffer(real, size*size);

	for (int i = 0; i<size; i++)
	{
		for (int j = 0; j<size; j++)
		{
			fk->Kernel[i*size+j] = (1 / (2 * PI*sigma*sigma))*exp(-((i - center)*(i - center) + (j - center)*(j - center)) / (2 * sigma*sigma));
		}
	}

	fk->Size = size;
}



