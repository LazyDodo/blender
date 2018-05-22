#pragma once

#include "NUL4.h"
#include "tinycthread.h"
#include <float.h>

#define TNS_PI 3.1415926535897932384626433832795
#define deg(r) r/TNS_PI*180.0
#define rad(d) d*TNS_PI/180.0


//typedef real tnsMatrix33d[9];

typedef float tnsMatrix44f[16];

typedef real tnsMatrix44d[16];
typedef real tnsVector2d[2];
typedef real tnsVector3d[3];
typedef real tnsVector4d[4];

typedef int tnsVector2i[2];

//typedef double tnsMatrix33d[9];
//typedef double tnsMatrix44d[16];
//
typedef struct _tnsMatrixStackItem tnsMatrixStackItem;
struct _tnsMatrixStackItem{
	tnsMatrix44d view;
	tnsMatrix44d model;
	tnsMatrix44d projection;
};
typedef struct _tnsMatrixStack tnsMatrixStack;
struct _tnsMatrixStack{
	tnsMatrixStackItem* level;
	int                 current_level;
	int                 max_level;
};

typedef struct _tnsShader tnsShader;
struct _tnsShader{
	nListItem          Item;
	int                vtShaderID;
	int                fgShaderID;
	int                glProgramID;

	int                CustomID;

	int                modelIndex;
	int                projectionIndex;
	int                viewIndex;
	int                projectionInverseIndex;

	int                vertexIndex;
	int                normalIndex;
	int                colorIndex;
	int                uvIndex;
	int                texture0Index;
	int                texture1Index;
	int                texture2Index;
	int                texture3Index;
	int                texture4Index;

	int                uniform0Index;
	int                uniform1Index;
	int                uniform2Index;
	int                uniform3Index;
	int                uniform4Index;
	int                uniform5Index;
	int                uniform6Index;
	int                uniform7Index;
};
typedef struct _tnsTexture tnsTexture;
typedef struct _tnsCommand tnsCommand;
struct _tnsCommand{
	GLenum Mode;
	short  Dimentions;// 2 or 3
	char   UseVert;
	char   UseColor;
	char   UseNormal;
	char   UseTexCoord;
	char   UseIndex;
	GLenum PolyMode;//0-solid 1-wire
	GLenum Shade;//0-falt 1-smooth
	GLfloat UniformColor[4];

	short NumVert;
	short NumIndex;
	short VertBegin;
	short VertEnd;//'END'is the next one after the last one.
	short ColorBegin;
	short ColorEnd;
	short TexCoordBegin;
	short TexCoordEnd;
    short NormalBegin;
	short NormalEnd;
	GLushort IndexBegin;
	GLushort IndexEnd;

	tnsShader*  ReplaceShader;
	tnsTexture* Texture0;
	tnsTexture* Texture1;
	tnsTexture* Texture2;
};
typedef struct _tnsMain tnsMain;
typedef struct _tnsScene tnsScene;
typedef struct _tnsWorld tnsWorld;
struct _tnsWorld{
	nHyperItem  Item;
	nListHandle Scenes;
	tnsScene*   ActiveScene;

	//Basic Properties;
	u16bit      TimeYear;
	u8bit       TimeMonth;
	u8bit       TimeDay;

	nListHandle BezierCurves;
	nListHandle  Meshes;
};

NEED_STRUCTURE(tnsLoopItem);
NEED_STRUCTURE(tnsRenderLine);
NEED_STRUCTURE(tnsRenderBuffer);

STRUCTURE(tnsFilterKernel) {
	int Size;
	real* Kernel;
};

STRUCTURE(tnsTriangulateNode) {
	char  Picked;
	real Angle;//rad
	tnsLoopItem* LoopItem;
	tnsRenderLine* FowardRL;
	tnsRenderLine* BackwardRL;
};
STRUCTURE(tnsTriangulateEdgeNode) {
	tnsRenderLine* RL;

};

struct _tnsMain{
	nHyperItem     Item;

	tnsMatrixStack stack;

	nListHandle    Shaders;
	int            NextShaderIndex;
	tnsShader*     CurrentShader;
	tnsShader*     BindedShader;

	//int            MatrixMode;
	int            IsOffscreen;

	char*          GLVersionStr;
	char*          GLVendorStr;
	char*          GLRendererStr;
	char*          GLSLVersionStr;


	tnsShader*     uiShader;
	tnsShader*     stringShader;
	tnsShader*     TextureShader;
	tnsShader*     RectangleTextureShader;
	tnsShader*     TextureMultiplyShader;
	tnsShader*     MSTextureShader;
	tnsShader*     MSATextureShader;
	tnsShader*     TEST_MatcapShader;
	tnsShader*     TransparentGridShader;
	tnsShader*     SobelColorShader;
	tnsShader*     SobelShader;
	tnsShader*     ExtraBuffersShader;

	tnsShader*     AtlasTransformShader;
	tnsShader*     AtlasPreviewShader;

	tnsShader*     ImagePeelShader;
	tnsShader*     LineConnectionShader;

	nListHandle    Textures;
	tnsTexture*    PreviewTexture;

	tnsFilterKernel EdgeGaussFilter;

	GLenum         GlTextureSets;

	GLenum         StateShadeMode;
	GLenum         StatePolyMode;
	tnsTexture*    StateTexture0;
	tnsTexture*    StateTexture1;
	tnsTexture*    StateTexture2;
	tnsTexture*    Texture0Enabled;
	tnsTexture*    Texture1Enabled;
	tnsTexture*    Texture2Enabled;
	tnsTexture*    Texture3Enabled;
	tnsTexture*    Texture4Enabled;

	GLfloat        StateColor[4];
	
	GLfloat        Vert[8192];
	GLuint         VertBufObject;
	short          NextVert;
	
	GLfloat        Color[8192];
	GLuint         ColorBufObject;
	short          NextColor;

	GLfloat        Normal[8192];
	GLuint         NormalBufObject;
	short          NextNormal;

	GLfloat        TexCoord[8192];
	GLuint         TexCoordBufObject;
	short          NextTexCoord;


	GLushort       Index[1024];
	GLuint         IndexBufObject;
	GLushort       NextIndex;

	tnsTriangulateNode SharedTN[4096];
	u32bit             SharedTEBuf[4096];
	tnsRenderLine*     SharedRLBuf[4096];

	tnsCommand     DrawingCommand[128];
	tnsCommand*    UsingCommand;

	tnsWorld       World;

	nListHandle      RenderBuffers;
	tnsRenderBuffer* ActiveRenderBuffer;
};

typedef struct _tnsLineStripPoint tnsLineStripPoint;
struct _tnsLineStripPoint {
	nListItem  Item;
	tnsVector3d P;
};

typedef struct _tnsLineStrip tnsLineStrip;
struct _tnsLineStrip {
	nListItem Item;
	nListHandle Points;
	int PointCount;
};

#define TNS_SNAKE_STRENGTH_LIMIT 30
#define TNS_SNAKE_STEP_LENGTH 5
#define TNS_SNAKE_STRENGTH_MARCHING_LIMIT 30

#define TNS_SNAKE_EDGE_WIDTH 2
#define TNS_SNAKE_FILTER_SIZE 6
#define TNS_SNAKE_ANGLE_DEVIATE 0.5
#define TNS_SNAKE_ANGLE_PRICISION 1
#define TNS_SNAKE_STEP1 3
#define TNS_SNAKE_STEP2 5


typedef struct _tnsTextureSample tnsTextureSample;
struct _tnsTextureSample {
	nListItem Item;
	u8bit     Sample;
	int       X,Y;
	real      Z;
};

typedef struct _tnsTexture tnsTexture;
struct _tnsTexture{
	nListItem Item;

	int    IsRenderBuffer;

	GLuint GLTexHandle;
	GLuint GLTexBitsType;//like GL_RGBA
	GLuint DataFormat;
	GLuint DataType;
	GLuint GLTexType;    //like GL_TEXTURE_2D
	int    Width;
	int    Height;

	void*  DrawData;


	int    RBWidth, RBHeight;
	int    ElemSize;
	void*  TextureReadBack;
	void** SamplePtr;
	nListHandle PendingSamples;
	nListHandle ErasedSamples;

	nListHandle LineStrips;
};

typedef struct _tnsOffscreen tnsOffscreen;
struct _tnsOffscreen{
	nListItem  Item;

	tnsTexture* pColorTextures[16];
	tnsTexture* pDepthTexture;
	tnsTexture* pStencilTexture;

	GLuint      FboHandle;

	int         UseSecondary;
};

#define TNS_PROJECTION_MATRIX 1
#define TNS_MODEL_MATRIX  2
#define TNS_VIEW_MATRIX   3
#define TNS_TEXTURE_MATRIX    4

//==========================================================================[FT]
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include "ft2build.h"
#include "freetype/ftglyph.h"
#include "freetype/ftoutln.h"
#include "freetype/fttrigon.h"


extern char TNS_VERTEX_SHADER_SRC_UI_DEFAULT_130[];
extern char TNS_FRAGMENT_SHADER_SRC_UI_DEFAULT_130[];
extern char TNS_VERTEX_SHADER_SRC_UI_TEXT_130[];
extern char TNS_FRAGMENT_SHADER_SRC_UI_TEXT_130[];
extern char TNS_FRAGMENT_SHADER_SRC_SINGLE_TEXTURE_130[];


NEED_STRUCTURE(tnsTexture);

typedef struct _tnsFontSingleCharacter tnsFontSingleCharacter;
struct _tnsFontSingleCharacter{
	tnsTexture* Tex;
	int        Generated;
	wchar_t    charID;
	int        width;
	int        height;
	int        advx;
	int        advy;
	int        deltax;
	int        deltay;

	int        bufferx;
	int        buffery;
};

typedef struct _tnsFont tnsFont;
struct _tnsFont{
	nListItem Item;

	char*    fontName;
	char*    IconFontName;

	tnsFontSingleCharacter characters[65536];
	FT_Library             ftlib;
	FT_Face                ftface;

	tnsFontSingleCharacter icons[4096];
	FT_Library             iconftlib;
	FT_Face                iconftface;

	unsigned int           height;

	tnsTexture             TexBuffer;
	int                    CurrentX;
	int                    CurrentY;

};

typedef struct _tnsFontManager tnsFontManager;
struct _tnsFontManager{
	nListHandle Fonts;
	tnsFont*    UsingFont;

	tnsFont*    VectorsGrapghs;

	unsigned    LastDlst;

};

typedef struct _tnsFontBoundBox tnsFontBoundBox;
struct _tnsFontBoundBox{
	int x, y, w, h;/* Upper Left and width,height */
};

#define TNS_FONT_BUFFER_W 2048
#define TNS_FONT_BUFFER_H 2048

#define TNS_FONT_ALIGN_LEFT   1
#define TNS_FONT_ALIGN_CENTER 2
#define TNS_FONT_ALIGN_RIGHT  4
#define TNS_FONT_ALIGN_LEFT_PROTECT 8
#define TNS_FONT_ALIGN_CENTER_FIRST (2|8)


#define CSS_FLOW_DIRECTION_LR 1
#define CSS_FLOW_DIRECTION_TD 0

//typedef struct _tnsCssBorder tnsCssBorder;
//struct _tnsCssBorder{
//	real Px;//thickness
//	real Color[4];//rgba
//	int   Style;//Like________ or ....... or __.__.__. or ======== or ~~~~~~~
//	int   Padding;
//	int   Margin;
//};
//
//typedef struct _tnsCssBackground tnsCssBackground;
//struct _tnsCssBackground{
//	real Color[4];
//	int   ImageID;//NUL SPECIFIC
//	int   Repeat;
//	int   PositionFixed;
//	real PositionRatio;
//};
//
//typedef struct _tnsCssText tnsCssText;
//struct _tnsCssText{
//	real Color[4];
//	int LineSpacing;
//	int LetterSpacing;
//	int Align;//0-LR 1-C 2-L 3-R all support number-ready display
//	int Decoration;
//	int Underline;
//	int Shadow;
//};
//
//typedef struct _tnsCssMisc tnsCssMisc;
//struct _tnsCssMisc{
//	real       Color[4];
//	int         Size;
//	const char* Name;
//	const char* Description;
//};
//
//typedef struct _tnsCssState tnsCssState;
//struct _tnsCssState{
//	char             NameReplace[36];
//	tnsCssBackground Bkg;
//	tnsCssBorder	 BorderLeft;
//	tnsCssBorder	 BorderRight;
//	tnsCssBorder	 BorderTop;
//	tnsCssBorder	 BorderBottom;
//	tnsCssText       Text;
//	tnsCssMisc       Misc1;
//	tnsCssMisc       Misc2;
//	tnsCssMisc       Misc3;
//	tnsCssMisc       Misc4;
//};

#define TNS_STATE_NORMAL           0
#define TNS_STATE_PUSHED           1
#define TNS_STATE_HIGHLIGHT        2
#define TNS_STATE_HIGHLIGHT_SELECT 3
#define TNS_STATE_KEYING           4
#define TNS_STATE_INTERPOLATING    5
#define TNS_STATE_LIMITED          6
//
//typedef struct _tnsCssStyle tnsCssStyle;
//struct _tnsCssStyle{
//	tnsCssState Normal;
//	tnsCssState Pushed;
//	tnsCssState Highlight;
//	tnsCssState HighlightSelect;
//	tnsCssState Keying;
//	tnsCssState Interpolating;
//	tnsCssState Limited;
//};

//=====================================================[3d comp]

typedef struct _tnsWorld tnsWorld;
typedef struct _tns3DObject tns3DObject;

NEED_STRUCTURE(tnsTexture);
NEED_STRUCTURE(tnsMaterial);
typedef struct _tnsScene tnsScene;
struct _tnsScene{
	nHyperItem   Item;

	int         ID;
	nSafeString* Name;

	//List [ROOT]s

	nListHandle  Objects;

	nListHandle  AllObjects;

	nListHandle  Materials;
	tnsMaterial* MaterialPointers[2048];

	nListHandle Groups;

	tns3DObject* ActiveCamera;

	int          AtlasCount;
	tnsTexture*  AtlasPointsL;
	tnsTexture*  AtlasPointsR;
	tnsTexture*  AtlasFaceNormalL;
	tnsTexture*  AtlasFaceNormalR;
};

#define TNS_ROTATION_XYZ_EULER  0
#define TNS_ROTATION_XZY_EULER  1
#define TNS_ROTATION_YXZ_EULER  2
#define TNS_ROTATION_YZX_EULER  3
#define TNS_ROTATION_ZXY_EULER  4
#define TNS_ROTATION_ZYX_EULER  5
#define TNS_ROTATION_QUATERNION 6

#define TNS_OBJECT_PLACEHOLDER  0
#define TNS_OBJECT_CAMERA       1
#define TNS_OBJECT_LAMP         2
#define TNS_OBJECT_MESH         4
#define TNS_OBJECT_BEZIER_CURVE 8

//typedef struct _tns3DObjectLinker tns3DObjectLinker;
//struct _tns3DObjectLinker {
//	nListItem    Item;
//	tns3DObject* o;
//};

NEED_STRUCTURE(tnsRenderTriangle);
NEED_STRUCTURE(tnsRenderVert);

typedef struct _tns3DObject tns3DObject;
struct _tns3DObject{
	nHyperItem    Item;

	int          Type;
	u64bit       ID;
	nSafeString* Name;
	int          DrawMode;

	int          Show;
	int          ShowOnRender;
	int          HideChildren;

	real         GLocation[3];
	real         GRotation[3];
	real         GScale;//Not Basically Used

	real         Location[3];//Loc Rot Scale All Local
	real         Rotation[4];//vector value also used to determing lamp directons.
	real         Scale;//Only allow uniform scale
	u8bit        RotationMode;

	tnsMatrix44d GlobalTransform;
	tnsMatrix44d SelfTransform;
	//tnsMatrix44d CameraTransform;

	tnsScene*    InScene;
	tns3DObject* ParentObject;
	nListHandle  ChildObjects;
	nListHandle  ChildReadTemp;

	tns3DObject* ConstraintObject;

	tnsRenderVert*     RenderVertices;
	tnsRenderTriangle* RenderTriangles;

	nListHandle  InGroups;
	int          LineRenderingDone;//'exclusive' mark
	int          MaterialRenderingDone;//'exclusive' mark

	void*        ExtraRead;
};

STRUCTURE(tnsChildObjectReadTemp) {
	nListItem    Item;
	char         Name[128];
};

#define TNS_PRESPECTIVE_CAMERA 0
#define TNS_ORTHOGRAPHICAL_CAMERA 1
#define TNS_FISHEYE_CAMERA 2

NEED_STRUCTURE(tnsBatch)

#define TNS_MATERIAL_DRAW_MODE_SOLID 0
#define TNS_MATERIAL_DRAW_MODE_WIRE  1

typedef struct _tnsMaterial tnsMaterial;
struct _tnsMaterial{
	nHyperItem Item;

	nSafeString* Name;

	int       ID;

	int       DrawMode;

	real      Color[4];
	real      LinearColor[4];
	real      SpectacularColor[4];
	real      ReflexThreshold;
	real      ReflexStrengh;
	real      ReflexSharpeness;

	real      LowBrightnessColor[4];//Windows are lighten up in the evenings
	u8bit     HaloMode;

	//tnsBatch* PreviewBatch;
	//long      PreviewVCount;
};

#define TNS_CAMERA_PERSPECTIVE 0
#define TNS_CAMERA_ORTHO 1
#define TNS_CAMERA_FISHEYE 2


typedef struct _tnsCamera tnsCamera;
struct _tnsCamera{
	tns3DObject Base;
	
	int   CameraType;
	real  FOV;
	real  ZMin, ZMax;
	real  FocusDistance;

	real  OrthScale;
	tnsVector3d RenderViewDir;
};

typedef struct _tnsLamp tnsLamp;
struct _tnsLamp{
	tns3DObject  Base;

	int          LampType;
	u8bit        UniformDirectional;
	u8bit        Running;

	tnsMaterial* Material;
};

STRUCTURE(tnsGroup) {
	nHyperItem   Item;

	nListHandle  ObjectLinks;

	nSafeString* Name;

	tnsVector4d  Color1;
	tnsVector4d  Color2;

	int          ExcludeFromCalculation;
};

STRUCTURE(tnsVert) {
	nListItem Item;
	real      P[3];//Position
	real      N[3];//Normal
	//real      C[4];//Color
	//real      UV[3];//UV
	char      Selected;
	u32bit    I;//Index
	nListHandle EdgeItems;
	nListHandle FaceItems;
	tnsRenderVert* RV;
};
NEED_STRUCTURE(tnsFace);
NEED_STRUCTURE(tnsRenderLine);
NEED_STRUCTURE(tnsRenderVert);
STRUCTURE(tnsEdge) {
	nListItem Item;
	tnsVert*  VL;
	tnsVert*  VR;
	tnsFace*  FL;
	tnsFace*  FR;
	char      Selected;
	tnsRenderLine* RenderLine;
	tnsRenderLine* RenderLinePositiveSide;//also no intersection
	tnsRenderLine* RenderLineNegativeSide;
	tnsRenderVert* Intersected;
};
STRUCTURE(tnsEdgeItem) {
	nListItem Item;
	tnsEdge*  E;
};
STRUCTURE(tnsLoopItem) {
	nListItem Item;
	tnsEdge*  E;
	tnsVert*  Begin;
};
STRUCTURE(tnsFace) {
	nListItem   Item;
	nListHandle Loop; // MUST IN CCW ORDER!
	tnsVector3d FaceNormal;
	tnsVector3d GNormal;
	tnsVector3d Center;
	short       TriangleCount;
	char        Selected;
	short       MaterialID;
};

STRUCTURE(tnsBranchedCommand) {
	nListItem Item;
	GLuint    EBO;//elem
	u32bit    ElementCount;
	GLenum    DrawAs;
	int       InternalMode;
	int       Dimention;
};
STRUCTURE(tnsBatch) {
	nListItem   Item;
	GLuint      VBO;//vert
	GLuint      NBO;//normal
	GLuint      CBO;//color
	u32bit      NumVert;
	int         Dimention;
	nListHandle Branches;
	GLuint      BeginElementOffset;
};

//STRUCTURE(tnsMesh) {
//	nHyperItem  Item;
//
//	nSafeString Name;
//
//	nListHandle V; u32bit numV;
//	nListHandle E; u32bit numE;
//	nListHandle F; u32bit numF;
//	u32bit      TriangleCount;
//	u32bit      TriangulatedEdgeCount;
//	tnsBatch*   Batch;
//
//};

STRUCTURE(tnsMeshObject) {
	tns3DObject Base;
	nListHandle V; u32bit numV;
	nListHandle E; u32bit numE;
	nListHandle F; u32bit numF;
	u32bit      TriangleCount;
	u32bit      TriangulatedEdgeCount;
	tnsBatch*   Batch;
	tnsBatch*   AtlasTriggerBatch;
	//nListHandle PointerPool;
	//nListHandle LoopPool;
};

NEED_STRUCTURE(tnsBezierCurve);
STRUCTURE(tnsBezierObject) {
	tns3DObject     Base;
	int             Subdiv;
	tnsBezierCurve* Curve;
	tnsBatch*       Batch;
};


//==================================================================[curve/surface types]


STRUCTURE(tnsBezierHandle) { //in TNS bezier curves are all cubic
	nListItem Item;
	u8bit     LMode;
	u8bit     RMode;
	real      P[3];
	real      LP[3];
	real      RP[3];
	real      Thickness;
	real      Tilt;
};
STRUCTURE(tnsBezierCurve) { //in TNS bezier curves are all cubic
	nHyperItem   Item;
	int          Segments;
	nListHandle  Handles;
};


//=======================================================================[software render]


STRUCTURE(tnsRenderVert) {
	nListItem   Item;
	tnsVector4d GLocation;
	tnsVector4d FrameBufferCoord;
	tnsVector2i FrameBufferCoordi;
	tnsVert*    V;           //Used As R When Intersecting
	tnsRenderLine*     IntersectingLine;
	tnsRenderVert*     IntersectintLine2;
	tnsRenderTriangle* IntersectWith;     //   Positive 1         Negative 0
	//tnsRenderTriangle* IntersectingOnFace;//     <|                  |>
	char        Positive;  //                 L---->|----->R	 L---->|----->R
	char        EdgeUsed;  //                      <|		           |>
};

#define TNS_CULL_DISCARD 2
#define TNS_CULL_USED    1

STRUCTURE(tnsRenderTriangle) {
	nListItem      Item;
	tnsRenderVert* V[3];
	tnsRenderLine* RL[3];
	tnsVector3d    GN;
	tnsVector3d    GC;
	tnsFace*       F;
	nListHandle    IntersectingVerts;
	char           CullStatus;
	tnsRenderLine* Testing;	//Should Be tRT** Testing[NumOfThreads]
};

STRUCTURE(tnsRenderTriangleThread) {
	tnsRenderTriangle Base;
	tnsRenderLine*    Testing[128]; //max thread support;
};

STRUCTURE(tnsRenderElementLinkNode) {
	nListItem Item;
	void*     Pointer;
	int       ElementCount;
	void*     ObjectRef;
	char      Additional;
};

STRUCTURE(tnsRenderLineSegment) {
	nListItem Item;
	//real     Begin, End;  // 0->At[L] 1->At[R]
	real      at;
	u8bit     OccludeLevel;//after
	int       PreviewIndex;
};

STRUCTURE(tnsRenderLine) {
	nListItem      Item;
	tnsRenderVert *L, *R;
	tnsRenderTriangle *TL, *TR;
	nListHandle    Segments;
	//tnsEdge*       Edge;//should be edge material
	//tnsRenderTriangle* Testing;//Should Be tRT** Testing[NumOfThreads]
	char            MinOcclude;
	tns3DObject*    ObjectRef;
	//char            IgnoreConnectedFace;
	//char            CullStatus;
};

STRUCTURE(tnsBoundingArea) {
	real        L, R, U, B;
	real        CX, CY;
	
	tnsBoundingArea* Child;//1,2,3,4 quadrant

	nListHandle LP;
	nListHandle RP;
	nListHandle UP;
	nListHandle BP;

	int         TriangleCount;
	nListHandle AssociatedTriangles;
};

STRUCTURE(tnsRenderSubPixel) {
	real               Depth;
	tnsRenderTriangle* BelongTo;
	tnsVector3d        Weight;  //belongto->vp 1 2 3
};

STRUCTURE(tnsRenderTile) {
	int                Row, Column;
	int                SubX, SubY, SubXLim, SubYLim;//lower Left Corner As 0
	real               FX, FY, FXLim, FYLim;  //ratio;
	tnsRenderSubPixel* FirstPixel;            //lower Left Corner As 0
	nListHandle        AssociatedTriangles;   //lstptrs
	nListHandle        AssociatedLines;       //lstptrs
	char               Rendered;
};

STRUCTURE(tnsFrameBuffer) {
	char Offset;
	int W, H;
	int SubPixelSample;//1,2,3,4, Use Squared Value.
	int TileSizeW, TileSizeH;
	int TileCountX, TileCountY;
	real  WidthPerTile, HeightPerTile;
	tnsRenderSubPixel* Pixels;
	tnsRenderTile*  Tiles;//[Row*CC+Colum]
	tnsMatrix44d    ViewProjection;
	tnsMatrix44d    VPInverse;

	nSafeString*    OutputFolder;//end with slash;
	nSafeString*    ImagePrefix;
	nSafeString*    ImageNameConnector;

	int             OutputMode;
	int             OutputAALevel;
};

#define TNS_OUTPUT_MODE_COMBINED  0
#define TNS_OUTPUT_MODE_PER_LAYER 1

#define TNS_OUTPUT_AA_1 1
#define TNS_OUTPUT_AA_2 2
#define TNS_OUTPUT_AA_4 4
#define TNS_OUTPUT_AA_8 8
#define TNS_OUTPUT_AA_16 16


STRUCTURE(tnsRenderBufferPreviewNode) {
	u8bit           Used;
	u8bit           Level;
	u32bit          VertCount;
	GLuint          VBO;
};

#define TNS_THREAD_LINE_COUNT 10000

#define TNS_OVERRIDE_DISPLAY_NULL 0
#define TNS_OVERRIDE_DISPLAY_SHOW 1
#define TNS_OVERRIDE_DISPLAY_HIDE 2
#define TNS_OVERRIDE_DISPLAY_SWAP 3

STRUCTURE(tnsRenderBuffer) {
	nHyperItem      Item;

	nSafeString*    Name;

	tnsFrameBuffer* FrameBuffer;

	tnsBoundingArea* InitialBoundingAreas;
	u32bit           BoundingAreaCount;
	u32bit           BaVBO;
	//u32bit           BaFillVBO;

	nListHandle     VertexBufferPointers;
	nListHandle     TriangleBufferPointers;
	
	nListHandle     AllRenderLines;

	//nListHandle     IntersectingVertexBuffer;


	nStaticMemoryPool RenderDataPool;

	//render status

	tnsVector3d     ViewVector;

	int             TriangleSize;

	u32bit          ContourCount;
	u32bit          ContourProcessed;
	nListItemPointer* ContourManaged;
	nListHandle     Contours;

	u32bit          IntersectionCount;
	u32bit          IntersectionProcessed;
	nListItemPointer* IntersectionManaged;
	nListHandle     IntersectionLines;

	u32bit          CreaseCount;
	u32bit          CreaseProcessed;
	nListItemPointer* CreaseManaged;
	nListHandle     CreaseLines;

	u32bit          MaterialLineCount;
	u32bit          MaterialProcessed;
	nListItemPointer* MaterialManaged;
	nListHandle     MaterialLines;

	CRITICAL_SECTION csInfo;
	CRITICAL_SECTION csData;
	CRITICAL_SECTION csManagement;

	//settings

	int             OutputTransparent;
	tnsVector4d     BackgroundColor;

	int             MaxOccludeLevel;
	real            CreaseAngle;
	real            CreaseCos;
	int             CreaseAllowOverride;
	int             ThreadCount;
	
	real            OverallProgress;
	int             CalculationStatus;

	int             DrawMaterialPreview;
	real            MaterialTransparency;

	int             ShowLine;
	int             ShowFast;
	int             ShowMaterial;
	int             OverrideDisplay;

	nListHandle     DrawCommands;

	tnsRenderBufferPreviewNode RenderPreview[32];

	tnsScene*  Scene;
	//tnsCamera* Camera;

	//tnsRenderTriangles are in mesh object.
};

#define TNS_COMMAND_LINE     0
#define TNS_COMMAND_MATERIAL 1
#define TNS_COMMAND_EDGE     2


#define TNS_TRANSPARENCY_DRAW_SIMPLE  0
#define TNS_TRANSPARENCY_DRAW_LAYERED 1

#define TNS_OVERRIDE_ONLY                     0
#define TNS_OVERRIDE_EXCLUDE                  1
//#define TNS_OVERRIDE_ALL_OTHERS_OUTSIDE_GROUP 2
//#define TNS_OVERRIDE_ALL_OTHERS_IN_GROUP      3
//#define TNS_OVERRIDE_ALL_OTHERS               4

STRUCTURE(tnsRenderDrawCommand) {
	nHyperLevel1     Item;
	tnsRenderBuffer* ParentRB;

	int              Type;

	nSafeString*     Name;

	tnsVector4d Color;
	real        Thickness;
	int         OccludeBegin, OccludeEnd;

	int         UseStipple;
	u16bit      StipplePattern;
	u8bit       StippleSize;
	int         DrawThisCommand;

	int         DrawContour;
	int         DrawCrease;
	int         DrawIntersections;
	int         DrawMaterialLines;

	GLuint      VBO;
	GLuint      NBO;
	int         VertCount;

	int          OverrideColor;
	tnsMaterial* MaterialRef;
	nSafeString* ReadMaterialName;

	tnsGroup*    OverrideGroup;
	nSafeString* ReadGroupName;
	int          ExcludeGroup;

	real         NormalEdgeClamp;
	real         NormalEdgeStrength;
	real         DepthEdgeClamp;
	real         DepthEdgeStrength;

	int          ClearDepthBuffer;

	int          DepthTest;

	int          TransparencyMode;
	real         Transparency;
};

#define TNS_CALCULATION_IDLE         0
#define TNS_CALCULATION_GEOMETRY     1
#define TNS_CALCULATION_CONTOUR      2
#define TNS_CALCULATION_INTERSECTION 3
#define TNS_CALCULATION_OCCLUTION    4
#define TNS_CALCULATION_FINISHED     100


STRUCTURE(tnsRenderTaskInfo) {
	thrd_t           ThreadHandle;

	tnsRenderBuffer* RenderBuffer;
	int              ThreadID;

	nListItemPointer* Contour;
	nListHandle       ContourPointers;

	nListItemPointer* Intersection;
	nListHandle       IntersectionPointers;

	nListItemPointer* Crease;
	nListHandle       CreasePointers;

	nListItemPointer* Material;
	nListHandle       MaterialPointers;

};





#define TNS_TILE(tile,r,c,CCount)\
tile[r*CCount+c]

#define TNS_CLAMP(a,Min,Max)\
a=a<Min?Min:(a>Max?Max:a)

#define TNS_MAX2(a,b)\
(a>b?a:b)

#define TNS_MIN2(a,b)\
(a<b?a:b)

#define TNS_MAX3(a,b,c)\
(a>TNS_MAX2(b,c)?a:TNS_MAX2(b,c))

#define TNS_MIN3(a,b,c)\
(a<TNS_MIN2(b,c)?a:TNS_MIN2(b,c))

#define TNS_MAX2_INDEX(a,b)\
(a>b?0:1)

#define TNS_MIN2_INDEX(a,b)\
(a<b?0:1)

#define TNS_MAX3_INDEX(a,b,c)\
(a>b?(a>c?0:(b>c?1:2)):(b>c?1:2))

#define TNS_MIN3_INDEX(a,b,c)\
(a<b?(a<c?0:(b<c?1:2)):(b<c?1:2))

#define TNS_MAX3_INDEX_ABC(x,y,z)\
(x>y?(x>z?a:(y>z?b:c)):(y>z?b:c))

#define TNS_MIN3_INDEX_ABC(x,y,z)\
(x<y?(x<z?a:(y<z?b:c)):(y<z?b:c))

#define TNS_ABC(index)\
(index==0?a:(index==1?b:c))


#define TNS_DOUBLE_CLOSE_ENOUGH(a,b)\
(((a)+DBL_EDGE_LIM)>=(b) && ((a)-DBL_EDGE_LIM)<=(b))

//#define TNS_DOUBLE_CLOSE_ENOUGH(a,b)\
//(((a)+0.00000000001)>=(b) && ((a)-0.0000000001)<=(b))


#define TNS_FLOAT_CLOSE_ENOUGH_WIDER(a,b)\
(((a)+0.0000001)>=(b) && ((a)-0.0000001)<=(b))



#define TNS_FRAMEBUFFER_PIXEL(FrameBuffer,Row,Column)\
&((FrameBuffer)->Pixels[Row*FrameBuffer->TileSizeW*FrameBuffer->W*FrameBuffer->SubPixelSample + Column*FrameBuffer->H*FrameBuffer->TileSizeH*FrameBuffer->SubPixelSample])

#define TNS_IN_TILE_X(RenderTile,Fx)\
(RenderTile->FX<=Fx && RenderTile->FXLim>=Fx)

#define TNS_IN_TILE_Y(RenderTile,Fy)\
(RenderTile->FY<=Fy && RenderTile->FYLim>=Fy)


#define TNS_IN_TILE(RenderTile,Fx,Fy)\
(TNS_IN_TILE_X(RenderTile,Fx) && TNS_IN_TILE_Y(RenderTile,Fy))


#define TNS_RENDERBUFFER_INCOMPLETE            0
#define TNS_RENDERBUFFER_GEOMETRY_COMPLETE     1
#define TNS_RENDERBUFFER_RASTERIZER_COMPLETE   2
#define TNS_RENDERBUFFER_COMPLETE  (TNS_RENDERBUFFER_GEOMETRY_COMPLETE|TNS_RENDERBUFFER_RASTERIZER_COMPLETE)


#define TNS_DISPLAY_MODE_WIRE     GL_LINE_LOOP
#define TNS_DISPLAY_MODE_SOLID    GL_TRIANGLES
#define TNS_INTERNAL_MODE_WIRE    1
#define TNS_INTERNAL_MODE_SOLID   2

#define TNS_ARROW_L 1
#define TNS_ARROW_R 2
#define TNS_ARROW_U 3
#define TNS_ARROW_D 4


#define TNS_KEYWORD_OBJECT_COUNT          1
#define TNS_KEYWORD_OBJECT_NAME           2
#define TNS_KEYWORD_OBJECT_TYPE           3
#define TNS_KEYWORD_OBJECT_LOCATION       4
#define TNS_KEYWORD_OBJECT_ROTATION       5
#define TNS_KEYWORD_OBJECT_SCALE          6
#define TNS_KEYWORD_OBJECT_CHILDREN_COUNT 7
#define TNS_KEYWORD_OBJECT_CHILDREN_NAMES 8
//#define TNS_KEYWORD_MESH_COUNT            9
//#define TNS_KEYWORD_MESH_NAME             10
#define TNS_KEYWORD_MESH_VERTEX_COUNT     11
#define TNS_KEYWORD_MESH_VERTICES         12
#define TNS_KEYWORD_MESH_LOOP_COUNT       13
#define TNS_KEYWORD_MESH_TOPOLOGY         14
#define TNS_KEYWORD_MESH_UV_COUNT         15
#define TNS_KEYWORD_MESH_UV_NAME          16
#define TNS_KEYWORD_MATERIAL_COUNT        17
#define TNS_KEYWORD_MATERIAL_NAME         18
#define TNS_KEYWORD_MATERIAL_COLOR        19
#define TNS_KEYWORD_MATERIAL_END          20

#define TNS_KEYWORD_CAMERA_FOV            20
#define TNS_KEYWORD_CAMERA_IS_ACTIVE      21
#define TNS_KEYWORD_CAMERA_NEAR           22
#define TNS_KEYWORD_CAMERA_FAR            23

#define TNS_KEYWORD_GROUP_COUNT           24
#define TNS_KEYWORD_GROUP_NAME            25
#define TNS_KEYWORD_OBJECT_GROUP_COUNT    26




void tnsRegisterAllActuators();

void tnsSetuptnsFontManager();

tnsShader* tnsNewShaderProgram(int VertexShaderID, int FragmentShaderID, int GeometryShaderID);
int tnsNewFragmentShaderStringBased(char* Content);
int tnsNewVertexShaderStringBased(char* Content);

int tnsNextPowOf2(int i);

int tnsNewShader(char* vtsfilename, char* fgsfilename, int* result_index);
int tnsNewShaderdv(char* vtsfilename, char* fgsfilename, tnsShader** result_shader);
int tnsEnableShader(int index);
int tnsEnableShaderv(tnsShader* shader);
int tnsUseShader(tnsShader* shader);
void tnsUseUiShader();
void tnsUseStringShader();
void tnsUseTextureShader();
void tnsUseMultiplyTextureShader();
void tnsUseTextureMultisampleShader();
void tnsUseTextureMultisampleAlphaShader();
void tnsUseTransparentGridShader();
void tnsUseTextureSobelColorMSShader();
void tnsUseTextureSobelMSShader();
void tnsUseFixedShader();

void tnsUseImagePeelShader();


void tnsInitRenderKernel(int matrixStackLevel);
void tnsInitBuiltinShaders();
void tnsInitWindowDefaultRenderConfig();
void tnsQuit();

real* tnsGetModelMatrix();
real* tnsGetViewMatrix();
real* tnsGetProjectionMatrix();
void tnsResetModelMatrix();
void tnsResetViewMatrix();
void tnsResetProjectionMatrix();

void tnsOrtho(real xMin, real xMax, real yMin, real yMax, real zMin, real zMax);
void tnsPerspective(real fFov_rad, real fAspect, real zMin, real zMax);

void tnsPopMatrix();
void tnsPushMatrix();

void tnsTranslate3d(real x, real y, real z);
void tnsPreTranslate3d(real x, real y, real z);
void tnsRotate4d(real degrees, real x, real y, real z);
void tnsPreRotate4d(real degrees, real x, real y, real z);
void tnsScale3d(real x, real y, real z);
void tnsPreScale3d(real x, real y, real z);

void tnsColor4d(real r, real g, real b, real a);
void tnsColor4dv(real* rgba);

void tnsPolygonMode(GLenum PolyMode);
void tnsShadeMode(GLenum ShadeMode);


tnsBatch* tnsCreateBatch(u32bit NumVert, int Dimention, float* Data, float* Normal);
tnsBatch* tnsCreateBatchi(u32bit NumVert, int Dimention, int* Data);
tnsBranchedCommand* tnsCreateCommand(tnsBatch* b, u32bit ElementCount, int Dimention, GLenum DrawAs, int InternalMode, u32bit* Elements);
void tnsDeleteBatch(tnsBatch* b);

void tnsVertex3d(real x, real y, real z);
void tnsVertex2d(real x, real y);
void tnsVertexArray2d(real* verts, int amount);
void tnsVertexArray3d(real* verts, int amount);
void tnsColorArray4d(real* colors, int amount);
void tnsNormalArray3d(real* normals, int amount);
void tnsTexCoordArray2d(real* coords, int amount);
void tnsIndexArray(GLushort* index, short amount);
void tnsPackAs(GLenum Mode);
void tnsFlush();

#define tnsLinearItp(L,R,T)\
((L)*(1.0f - (T)) + (R)*(T))

__inline void tMatConvert44df(tnsMatrix44d from, tnsMatrix44f to) {
	to[0] = from[0];
	to[1] = from[1];
	to[2] = from[2];
	to[3] = from[3];
	to[4] = from[4];
	to[5] = from[5];
	to[6] = from[6];
	to[7] = from[7];
	to[8] = from[8];
	to[9] = from[9];
	to[10] = from[10];
	to[11] = from[11];
	to[12] = from[12];
	to[13] = from[13];
	to[14] = from[14];
	to[15] = from[15];
}

__inline int tRdrTrangleLineBoundBoxTest(tnsRenderTriangle* rt, tnsRenderLine* rl) {
	if (TNS_MAX3(rt->V[0]->FrameBufferCoord[2], rt->V[1]->FrameBufferCoord[2], rt->V[2]->FrameBufferCoord[2]) > TNS_MIN2(rl->L->FrameBufferCoord[2], rl->R->FrameBufferCoord[2])) return 0;
	if (TNS_MAX3(rt->V[0]->FrameBufferCoord[0], rt->V[1]->FrameBufferCoord[0], rt->V[2]->FrameBufferCoord[0]) < TNS_MIN2(rl->L->FrameBufferCoord[0], rl->R->FrameBufferCoord[0])) return 0;
	if (TNS_MIN3(rt->V[0]->FrameBufferCoord[0], rt->V[1]->FrameBufferCoord[0], rt->V[2]->FrameBufferCoord[0]) > TNS_MAX2(rl->L->FrameBufferCoord[0], rl->R->FrameBufferCoord[0])) return 0;
	if (TNS_MAX3(rt->V[0]->FrameBufferCoord[1], rt->V[1]->FrameBufferCoord[1], rt->V[2]->FrameBufferCoord[1]) < TNS_MIN2(rl->L->FrameBufferCoord[1], rl->R->FrameBufferCoord[1])) return 0;
	if (TNS_MIN3(rt->V[0]->FrameBufferCoord[1], rt->V[1]->FrameBufferCoord[1], rt->V[2]->FrameBufferCoord[1]) > TNS_MAX2(rl->L->FrameBufferCoord[1], rl->R->FrameBufferCoord[1])) return 0;
	return 1;
}

double tMatGetLinearRatio(real L, real R, real FromL);
__inline int tRdrLineIntersectTest2d(tnsVector2d a1, tnsVector2d a2, tnsVector2d b1, tnsVector2d b2, double* aRatio)  {
	double k1, k2;
	double x;
	double y;
	double Ratio;
	double xDiff = (a2[0] - a1[0]);// +DBL_EPSILON;
	double xDiff2 = (b2[0] - b1[0]);

	if (xDiff == 0) {
		if (xDiff2 == 0) {
			*aRatio = 0;
			return 0;
		}
		double r2 = tMatGetLinearRatio(b1[0], b2[0], a1[0]);
		y = tnsLinearItp(b1[1], b2[1], r2);
		*aRatio = Ratio = tMatGetLinearRatio(a1[1],a2[1],y);
	} else {
		if (xDiff2 == 0) {
			Ratio = tMatGetLinearRatio(a1[0], a2[0], b1[0]);
			//y = tnsLinearItp(a1[1], a2[1], r2);
			*aRatio = Ratio;
		}else {
			k1 = (a2[1] - a1[1]) / xDiff;
			k2 = (b2[1] - b1[1]) / xDiff2;

			if ((k1 == k2)) 
				return 0;

			x = (a1[1] - b1[1] - k1*a1[0] + k2*b1[0]) / (k2 - k1);

			Ratio = (x - a1[0]) / xDiff;

			*aRatio = Ratio;
		}
	}

	

	if (b1[0] == b2[0]) {
		y = tnsLinearItp(a1[1], a2[1], Ratio);
		if (y > TNS_MAX2(b1[1], b2[1]) || y < TNS_MIN2(b1[1], b2[1])) return 0;
	}else
		if(Ratio <= 0 || Ratio>1 ||
		(b1[0]>b2[0] && x>b1[0]) ||
		(b1[0]<b2[0] && x<b1[0]) ||
		(b2[0]>b1[0] && x>b2[0]) ||
		(b2[0]<b1[0] && x<b2[0]))		
		return 0;

	return 1;
}
__inline double tRdrGetLineZ(tnsVector3d L, tnsVector3d R, real Ratio) {
	//double z = 1 / tnsLinearItp(1 / L[2], 1 / R[2], Ratio);
	double z = tnsLinearItp(L[2], R[2], Ratio);
	return z;
}
__inline double tRdrGetLineZPoint(tnsVector3d L, tnsVector3d R, tnsVector3d FromL) {
	double r = (FromL[0] - L[0]) / (R[0] - L[0]);
	return tnsLinearItp(L[2], R[2], r);
	//return 1 / tnsLinearItp(1 / L[2], 1 / R[2], r);
}
__inline double tRdrGetLinearRatio(tnsVector3d L, tnsVector3d R, tnsVector3d FromL) {
	double r = (FromL[0] - L[0]) / (R[0] - L[0]);
	return r;
}

__inline double tMatGetLinearRatio(real L, real R, real FromL) {
	double r = (FromL - L) / (R - L);
	return r;
}
__inline void tMatVectorMinus2d(tnsVector2d result, tnsVector2d l, tnsVector2d r) {
	result[0] = l[0] - r[0];
	result[1] = l[1] - r[1];
}

__inline void tMatVectorMinus3d(tnsVector3d result, tnsVector3d l, tnsVector3d r) {
	result[0] = l[0] - r[0];
	result[1] = l[1] - r[1];
	result[2] = l[2] - r[2];
}
__inline void tMatVectorSubtract3d(tnsVector3d l, tnsVector3d r) {
	l[0] = l[0] - r[0];
	l[1] = l[1] - r[1];
	l[2] = l[2] - r[2];
}
__inline void tMatVectorPlus3d(tnsVector3d result, tnsVector3d l, tnsVector3d r) {
	result[0] = l[0] + r[0];
	result[1] = l[1] + r[1];
	result[2] = l[2] + r[2];
}
__inline void tMatVectorAccum3d(tnsVector3d l, tnsVector3d r) {
	l[0] = l[0] + r[0];
	l[1] = l[1] + r[1];
	l[2] = l[2] + r[2];
}
__inline void tMatVectorAccum2d(tnsVector2d l, tnsVector2d r) {
	l[0] = l[0] + r[0];
	l[1] = l[1] + r[1];
}
__inline void tMatVectorNegate3d(tnsVector3d result, tnsVector3d l) {
	result[0] = -l[0];
	result[1] = -l[1];
	result[2] = -l[2];
}
__inline void tMatVectorNegateSelf3d(tnsVector3d l) {
	l[0] = -l[0];
	l[1] = -l[1];
	l[2] = -l[2];
}
__inline void tMatVectorCopy2d(tnsVector2d from, tnsVector2d to) {
	to[0] = from[0];
	to[1] = from[1];
}
__inline void tMatVectorCopy3d(tnsVector3d from, tnsVector3d to) {
	to[0] = from[0];
	to[1] = from[1];
	to[2] = from[2];
}
__inline void tMatVectorCopy4d(tnsVector4d from, tnsVector4d to) {
	to[0] = from[0];
	to[1] = from[1];
	to[2] = from[2];
	to[3] = from[3];
}
__inline void tMatVectorMultiSelf4d(tnsVector3d from, real num) {
	from[0] *= num;
	from[1] *= num;
	from[2] *= num;
	from[3] *= num;
}
__inline void tMatVectorMultiSelf3d(tnsVector3d from, real num) {
	from[0] *= num;
	from[1] *= num;
	from[2] *= num;
}
__inline void tMatVectorMultiSelf2d(tnsVector3d from, real num) {
	from[0] *= num;
	from[1] *= num;
}

__inline real tMatDirectionToRad(tnsVector2d Dir) {
	real arcc = acos(Dir[0]);
	real arcs = asin(Dir[1]);

	if (Dir[0] >= 0) {
		if(Dir[1]>=0) return arcc;
		else return TNS_PI * 2 - arcc;
	}else {
		if (Dir[1] >= 0) return arcs + TNS_PI / 2;
		else return TNS_PI + arcs;
	}
}


#define TNS_ATLAS_DEFAULT_INPUT_WIDTH 2048


void tSdrApplyProjectionInverse(tnsShader* ts, tnsMatrix44d m);
void tSdrApplyProjection(tnsShader* ts, tnsMatrix44d m);
void tSdrApplyModel(tnsShader* ts, tnsMatrix44d m);
void tSdrApplyView(tnsShader* ts, tnsMatrix44d m);
void tSdrApplyNormalScaler(tnsShader* ts, tnsMatrix44d m);



real tMatDistIdv2(real x1, real y1, real x2, real y2);
real tMatDist3dv(tnsVector3d l, tnsVector3d r);
real tMatDist2dv(tnsVector2d l, tnsVector2d r);

real tMatLength3d(tnsVector3d l); real tMatLength2d(tnsVector3d l);
void tMatNormalize3d(tnsVector3d result, tnsVector3d l);
void tMatNormalizeSelf3d(tnsVector3d result);
real tMatDot3d(tnsVector3d l, tnsVector3d r, int normalize);
real tMatDot2d(tnsVector2d l, tnsVector2d r, int normalize);
real tMatVectorCross3d(tnsVector3d result, tnsVector3d l, tnsVector3d r);
real tMatAngleRad3d(tnsVector3d from, tnsVector3d to, tnsVector3d PositiveReference);
void tMatApplyRotation33d(tnsVector3d result, tnsMatrix44d mat, tnsVector3d v);
void tMatApplyRotation43d(tnsVector3d result, tnsMatrix44d mat, tnsVector3d v);
void tMatApplyTransform43d(tnsVector3d result, tnsMatrix44d mat, tnsVector3d v);
void tMatApplyNormalTransform43d(tnsVector3d result, tnsMatrix44d mat, tnsVector3d v);
void tMatApplyTransform44d(tnsVector4d result, tnsMatrix44d mat, tnsVector4d v);
void tMatApplyTransform44dTrue(tnsVector4d result, tnsMatrix44d mat, tnsVector4d v);


void tMatLoadIdentity44d(tnsMatrix44d m);
void tMatMakeOrthographicMatrix44d(tnsMatrix44d mProjection, real xMin, real xMax, real yMin, real yMax, real zMin, real zMax);
void tMatMakePerspectiveMatrix44d(tnsMatrix44d mProjection, real fFov_rad, real fAspect, real zMin, real zMax);
void tMatMakeTranslationMatrix44d(tnsMatrix44d mTrans, real x, real y, real z);
void tMatMakeRotationMatrix44d(tnsMatrix44d m, real angle_rad, real x, real y, real z);
void tMatMakeScaleMatrix44d(tnsMatrix44d m, real x, real y, real z);
void tMatMakeViewportMatrix44d(tnsMatrix44d m, real w, real h);
void tMatMultiply44d(tnsMatrix44d result, tnsMatrix44d l, tnsMatrix44d r);
void tMatMakeRotationXMatrix44d(tnsMatrix44d m, real angle_rad);
void tMatMakeRotationYMatrix44d(tnsMatrix44d m, real angle_rad);
void tMatMakeRotationZMatrix44d(tnsMatrix44d m, real angle_rad);
void tMatRemoveTranslation44d(tnsMatrix44d result, tnsMatrix44d mat);
void tMatClearTranslation44d(tnsMatrix44d mat);

real tMatAngleRad3d(tnsVector3d from, tnsVector3d to, tnsVector3d PositiveReference);
real tMatLength3d(tnsVector3d l);
void tMatNormalize2d(tnsVector2d result, tnsVector2d l);
void tMatNormalize3d(tnsVector3d result, tnsVector3d l);
void tMatNormalizeSelf3d(tnsVector3d result);
real tMatDot3d(tnsVector3d l, tnsVector3d r, int normalize);
real tMatVectorCross3d(tnsVector3d result, tnsVector3d l, tnsVector3d r);
void tMatVectorCrossOnly3d(tnsVector3d result, tnsVector3d l, tnsVector3d r);


void tMatExtractXYZEuler44d(tnsMatrix44d mat, real* x_result, real* y_result, real* z_result);

void tMatPrintMatrix44d(tnsMatrix44d l);




int tRdrPointInsideTriangle3d(tnsVector3d v, tnsVector3d v0, tnsVector3d v1, tnsVector3d v2);



void tObjExtractSelfEulerRotation(tns3DObject* o);
void tObjExtractSelfLocation(tns3DObject* o);
void tObjExtractSelfScaling(tns3DObject* o);
void tObjApplyGlobalTransformMatrix(tns3DObject* o);
void tObjApplyGlobalAndSelfTransformMatrix(tns3DObject* o);
void tObjApplyGlobalTransformMatrixRecursive(tns3DObject* o);
void tObjApplyGlobalTransformMatrixReverted(tns3DObject* o);
void tObjApplySelfTransformMatrix(tns3DObject* o, real* inverse_result);
void tObjApplySelfTransformValue(tns3DObject* o);
void tObjApplyGlobalTransformValue(tns3DObject* o);
void tObjApplySameGlobalTransform(tns3DObject* to, tns3DObject* from);
void tObjGetSelfTransformFromGlobal(tns3DObject* o);
void tObjAssignListLinker(tns3DObject* o, tnsScene* Scene);
void tObjDestroyListLinker(tns3DObject* o, tnsScene* Scene);
void tObjInitObjectBase(tns3DObject* o, tnsScene* Scene, char* Name, int Type,
	real AtX, real AtY, real AtZ,
	real RotX, real RotY, real RotZ, real RotW, u8bit RotationMode,
	real Scale);

void tObjSimpleTriangulateRender(tnsFace* F, tnsRenderTriangle* TBuf, int TriangleSize, tnsRenderVert* VBuf, nListHandle* LBuf, tnsRenderTriangle** Next, nMemoryPool* RenderData);


void tObjRefreshBezierIndex(tnsBezierObject* bo);
void tObjRefreshGeometeryIndex(tnsMeshObject* mo);


tnsMeshObject* tnsLoadObjectFromFile(char* FileName);

int tnsLoadExchange(char* FileName);


void tnsDrawMeshObject(tnsMeshObject* mo, int OverrideDisplayMode);

tns3DObject* tnsFindRootObject(tnsScene* s, char* Name);
tns3DObject* tnsFindObject(tnsScene* s, char* Name, tns3DObject** From);

void tnsRotateObjectSelf(tns3DObject* o, real x, real y, real z);
void tnsRotateObjectGlobalXYZ(tns3DObject* o, real x, real y, real z);
void tnsRotateObjectGlobalBaseXYZ(tns3DObject* o, tnsVector3d* BaseLocation, real x, real y, real z);
void tnsTranslateObjectGlobal(tns3DObject* o, real x, real y, real z);
void tnsTranslateObjectLocal(tns3DObject* o, real x, real y, real z);
void tnsMoveObjectPrior(tns3DObject* o, real x, real y, real z);
void tnsMoveObjectGlobal(tns3DObject* o, real x, real y, real z);
void tnsScaleObject(tns3DObject* o, real fac);
void tnsSetObjectLocation(tns3DObject* o, real x, real y, real z);
void tnsParentObject(tns3DObject* child, tns3DObject* parent);
void tnsUnparentObject(tns3DObject* o);
void tnsZoomViewingCamera(tnsCamera* c, real Ratio);
void tnsRotateViewingCamera(tnsCamera* c, real x, real z);
void tnsTranslateViewingCamera(tnsCamera* c, int ViewportW, int ViewportH, real x, real y);


tnsScene* tnsCreateScene(char* name);
void tnsDestroyScene(tnsScene* s);
void tnsActiveScene(tnsScene* s);
tnsCamera* tnsCreateXYZCamera(tnsScene* Scene, char* Name, real FOV,
	real AtX, real AtY, real AtZ,
	real RotX, real RotY, real RotZ,
	real FocusDistance);
tns3DObject* tnsCreateXYZEmpty(tnsScene* Scene, char* Name,
	real AtX, real AtY, real AtZ);

tns3DObject* tnsCreateMeshObjectPlane(tnsScene* Scene, char* Name, real AtX, real AtY, real AtZ);

void tnsApplyCameraView(int W, int H, tnsCamera* Camera);
void tnsSetActiveCamera(tnsScene* s, tnsCamera* Camera);


void tnsApplyObjectMatrix(tns3DObject* Object);
void tnsDestroyObjectRecursive(tns3DObject* Object, tnsScene* Scene);

void tnsTrackZTo(tns3DObject* o, tnsVector3d Target);


void tnsDrawThisObject(tns3DObject* o, tnsScene* s);
void tnsDrawObjectTree(tnsScene* Scene, tns3DObject* from);
void tnsDrawScene(int W, int H, tnsScene* Scene);
void tnsDrawWorld(int W, int H);



tnsMaterial* tnsCreateMaterial(tnsScene* s, char* Name);
tnsMaterial* tnsFindMaterialByIndex(tnsScene* s, int index);
tnsMaterial* tnsFindMaterial(tnsScene* s, char* name);
__inline tnsMaterial* tnsGetIndexedMaterial(tnsScene* s, int index) {
	if (s->MaterialPointers[index]) return s->MaterialPointers[index];
	else return (s->MaterialPointers[index] = tnsFindMaterialByIndex(s, index));
}


tnsGroup* tnsCreateGroup(tnsScene* s, char* Name);
tnsGroup* tnsFindGroup(tnsScene* s, char* Name);
void tnsAddToGroup(tns3DObject* o, tnsGroup* g);
void tnsRemoveFromGroup(tns3DObject* o, tnsGroup* g);
int tnsGroupHaveObject(tnsGroup* g, tns3DObject* o);
int tnsObjectInGroup(tns3DObject* o, tnsGroup* g);
void tnsObjectGroupMarkFinished(tns3DObject* o);
int tnsObjectGroupIsFinished(tns3DObject* o);
void tnsCleanObjectFinishMarks(tnsScene* s);

void tnsClearAll();
void tnsClearColorv(real* rgba);
void tnsClearColor(real r, real g, real b, real a);
void tnsSwitchToCurrentWindowContext(void* wnd);

void tObjRecalculateVertNormal(tnsVert* v);
void tObjRecalculateFaceAverageNormal(tnsFace* f);


//tnsSlicedMesh* tnsCreateEmptyMesh(tnsScene* Scene, char* Name, real FOV,
//	real AtX, real AtY, real AtZ,
//	real RotX, real RotY, real RotZ);
//tnsVert* tnsCreateVert(tnsSlicedMesh* sm, real x, real y, real z);
//tnsVert* tnsCreateVertV(tnsSlicedMesh* sm, real* xyz);
//int tnsVertInFace(tnsFace* f, tnsVert* V);
//tnsEdge* tnsCreateEdge(tnsSlicedMesh* sm, tnsVert* L, tnsVert* R);
//void tnsDestroyEdge(tnsSlicedMesh* sm, tnsEdge* e);
//
//void tnsSelectVertOnly(tnsSlicedMesh* sm, tnsVert* V);
//void tnsDeselectVertOnly(tnsSlicedMesh*m, tnsVert* V);
//void tnsSelectEdgeOnly(tnsEdge* edge);
//void tnsDeselectEdgeOnly(tnsEdge* edge);
//void tnsSelectFaceOnly(tnsFace* face);
//void tnsDeselectFaceOnly(tnsEdge* face);
//
//void tnsScanAndSelectEdgeInLoop(tnsSlicedMesh* sm, nListHandle* loop, int* AllSelected);
//void tnsScanAndSelectEdgeAndFace(tnsSlicedMesh* sm);
//void tnsSelectVert(tnsSlicedMesh* sm, tnsVert* V);
//void tnsDeselectVert(tnsSlicedMesh* sm, tnsVert* V);
//void tnsSelectEdge(tnsSlicedMesh* sm, tnsEdge* e);
//void tnsDeselectEdge(tnsSlicedMesh* sm, tnsEdge* e);
//void tnsSelectFace(tnsSlicedMesh* sm, tnsFace* f);
//void tnsDeselectFace(tnsSlicedMesh* sm, tnsFace* f);
//void tnsSelectSlice(tnsSlicedMesh* sm, tnsSlice* s);
//void tnsDeselectSlice(tnsSlicedMesh* sm, tnsSlice* s);
//
//real tnsCalcEdgeRelativeAngle(tnsEdge* E1, tnsEdge* E2);
//void tnsExtrudeSelection(tnsSlicedMesh* sm, int Mode, real Depth, int Dir);
//void tnsDestroyEdgeListNode(tnsSlicedMesh* sm, tnsFace* F, tnsEdgeListNode* eln);
//tnsEdgeListNode* tnsCreateEdgeListNode(tnsSlicedMesh* sm, tnsFace* F, tnsVert* V1, tnsVert* V2);
//tnsFace* tnsCreateFaceLoop(tnsSlicedMesh* sm, tnsVert* V1, tnsVert* V2, tnsVert* V3);
//tnsFace* tnsCreateFaceLoopElementList(tnsSlicedMesh* sm, nListHandle* VertElem);
//tnsFace* tnsExtendFaceLoop(tnsSlicedMesh* sm, tnsFace* F, tnsVert* V);
//tnsFace* tnsExtendFaceLoopElementList(tnsSlicedMesh* sm, tnsFace* F, nListHandle* VertElem);
//tnsFace* tnsMakeQuadDirectFromTwoEdges(tnsSlicedMesh* sm, tnsEdge* E1, tnsEdge* E2);
//void tnsBridgeSlices(tnsSlice* S1, tnsSlice* S2);
//tnsSlice* tnsCreateEmptySlice(tnsSlicedMesh* sm);
//tnsEdgeListNode* tnsBeginSlice(tnsSlicedMesh* sm, tnsSlice* s, tnsVert* V1, tnsVert* V2);
//tnsEdgeListNode* tnsExtendSlice(tnsSlicedMesh* sm, tnsSlice* s, tnsVert* V2);
//void tnsCloseSlice(tnsSlicedMesh* sm, tnsSlice* s);
//tnsSlice* tnsCreateRectangularSlice(tnsSlicedMesh* sm,real AtX, real AtY, real AtZ, real SizeX, real SizeY, int Fill);
//
//




int tnsLoadVectorGraphPackage(const char* name, unsigned int size);
int tnsLoadSystemFont(const char* name, const char* IconName, unsigned int size);
int tnsGetTextureMemoryComponetCount(tnsTexture* t);

int tnsInit2DRectTextureAdvanced(tnsTexture* t, GLuint InternalFormat, GLuint DataFormat, GLuint DataType, int w, int h, void* bits);
int tnsInit2DTexture(tnsTexture* t, GLuint BitsFormat, int w, int h, void* bits);
tnsTexture* tnsCreate2DTexture(GLuint BitsFormat, int w, int h, void* bits);
tnsTexture* tnsCreate2DRectTextureAdvanced(GLuint InternalFormat, GLuint DataFormat, GLuint DataType, int w, int h, void* bits);
tnsTexture* tnsCreate2DTextureMultisample(GLuint BitsFormat, int w, int h, void* bits);
tnsTexture* tnsCreate2DDepthBufferMultisample(GLuint BitsFormat, int w, int h, void* bits);
tnsTexture* tnsCreate2DTextureSupersample(GLuint BitsFormat, int w, int h, void* bits);
tnsTexture* tnsCreate2DDepthBufferSupersample(GLuint BitsFormat, int w, int h, void* bits);
void tnsFeed2DRectTextureData(tnsTexture* t, void* data);
void tnsFeed2DTextureData(tnsTexture* t, void* data);
void tnsConfigure2DTexture(tnsTexture* t, GLuint BitsFormat, int w, int h, void* bits);
tnsOffscreen* tnsCreate2DOffscreenBasic(int ColorBitsFormat, int w, int h, void* ColorData);
void tnsCopyScreenTo2DTexture(tnsTexture* target, int x_lower_left, int y_lower_left, int w, int h);
void tnsActiveTexture(GLenum tex);
void tnsBindTexture0(tnsTexture* t);
void tnsUseTexture0(tnsTexture* t);
void tnsUnbindTexture0(tnsTexture* t);
void tnsBindTexture1(tnsTexture* t);
void tnsUseTexture1(tnsTexture* t);
void tnsBindTexture2(tnsTexture* t);
void tnsBindTexture3(tnsTexture* t);
void tnsBindTexture4(tnsTexture* t);
void tnsUseTexture2(tnsTexture* t);
void tnsUnbindTexture1(tnsTexture* t);
void tnsUnbindTexture2(tnsTexture* t);
void tnsUnbindTexture3(tnsTexture* t);
void tnsUnbindTexture4(tnsTexture* t);
void tnsDeleteTexture(tnsTexture* t);
void tnsDraw2DTextureDirectly(tnsTexture* t, int x_upper_right, int y_upper_right, int w, int h);
void tnsDraw2DTextureArg(tnsTexture* t,
	real x_upper_right, real y_upper_right, int w, int h,
	real* MultiplyColor,
	real LPadding, real RPadding, real TPadding, real BPadding);

int tnsTextureMemorySize(tnsTexture* t);
void tnsCreateTextureReadbackBuffer(tnsTexture* t);
void tnsDeleteTextureReadbackBuffer(tnsTexture* t);
void tnsReadbackTexture(tnsTexture* t);

#define TNS_CLAMP_TEXTURE_W(t,Col)\
	{if (Col >= t->Width) Col = t->Width - 1; if (Col < 0) Col = 0;}

#define TNS_CLAMP_TEXTURE_H(t,Row)\
	{if (Row >= t->Height) Row = t->Height - 1;if (Row < 0) Row = 0;}

#define TNS_CLAMP_TEXTURE_CONTINUE(t,Col,Row)\
	{if (Col >= t->Width) continue; if (Col < 0) continue;\
     if (Row >= t->Height) continue; if (Row < 0)continue; }

__inline void* tnsTextureSampleSafe(tnsTexture* t, int Col, int Row) {
	int W = t->Width, H=t->Height;
	TNS_CLAMP_TEXTURE_W(t, Col);
	TNS_CLAMP_TEXTURE_H(t, Row);

	return ((char*)t->TextureReadBack) + (Row*W + Col)*t->ElemSize;
}

__inline int tnsTextureSampleWeightSafe(tnsTexture* t, int Col, int Row, real Size, tnsVector2d Direction) {
	tnsVector2d Dir = { 0 };
	tnsVector2d Weighted = { 0 };
	int r, c;

	if (!t->SamplePtr) return 0;

	tnsTextureSample* ts;

	for (r = Row - Size; r <= Row + Size; r++) {
		for (c = Col - Size; c <= Col + Size; c++) {
			TNS_CLAMP_TEXTURE_CONTINUE(t, c, r);

			if (tMatDistIdv2(r, c, Row, Col) > Size) continue;

			ts = t->SamplePtr[r*t->Width + c];
			if (!ts) continue;

			Weighted[0] = ts->X - Col; Weighted[1] = ts->Y - Row;
			tMatVectorMultiSelf2d(Weighted, (real)ts->Sample / (real)255);
			tMatVectorAccum2d(Dir, Weighted);
		}
	}

	if (tMatLength2d(Dir) < 0.000001) return 0;

	tMatNormalize2d(Direction, Dir);

	return 1;
}

__inline real tnsTextureSampleSafeGaussU8R(tnsTexture* t, tnsFilterKernel* fk, int Col, int Row) {
	int i,j;
	real result=0;
	int c = Col - fk->Size / 2;
	int r = Row - fk->Size / 2;
	for (i = 0; i < fk->Size; i++) {
		for (j = 0; j < fk->Size; j++) {
			result += (*((u8bit*)tnsTextureSampleSafe(t, c + i, r + j)))*fk->Kernel[i+j*fk->Size];
		}
	}
	return result;
}

__inline void tnsTextureSetBitSafeU8R(tnsTexture* t, int Col, int Row, int Bit) {
	int W = t->Width, H = t->Height;
	if (Row >= H) Row = H - 1;
	if (Row < 0) Row = 0;
	if (Col >= W) Col = W - 1;
	if (Col < 0) Col = 0;

	*(((char*)t->TextureReadBack) + (Row*W + Col)*t->ElemSize) = Bit;
}

__inline void tnsTextureSetRGBSafeU8R(tnsTexture* t, int Col, int Row, u8bit r, u8bit g, u8bit b) {
	int W = t->Width, H = t->Height;
	if (Row >= H) Row = H - 1;
	if (Row < 0) Row = 0;
	if (Col >= W) Col = W - 1;
	if (Col < 0) Col = 0;

	u8bit* rgb =(((u8bit*)t->TextureReadBack) + (Row*W + Col)*t->ElemSize);
	rgb[0] = r;
	rgb[1] = g;
	rgb[2] = b;

}

__inline u8bit* tnsTextureSampleU8(tnsTexture* t, int Col, int Row) {
	int W = t->Width;

	return ((char*)t->TextureReadBack) + (Row*W + Col)*t->ElemSize;
}


tnsOffscreen* tnsCreate2DOffscreenWithDepthMultisample(int ColorBitsFormat, int w, int h, void* ColorData, void* DepthData);
tnsOffscreen* tnsCreate2DOffscreenMultisample(int ColorBitsFormat, int w, int h, void* ColorData);
void tnsCreate2DOffscreenMSAttachmentExtraColor(tnsOffscreen* off);
void tnsCreate2DOffscreenMSAttachmentExtraNormal(tnsOffscreen* off);
tnsOffscreen* tnsCreate2DOffscreenWithDepthSupersample(int ColorBitsFormat, int w, int h, void* ColorData, void* DepthData);
tnsOffscreen* tnsCreate2DOffscreenSupersample(int ColorBitsFormat, int w, int h, void* ColorData);


void tnsDelete2DOffscreen(tnsOffscreen* o);
tnsOffscreen* tnsCreateOffscreenHandle();
void tnsAttach2DOffscreenBuffer(tnsOffscreen* target, GLuint attatchment, tnsTexture* use);
void tnsDetach2DOffscreenBuffer(tnsOffscreen* target, GLuint which_attach_point);

void tnsDrawToOffscreen(tnsOffscreen* toff, int HowMany, GLuint* AttachmentArray);
void tnsDrawToExtraColorAttachment(tnsOffscreen* toff);
void tnsDrawToExtraNormalAttachment(tnsOffscreen* toff);
void tnsDrawToAllExtraAttachments(tnsOffscreen* toff);
void tnsDrawToOffscreenOnlyBind(tnsOffscreen* toff, int HowMany, GLuint* AttachmentArray);
void tnsReadFromOffscreen(tnsOffscreen* toff);
void tnsPassColorBetweenOffscreens(tnsOffscreen* from, tnsOffscreen* to,
	GLint srcX0, GLint srcY0, GLint srcX1, GLint srcY1, GLint dstX0, GLint dstY0, GLint dstX1, GLint dstY1, GLenum FilterMode);

void tnsDrawToScreen();

//======================================================[atlas]

NEED_STRUCTURE(n3DViewUiExtra);
void tns_InitAtlasInput(tnsScene* s, int width);
void tns_DestroyAtlasInputs(tnsScene* s);
void tns_InitAtlasOutputs(n3DViewUiExtra* e, int width);
void tns_DestroyAtlasOutputs(n3DViewUiExtra* e);
void tns_FeedAtlasData(tnsScene* s);
void tns_TriggerAtlasTransform(tnsScene* s, n3DViewUiExtra* e);
void tns_MakeAtlasTriggerBatch(tnsScene* s);
void tns_AtlasDrawTransform(n3DViewUiExtra* e);
void tns_AtlasDrawEdgePreview(n3DViewUiExtra* e);


//==================================================[snake]

void tns_CreateSnakeOffscreen(n3DViewUiExtra* e, int W, int H);

void tns_BuildSnakes(n3DViewUiExtra* e);
void tns_DrawSnakes(n3DViewUiExtra* e);

//==============================================[STR]

int tnsUseFont(char* name);
int nulStopFontService();

NEED_STRUCTURE(nThemeState);

tnsFontSingleCharacter* tfntFetchVectorGraphTextureIDW(int ID);

int tnsStringGetWidthU(wchar_t* content, nThemeState* ts, int Count);
int tnsStringGetWidthM(char* content, nThemeState* ts, int Count);
void tnsDrawStringDirectU(wchar_t* content, nThemeState* ts, int L, int R, int T, int RevY);
void tnsDrawStringAutoM(char* content, nThemeState* ts, int L, int R, int T, int RevY, nStringSplitor* Instructions);
void tnsDrawStringAutoU(wchar_t* content, nThemeState* ts, int L, int R, int T, int RevY, nStringSplitor* Instructions);
void tnsDrawStringWithPriority(wchar_t* Label, wchar_t* MajorContent, nThemeState* ts, int L, int R, int T, int RevY);
void tnsDrawVectorGraphPackage(int ID, nThemeState* ts, int L, int R, int T, int B, int Align, int RevY);

void tnsDrawIcon(short ID, nThemeState* ts, int L, int T, int RevY);

void tnsUseFontCoord(real x, real y, real size);

///=================================




void tnsMakeTriangle(real* arr, real x1, real y1, real x2, real y2, real x3, real y3);
void tnsMakeQuad2d(real* arr, real x1, real y1, real x2, real y2, real x3, real y3, real x4, real y4);
void tnsMakeQuad3d(real* arr, real x1, real y1, real z1, real x2, real y2, real z2, real x3, real y3, real z3, real x4, real y4, real z4);
void tnsMakeQuad4d(real* arr, real x1, real y1, real z1, real w1, real x2, real y2, real z2, real w2, real x3, real y3, real z3, real w3, real x4, real y4, real z4, real w4);
void tnsMakeCircle2d(real* arr, int slices, real ctrX, real ctrY, real r);
void tnsMakeArc2d(real* arr, int slices, real ctrX, real ctrY, real r, real rad_begin, real rad_end);

void tnsMakeLinerGradient3d(real* arr, int num_points, real r0, real g0, real b0, real r1, real g1, real b1);
void tnsMakeLinerGradient4d(real* arr, int num_points, real r0, real g0, real b0, real a0, real r1, real g1, real b1, real a1);
void tnsMakeLinerGradient3dv(real* arr, int num_points, real* rgb0, real* rgb1);
void tnsMakeLinerGradient4dv(real* arr, int num_points, real* rgb0, real* rgb1);

void tnsMakeFoucsSquare(int L, int R, int U, int B, int W);
void tnsDrawFloor(int Size, int Span, int* ShowAxis);
void tnsDrawFloorWithTheme(int Size, int Span, int* ShowAxis, nThemeState* ts);

void tnsMakeIndexUShort(unsigned short* result, short num, ...);
void tnsMakeBridgedIndex(unsigned short* result, short num, int revert, short begin);

//render================================================================


tnsRenderBuffer* tnsCreateRenderBuffer();
void tnsTransformGeomtryToRenderBuffer(tnsRenderBuffer* rb, tnsScene* Scene, tnsCamera* Camera);
void tnsCreateRasterizerbuffer(tnsRenderBuffer* rb, int W, int H, int Samples, int TileW, int TileH);
void tnsMakeRenderGeometryBuffers(tnsScene* s, tnsCamera* c, tnsRenderBuffer* rb);



void tns_RegisterResourcesForSoftwareRender();


//================================

void DrawWireRect2dp(real x, real y, real x2, real y2);
void tnsViewportWithScissor(int x, int y, int w, int h);



void tnsSingleLinearToLog(real* a, real gamma);
void tnsSingleLogToLinear(real* a, real gamma);
void tnsLinearToLog(real* rgb, real gamma);
void tnsLogToLinear(real* rgb, real gamma);
void tnsRgbToLuminance(real* rgb);
void tnsHSVtoRGB(real* hsv, real* rgb);
void tnsRGBtoHSV(real* rgb, real* hsv);




//curve/surf============================================================

tnsBezierHandle* tnsCreateBezierHandle(
	real PX, real PY, real PZ,
	real PLX, real PLY, real PLZ,
	real PRX, real PRY, real PRZ,
	real Thickness, real Tilt,
	u8bit LMode, u8bit RMode);
tnsBezierCurve* tnsCreateBezierCurve();
void tnsAppendBezierSegment(tnsBezierCurve* Target, tnsBezierHandle* Handle);
real tnsLinearInterpolate(real L, real R, real T);
void tnsLinearInterpolate2dv(real* L, real* R, real T, real* Result);
void tnsLinearInterpolate3dv(real* L, real* R, real T, real* Result);
void tnsLinearInterpolatePerspective4dv(tnsVector4d LG, tnsVector4d RG, tnsVector4d L, tnsVector4d R, real T, tnsVector3d Result);
void tnsLinearInterpolate3dv(real* L, real* R, real T, real* Result);
int tnsGetBezierPointCount(tnsBezierCurve* From, int StepsEach);
void tnsGetBezierPoint3dLR(tnsBezierHandle* L, tnsBezierHandle* R, real t, real* Result);
void tnsGetBezierPoint2dLRV(tnsVector2d L, tnsVector2d LH, tnsVector2d RH, tnsVector2d R, real t, real* Result);
void tnsGetBezierPoint3d(tnsBezierHandle* L, real t, real* Result);
real tnsGetBezierPointX(tnsBezierHandle* L, real t);
real tnsGetBezierPointY(tnsBezierHandle* L, real t);
real tnsGetBezierPointYExt(tnsBezierHandle* L, tnsBezierHandle* R, real t);
void tnsGetBezierSubHandle3d(tnsBezierHandle* L, real t, tnsBezierHandle* Result);
void tnsGetBezierSegmentPointArray3d(tnsBezierHandle* L, int StepsEach, int WithTail, real* Result);
void tnsGetBezierPointArray3d(tnsBezierCurve* From, int StepsEach, real* Result);
void tnsGetBeizierPointAtX(tnsBezierHandle* L, real X, real* Result);
real tnsGetBezierRatioX(tnsBezierHandle* L, real X);
real tnsGetBezierYX(tnsBezierHandle* L, real X);
real tnsGetBezierYXExt(tnsBezierHandle* L, tnsBezierHandle* R, real X);
real tnsGetBezierNodeTagent2dL(tnsBezierHandle* BH);
real tnsGetBezierNodeTagent2dR(tnsBezierHandle* BH);
void tnsSetBezierHandles2dL(tnsBezierHandle* BH, real LengthL, real TangentL);
void tnsSetBezierHandles2dR(tnsBezierHandle* BH, real LengthR, real TangentR);
void tnsSetBezierHandles2d(tnsBezierHandle* BH, real LengthL, real LengthR, real TangentL, real TangentR);
void tnsSetBezierHandles2dDistXL(tnsBezierHandle* BH, real XDistL, real TangentL);
void tnsSetBezierHandles2dDistXR(tnsBezierHandle* BH, real XDistR, real TangentR);
void tnsSetBezierHandles2dDistX(tnsBezierHandle* BH, real XDistL, real XDistR, real TangentL, real TangentR);
void tnsSetBezierHandles2dXL(tnsBezierHandle* BH, real XL, real TangentL);
void tnsSetBezierHandles2dXR(tnsBezierHandle* BH, real XR, real TangentR);
void tnsSetBezierHandles2dX(tnsBezierHandle* BH, real XL, real XR, real TangentL, real TangentR);
void tnsSetBezierHandles2dXYL(tnsBezierHandle* BH, real X, real Y);
void tnsSetBezierHandles2dXYR(tnsBezierHandle* BH, real X, real Y);
void tnsSetBezierHandles2dXY(tnsBezierHandle* BH, real XL, real YL, real XR, real YR);
real tnsGetBezierHandleLength2dL(tnsBezierHandle* BH);
real tnsGetBezierHandleLength2dR(tnsBezierHandle* BH);
void tnsClampBezierHandleLY(tnsBezierHandle* BH, real MaxY, real MinY);
void tnsClampBezierHandleRY(tnsBezierHandle* BH, real MaxY, real MinY);
void tnsClampBezierHandleY(tnsBezierHandle* BH, real MaxY, real MinY);
void tnsMirrorBezierHandleToL(tnsBezierHandle* BH);
void tnsMirrorBezierHandleToR(tnsBezierHandle* BH);
void tnsCopyBezierHandle(tnsBezierHandle* From, tnsBezierHandle* To);
void tnsApplyBezierHandleTransformFrom(tnsBezierHandle* Result, tnsMatrix44f* mat, tnsBezierHandle* From);
void tnsApplyBezierHandleTransform(tnsBezierHandle* Result, tnsMatrix44f* mat);


tnsBezierObject* tnsCreateBezierObject(tnsScene* Scene, char* Name, tnsBezierCurve* bc, real AtX, real AtY, real AtZ);

void tnsDrawBezierObject(tnsBezierObject* bo, int OverrideDisplayMode);




tnsLineStrip* tnsCreateLineStrip();
tnsLineStripPoint* tnsAppendPoint(tnsLineStrip* ls, real X, real Y, real Z);
tnsLineStripPoint* tnsPushPoint(tnsLineStrip* ls, real X, real Y, real Z);
void tnsRemovePoint(tnsLineStrip* ls, tnsLineStripPoint* lsp);
void tnsDestroyLineStrip(tnsLineStrip* ls);



void tnsGetGaussianKernel(tnsFilterKernel* fk, int size, real sigma);

nPropContainer* _TNS_PROP_SCENE;


