/*
 * ***** BEGIN GPL LICENSE BLOCK *****
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 *
 * Copyright 2015, Blender Foundation.
 *
 * ***** END GPL LICENSE BLOCK *****
 */

/** \file blender/python/intern/gpu_draw.c
 *  \ingroup pythonintern
 *
 * This file defines the draw functionalities of the 'gpu' module
 * used for off-screen OpenGL rendering.
 */

#include <Python.h>

#include "BLI_utildefines.h"

#include "DNA_image_types.h"

#include "GPU_immediate.h"
#include "GPU_texture.h"

#include "../generic/py_capi_utils.h"

#include "gpu.h"

/* -------------------------------------------------------------------- */
/* GPU draw methods */

PyDoc_STRVAR(pygpu_draw_rect_doc,
"rect(x1, y1, x2, y2, r, g, b, a)\n"
);
static PyObject *pygpu_draw_rect(PyObject *UNUSED(self), PyObject *args, PyObject *kwds)
{
	static const char *kwlist[] = {"x1", "y1", "x2", "y2", "r", "g", "b", "a", NULL};
	float x1, y1, x2, y2, rgba[4];

	if (!PyArg_ParseTupleAndKeywords(
	        args, kwds, "ffffffff", (char **)(kwlist),
	        &x1, &y1, &x2, &y2, &rgba[0], &rgba[1], &rgba[2], &rgba[3]))
	{
		return NULL;
	}

	GPU_blend(true);
	GPU_blend_set_func_separate(GPU_SRC_ALPHA, GPU_ONE_MINUS_SRC_ALPHA, GPU_ONE, GPU_ONE_MINUS_SRC_ALPHA);

	uint pos = GPU_vertformat_attr_add(immVertexFormat(), "pos", GPU_COMP_F32, 2, GPU_FETCH_FLOAT);
	immBindBuiltinProgram(GPU_SHADER_2D_UNIFORM_COLOR);

	immUniformColor4fv(rgba);
	immRectf(pos, x1, y1, x2, y2);

	immUnbindProgram();

	GPU_blend(false);

	Py_RETURN_NONE;
}

PyDoc_STRVAR(pygpu_draw_image_doc,
"image(image, x1, y1, x2, y2)\n"
);
static PyObject *pygpu_draw_image(PyObject *UNUSED(self), PyObject *args, PyObject *kwds)
{
	static const char *kwlist[] = {"image", "x1", "y1", "x2", "y2", NULL};
	PyObject *py_image;
	Image *image;
	float x1, y1, x2, y2;

	if (!PyArg_ParseTupleAndKeywords(
	        args, kwds, "Offff", (char **)(kwlist),
	        &py_image, &x1, &y1, &x2, &y2))

	{
		return NULL;
	}

	if (!(image = PyC_RNA_AsPointer(py_image, "Image"))) {
		return NULL;
	}

	GPUTexture *tex = GPU_texture_from_blender(image, NULL, GL_TEXTURE_2D, false, 0.0f);

	if (tex == NULL) {
		return NULL;
	}

	GPU_texture_bind(tex, 0);

	GPUVertFormat *format = immVertexFormat();
	uint texcoord = GPU_vertformat_attr_add(format, "texCoord", GPU_COMP_F32, 2, GPU_FETCH_FLOAT);
	uint pos = GPU_vertformat_attr_add(format, "pos", GPU_COMP_F32, 2, GPU_FETCH_FLOAT);
	immBindBuiltinProgram(GPU_SHADER_2D_IMAGE);

	immUniform1i("image", 0);

	immBegin(GPU_PRIM_TRI_FAN, 4);

	immAttrib2f(texcoord, 0.0f, 0.0f);
	immVertex2f(pos, x1, y1);

	immAttrib2f(texcoord, 1.0f, 0.0f);
	immVertex2f(pos, x2, y1);

	immAttrib2f(texcoord, 1.0f, 1.0f);
	immVertex2f(pos, x2, y2);

	immAttrib2f(texcoord, 0.0f, 1.0f);
	immVertex2f(pos, x1, y2);

	immEnd();

	immUnbindProgram();

	GPU_texture_unbind(tex);

	Py_RETURN_NONE;
}

static struct PyMethodDef BPy_GPU_draw_methods[] = {
	{"rect", (PyCFunction)pygpu_draw_rect, METH_VARARGS | METH_KEYWORDS, pygpu_draw_rect_doc},
	{"image", (PyCFunction)pygpu_draw_image, METH_VARARGS | METH_KEYWORDS, pygpu_draw_image_doc},
	{NULL, NULL, 0, NULL}
};

PyDoc_STRVAR(BPy_GPU_draw_doc,
"This module provides access to drawing functions."
);
static PyModuleDef BPy_GPU_draw_module_def = {
	PyModuleDef_HEAD_INIT,
	"gpu.draw",                                  /* m_name */
	BPy_GPU_draw_doc,                            /* m_doc */
	0,                                           /* m_size */
	BPy_GPU_draw_methods,                        /* m_methods */
	NULL,                                        /* m_reload */
	NULL,                                        /* m_traverse */
	NULL,                                        /* m_clear */
	NULL,                                        /* m_free */
};

PyObject *BPyInit_gpu_draw(void)
{
	PyObject *submodule;

	submodule = PyModule_Create(&BPy_GPU_draw_module_def);

	return submodule;
}

#undef BPY_GPU_OFFSCREEN_CHECK_OBJ
