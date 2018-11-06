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
 * Contributor(s): Chingiz Dyussenov, Arystanbek Dyussenov, Nathan Letwory.
 *
 * ***** END GPL LICENSE BLOCK *****
 */

#ifndef __MATERIAL_H__
#define __MATERIAL_H__

#include <map>

extern "C" {
#include "BKE_context.h"
#include "BKE_node.h"
#include "BLI_listbase.h"
#include "DNA_material_types.h"
#include "DNA_node_types.h"
}

#include "COLLADAFWEffectCommon.h"

typedef enum BC_pbr_inputs {
	BC_PBR_DIFFUSE = 0,
	BC_PBR_METALLIC = 4,
	BC_PBR_IOR = 14
} BC_pbr_inputs;

typedef std::map<COLLADAFW::UniqueId, Image*> Image_map;

class MaterialNode {

private:
	bContext *mContext;
	Material *material;
	COLLADAFW::EffectCommon *effect;
	Image_map &uid_image_map;
	bNodeTree *ntree;

	bNode *shader_node;
	bNode *output_node;

	bNodeTree *prepare_material_nodetree();
	bNode *bc_add_node(int node_type, int locx, int locy, std::string label);
	void bc_node_add_link(bNode *from_node, int from_index, bNode *to_node, int to_index);
	void setShaderType();

public:
	MaterialNode(bContext *C, COLLADAFW::EffectCommon *ef, Material *ma, Image_map &uid_image_map);
	void set_diffuse(COLLADAFW::ColorOrTexture &cot);
	void set_reflectivity(float val);
	void set_ior(float val);

};

#endif
