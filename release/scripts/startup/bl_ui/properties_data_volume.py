# ##### BEGIN GPL LICENSE BLOCK #####
#
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 2
#  of the License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software Foundation,
#  Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
#
# ##### END GPL LICENSE BLOCK #####

# <pep8 compliant>
import bpy
from bpy.types import Panel, UIList
from rna_prop_ui import PropertyPanel


class DataButtonsPanel:
    bl_space_type = 'PROPERTIES'
    bl_region_type = 'WINDOW'
    bl_context = "data"

    @classmethod
    def poll(cls, context):
        engine = context.scene.render.engine
        return context.volume and (engine in cls.COMPAT_ENGINES)


class DATA_PT_context_volume(DataButtonsPanel, Panel):
    bl_label = ""
    bl_options = {'HIDE_HEADER'}
    COMPAT_ENGINES = {'BLENDER_RENDER', 'BLENDER_GAME'}

    def draw(self, context):
        layout = self.layout

        ob = context.object
        volume = context.volume
        space = context.space_data

        split = layout.split(percentage=0.65)

        if ob:
            split.template_ID(ob, "data")
        elif volume:
            split.template_ID(space, "pin_id")


class DATA_PT_volume(DataButtonsPanel, Panel):
    bl_label = "Volume"
    COMPAT_ENGINES = {'BLENDER_RENDER'}

    def draw(self, context):
        layout = self.layout

        volume = context.volume

        layout.prop(volume, "filepath", text="")


class VOLUME_UL_grids(UIList):
    def draw_item(self, context, layout, data, item, icon, active_data, active_propname, index):
        layout.prop(item, "name", text="", emboss=False, icon='OUTLINER_DATA_LATTICE')


class DATA_PT_volume_grids(DataButtonsPanel, Panel):
    bl_label = "Grids"
    COMPAT_ENGINES = {'BLENDER_RENDER'}

    def draw(self, context):
        layout = self.layout

        volume = context.volume

        layout.template_list("VOLUME_UL_grids", "grids", volume, "grids", volume.grids, "active_index", rows=2)

class DATA_PT_custom_props_volume(DataButtonsPanel, PropertyPanel, Panel):
    COMPAT_ENGINES = {'BLENDER_RENDER', 'BLENDER_GAME'}
    _context_path = "object.data"
    _property_type = bpy.types.Volume


classes = (
    DATA_PT_context_volume,
    DATA_PT_volume,
    DATA_PT_volume_grids,
    DATA_PT_custom_props_volume,
    VOLUME_UL_grids,
)

if __name__ == "__main__":  # only for live edit.
    from bpy.utils import register_class
    for cls in classes:
        register_class(cls)
