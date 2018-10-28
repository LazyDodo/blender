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
from bpy.types import Menu, Panel
from rna_prop_ui import PropertyPanel


class DataButtonsPanel:
    bl_space_type = 'PROPERTIES'
    bl_region_type = 'WINDOW'
    bl_context = "data"

    @classmethod
    def poll(cls, context):
        return context.hair_system


class DATA_PT_context_hair(DataButtonsPanel, Panel):
    bl_label = ""
    bl_options = {'HIDE_HEADER'}

    def draw(self, context):
        layout = self.layout
        ob = context.object
        hsys = context.hair_system
        space = context.space_data

        split = layout.split(factor=0.65)

        if ob:
            split.template_ID(ob, "data")
        elif hsys:
            split.template_ID(space, "pin_id")


class DATA_PT_hair_system(DataButtonsPanel, Panel):
    bl_label = "Hair System"

    def draw(self, context):
        layout = self.layout
        hsys = context.hair_system

        split = layout.split()

        #col = split.column()
        #col.alert = (hsys.scalp_object is None)
        #col.label("Scalp Object:")
        #col.prop(hsys, "scalp_object", "")


class DATA_PT_hair_draw_settings(DataButtonsPanel, Panel):
    bl_label = "Draw Settings"

    def draw(self, context):
        layout = self.layout
        hsys = context.hair_system
        ds = hsys.draw_settings

        layout.prop(hsys, "material_slot")

        col = layout.column(align=True)
        col.label("Follicles:")
        col.prop(ds, "follicle_mode", expand=True)

        col = layout.column(align=True)
        col.label("Guide Curves:")
        col.prop(ds, "guide_mode", expand=True)

        layout.prop(ds, "shape")

        col = layout.column(align=True)
        col.prop(ds, "root_radius")
        col.prop(ds, "tip_radius")

        col = layout.column()
        col.prop(ds, "radius_scale")
        col.prop(ds, "use_close_tip")


class DATA_PT_custom_props_hair(DataButtonsPanel, PropertyPanel, Panel):
    _context_path = "object.data"
    _property_type = bpy.types.HairSystem


classes = (
    DATA_PT_context_hair,
    DATA_PT_hair_system,
    DATA_PT_hair_draw_settings,
    DATA_PT_custom_props_hair,
)

if __name__ == "__main__":  # only for live edit.
    from bpy.utils import register_class
    for cls in classes:
        register_class(cls)
