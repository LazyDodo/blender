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
from .properties_grease_pencil_common import (
        GreasePencilDataPanel,
        GreasePencilOnionPanel,
        )

###############################
# Base-Classes (for shared stuff - e.g. poll, attributes, etc.)

class DataButtonsPanel:
    bl_space_type = 'PROPERTIES'
    bl_region_type = 'WINDOW'
    bl_context = "data"

    @classmethod
    def poll(cls, context):
        return context.object and context.object.type == 'GPENCIL'


class LayerDataButtonsPanel:
    bl_space_type = 'PROPERTIES'
    bl_region_type = 'WINDOW'
    bl_context = "data"

    @classmethod
    def poll(cls, context):
        return (context.object and
                context.object.type == 'GPENCIL' and
                context.active_gpencil_layer)


###############################
# GP Object Properties Panels and Helper Classes

class DATA_PT_gpencil(DataButtonsPanel, Panel):
    bl_label = ""
    bl_options = {'HIDE_HEADER'}

    def draw(self, context):
        layout = self.layout

        # Grease Pencil data selector
        gpd_owner = context.gpencil_data_owner
        gpd = context.gpencil_data

        layout.template_ID(gpd_owner, "data", new="gpencil.data_add", unlink="gpencil.data_unlink")


class DATA_PT_gpencil_datapanel(GreasePencilDataPanel, Panel):
    bl_space_type = 'PROPERTIES'
    bl_region_type = 'WINDOW'
    bl_context = "data"
    bl_label = "Layers"

    # NOTE: this is just a wrapper around the generic GP Panel


class DATA_PT_gpencil_layer_optionpanel(LayerDataButtonsPanel, Panel):
    bl_space_type = 'PROPERTIES'
    bl_region_type = 'WINDOW'
    bl_context = "data"
    bl_label = "Layer Adjustments"
    bl_options = {'DEFAULT_CLOSED'}

    def draw(self, context):
        layout = self.layout
        layout.use_property_split = True

        gpl = context.active_gpencil_layer
        layout.active = not gpl.lock

        # Layer options
        # Offsets - Color Tint
        layout.enabled = not gpl.lock
        col = layout.column(align=True)
        col.prop(gpl, "tint_color")
        col.prop(gpl, "tint_factor", slider=True)

        # Offsets - Thickness
        row = layout.row(align=True)
        row.prop(gpl, "line_change", text="Thickness")
        row.operator("gpencil.stroke_apply_thickness", icon='STYLUS_PRESSURE', text="")

        layout.prop(gpl, "use_stroke_location", text="Draw On Stroke Location")


class DATA_PT_gpencil_onionpanel(Panel):
    bl_space_type = 'PROPERTIES'
    bl_region_type = 'WINDOW'
    bl_context = "data"
    bl_label = "Onion Skinning"
    bl_options = {'DEFAULT_CLOSED'}

    @classmethod
    def poll(cls, context):
        return bool(context.active_gpencil_layer)

    @staticmethod
    def draw_header(self, context):
        self.layout.prop(context.gpencil_data, "use_onion_skinning", text="")

    def draw(self, context):
        gpd = context.gpencil_data

        layout = self.layout
        layout.enabled = gpd.use_onion_skinning

        GreasePencilOnionPanel.draw_settings(layout, gpd)


class DATA_PT_gpencil_layer_onionpanel(Panel):
    bl_space_type = 'PROPERTIES'
    bl_region_type = 'WINDOW'
    bl_context = "data"
    bl_label = "Onion Skinning (Layer Override)"
    bl_options = {'DEFAULT_CLOSED'}

    @classmethod
    def poll(cls, context):
        return bool(context.active_gpencil_layer)

    def draw_header(self, context):
        self.layout.prop(context.active_gpencil_layer, "override_onion", text="")

    def draw(self, context):
        gpd = context.gpencil_data
        gpl = context.active_gpencil_layer

        layout = self.layout
        layout.enabled = gpd.use_onion_skinning and gpl.override_onion

        GreasePencilOnionPanel.draw_settings(layout, gpl)


class DATA_PT_gpencil_parentpanel(LayerDataButtonsPanel, Panel):
    bl_space_type = 'PROPERTIES'
    bl_region_type = 'WINDOW'
    bl_context = "data"
    bl_label = "Layer Relations"
    bl_options = {'DEFAULT_CLOSED'}

    def draw(self, context):
        layout = self.layout
        scene = context.scene
        gpl = context.active_gpencil_layer
        row = layout.row()

        col = row.column(align=True)
        col.active = not gpl.lock
        col.label(text="Parent:")
        col.prop(gpl, "parent", text="")

        sub = col.column()
        sub.prop(gpl, "parent_type", text="")
        parent = gpl.parent
        if parent and gpl.parent_type == 'BONE' and parent.type == 'ARMATURE':
            sub.prop_search(gpl, "parent_bone", parent.data, "bones", text="")


class GPENCIL_UL_vgroups(UIList):
    def draw_item(self, context, layout, data, item, icon, active_data, active_propname, index):
        vgroup = item
        if self.layout_type in {'DEFAULT', 'COMPACT'}:
            layout.prop(vgroup, "name", text="", emboss=False, icon_value=icon)
            # icon = 'LOCKED' if vgroup.lock_weight else 'UNLOCKED'
            # layout.prop(vgroup, "lock_weight", text="", icon=icon, emboss=False)
        elif self.layout_type == 'GRID':
            layout.alignment = 'CENTER'
            layout.label(text="", icon_value=icon)


class DATA_PT_gpencil_vertexpanel(DataButtonsPanel, Panel):
    bl_space_type = 'PROPERTIES'
    bl_region_type = 'WINDOW'
    bl_context = "data"
    bl_label = "Vertex Groups"
    bl_options = {'DEFAULT_CLOSED'}

    def draw(self, context):
        layout = self.layout

        ob = context.object
        group = ob.vertex_groups.active

        rows = 2
        if group:
            rows = 4

        row = layout.row()
        row.template_list("GPENCIL_UL_vgroups", "", ob, "vertex_groups", ob.vertex_groups, "active_index", rows=rows)

        col = row.column(align=True)
        col.operator("object.vertex_group_add", icon='ZOOMIN', text="")
        col.operator("object.vertex_group_remove", icon='ZOOMOUT', text="").all = False

        if ob.vertex_groups:
            row = layout.row()

            sub = row.row(align=True)
            sub.operator("gpencil.vertex_group_assign", text="Assign")
            sub.operator("gpencil.vertex_group_remove_from", text="Remove")

            sub = row.row(align=True)
            sub.operator("gpencil.vertex_group_select", text="Select")
            sub.operator("gpencil.vertex_group_deselect", text="Deselect")

            layout.prop(context.tool_settings, "vertex_group_weight", text="Weight")


class DATA_PT_gpencil_infopanel(DataButtonsPanel, Panel):
    bl_space_type = 'PROPERTIES'
    bl_region_type = 'WINDOW'
    bl_context = "data"
    bl_label = "Information"
    bl_options = {'DEFAULT_CLOSED'}

    def draw(self, context):
        layout = self.layout
        gpd = context.gpencil_data

        split = layout.split(percentage=0.5)

        col = split.column(align=True)
        col.label("Layers:", icon="LAYER_ACTIVE")
        col.label("Frames:", icon="LAYER_ACTIVE")
        col.label("Strokes:", icon="LAYER_ACTIVE")
        col.label("Points:", icon="LAYER_ACTIVE")

        col = split.column(align=True)
        col.label(str(gpd.info_total_layers))
        col.label(str(gpd.info_total_frames))
        col.label(str(gpd.info_total_strokes))
        col.label(str(gpd.info_total_points))


class DATA_PT_gpencil_display(DataButtonsPanel, Panel):
    bl_label = "Viewport Display"
    bl_options = {'DEFAULT_CLOSED'}

    def draw(self, context):
        layout = self.layout

        ob = context.object

        gpd = context.gpencil_data
        gpl = context.active_gpencil_layer


        layout.prop(gpd, "xray_mode", text="Depth Ordering")
        layout.prop(ob, "empty_draw_size", text="Marker Size")

        layout.separator()

        if gpl:
            col = layout.column(align=True)
            col.prop(gpd, "show_stroke_direction", text="Show Stroke Directions")

        layout.separator()

        col = layout.column(align=True)
        col.prop(gpd, "show_constant_thickness")
        sub = col.column()
        sub.active = not gpd.show_constant_thickness
        sub.prop(gpd, "pixfactor", text="Scale")

        layout.separator()

        col = layout.column()
        col.prop(gpd, "show_edit_lines", text="Show Edit Lines")
        col.prop(gpd, "edit_line_color", text="")
        col.prop(gpd, "show_multiedit_line_only", text="Only Lines in MultiEdit")


###############################

classes = (
    DATA_PT_gpencil,
    DATA_PT_gpencil_datapanel,
    DATA_PT_gpencil_onionpanel,
    DATA_PT_gpencil_layer_onionpanel,
    DATA_PT_gpencil_layer_optionpanel,
    DATA_PT_gpencil_parentpanel,
    DATA_PT_gpencil_vertexpanel,
    DATA_PT_gpencil_display,
    DATA_PT_gpencil_infopanel,

    GPENCIL_UL_vgroups,
)

if __name__ == "__main__":  # only for live edit.
    from bpy.utils import register_class
    for cls in classes:
        register_class(cls)
