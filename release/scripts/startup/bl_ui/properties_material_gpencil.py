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


class GPENCIL_UL_matslots(UIList):
    def draw_item(self, context, layout, data, item, icon, active_data, active_propname, index):
        slot = item
        ma = slot.material
        if (ma is not None) and (ma.grease_pencil is not None):
            gpcolor = ma.grease_pencil

            if self.layout_type in {'DEFAULT', 'COMPACT'}:
                if gpcolor.lock:
                    layout.active = False

                split = layout.split(percentage=0.25)
                row = split.row(align=True)
                row.enabled = not gpcolor.lock
                row.prop(gpcolor, "color", text="", emboss=gpcolor.is_stroke_visible)
                row.prop(gpcolor, "fill_color", text="", emboss=gpcolor.is_fill_visible)
                split.prop(ma, "name", text="", emboss=False)

                row = layout.row(align=True)
                row.prop(gpcolor, "lock", text="", emboss=False)
                row.prop(gpcolor, "hide", text="", emboss=False)
                if gpcolor.ghost is True:
                    icon = 'GHOST_DISABLED'
                else:
                    icon = 'GHOST_ENABLED'
                row.prop(gpcolor, "ghost", text="", icon=icon, emboss=False)

            elif self.layout_type == 'GRID':
                layout.alignment = 'CENTER'
                layout.label(text="", icon_value=icon)


class GPMaterialButtonsPanel:
    bl_space_type = 'PROPERTIES'
    bl_region_type = 'WINDOW'
    bl_context = "material"

    @classmethod
    def poll(cls, context):
        ob = context.object
        return (ob and ob.type == 'GPENCIL' and
                ob.active_material and
                ob.active_material.grease_pencil)



class MATERIAL_PT_gpencil_slots(Panel):
    bl_label = "Grease Pencil Material Slots"
    bl_space_type = 'PROPERTIES'
    bl_region_type = 'WINDOW'
    bl_context = "material"
    bl_options = {'HIDE_HEADER'}

    @classmethod
    def poll(cls, context):
        ob = context.object
        return ob and ob.type == 'GPENCIL'

    @staticmethod
    def draw(self, context):
        layout = self.layout
        gpd = context.gpencil_data

        mat = context.object.active_material
        ob = context.object
        slot = context.material_slot
        space = context.space_data

        if ob:
            is_sortable = len(ob.material_slots) > 1
            rows = 1
            if (is_sortable):
                rows = 4

            row = layout.row()

            row.template_list("GPENCIL_UL_matslots", "", ob, "material_slots", ob, "active_material_index", rows=rows)

            col = row.column(align=True)
            col.operator("object.material_slot_add", icon='ZOOMIN', text="")
            col.operator("object.material_slot_remove", icon='ZOOMOUT', text="")

            col.menu("GPENCIL_MT_color_specials", icon='DOWNARROW_HLT', text="")

            if is_sortable:
                col.separator()

                col.operator("object.material_slot_move", icon='TRIA_UP', text="").direction = 'UP'
                col.operator("object.material_slot_move", icon='TRIA_DOWN', text="").direction = 'DOWN'

                col.separator()

                sub = col.column(align=True)
                sub.operator("gpencil.color_isolate", icon='LOCKED', text="").affect_visibility = False
                sub.operator("gpencil.color_isolate", icon='RESTRICT_VIEW_OFF', text="").affect_visibility = True

            if gpd.use_stroke_edit_mode:
                row = layout.row(align=True)
                row.operator("gpencil.stroke_change_color", text="Assign")
                row.operator("gpencil.color_select", text="Select")
                #row.operator("gpencil.color_deselect", text="Deselect")

        split = layout.split(percentage=0.65)

        if ob:
            split.template_ID(ob, "active_material", new="material.new")
            row = split.row()

            if slot:
                row.prop(slot, "link", text="")
            else:
                row.label()
        elif mat:
            split.template_ID(space, "pin_id")
            split.separator()


class MATERIAL_PT_gpencil_strokecolor(GPMaterialButtonsPanel, Panel):
    bl_label = "Stroke"

    @staticmethod
    def draw(self, context):
        layout = self.layout
        layout.use_property_split = True

        ma = context.object.active_material
        if (ma is not None) and (ma.grease_pencil):
            gpcolor = ma.grease_pencil

            split = layout.split(percentage=1.0)
            split.active = not gpcolor.lock

            col = split.column(align=True)
            row = col.row(align=True)
            row.enabled = not gpcolor.lock
            row.prop(gpcolor, "mode")
            col.separator()

            col.enabled = not gpcolor.lock
            col.prop(gpcolor, "stroke_style", text="Style")

            if gpcolor.stroke_style == 'TEXTURE':
                row = layout.row()
                row.enabled = not gpcolor.lock
                col = row.column(align=True)
                col.template_ID(gpcolor, "stroke_image", open="image.open")
                col.prop(gpcolor, "pixel_size", text="UV Factor")
                col.prop(gpcolor, "use_pattern", text="Use As Pattern")

            if gpcolor.stroke_style == 'SOLID' or gpcolor.use_pattern is True:
                row = layout.row()
                col = row.column(align=True)
                col.prop(gpcolor, "color", text="Color")
                col.prop(gpcolor, "alpha", slider=True)

            # Options
            row = layout.row()
            row.active = not gpcolor.lock
            col = row.column(align=True)
            col.prop(gpcolor, "pass_index")


class MATERIAL_PT_gpencil_fillcolor(GPMaterialButtonsPanel, Panel):
    bl_label = "Fill"

    @staticmethod
    def draw(self, context):
        layout = self.layout
        layout.use_property_split = True

        ma = context.object.active_material
        if (ma is not None) and (ma.grease_pencil):
            gpcolor = ma.grease_pencil

            # color settings
            split = layout.split(percentage=1.0)
            split.active = not gpcolor.lock

            row = layout.row()
            col = row.column(align=True)
            col.enabled = not gpcolor.lock
            col.prop(gpcolor, "fill_style", text="Style")

            row = layout.row()
            col = row.column(align=True)

            if gpcolor.fill_style != 'TEXTURE':
                col.prop(gpcolor, "fill_color", text="Color")
                col.prop(gpcolor, "fill_alpha", text="Opacity", slider=True)
                col.separator()
                if gpcolor.texture_mix is True or gpcolor.fill_style in ('GRADIENT', 'RADIAL'):
                    col.prop(gpcolor, "mix_factor", text="Mix", slider=True)

            if gpcolor.fill_style in ('GRADIENT', 'RADIAL', 'CHESSBOARD'):
                if gpcolor.texture_mix is False or gpcolor.fill_style == 'CHESSBOARD':
                    col.prop(gpcolor, "mix_color", text="Mix Color")
                split = col.split(percentage=0.5)
                subcol = split.column(align=True)
                subcol.prop(gpcolor, "pattern_shift", text="Location")
                subrow = subcol.row(align=True)
                if gpcolor.fill_style == 'RADIAL':
                    subrow.enabled = False
                subrow.prop(gpcolor, "pattern_angle", text="Angle")
                subcol.prop(gpcolor, "flip", text="Flip")

                subcol = split.column(align=True)
                subcol.prop(gpcolor, "pattern_scale", text="Scale")
                subrow = subcol.row(align=True)
                if gpcolor.fill_style != 'RADIAL':
                    subrow.enabled = False
                subrow.prop(gpcolor, "pattern_radius", text="Radius")
                subrow = subcol.row(align=True)
                if gpcolor.fill_style != 'CHESSBOARD':
                    subrow.enabled = False
                subrow.prop(gpcolor, "pattern_boxsize", text="Box")

            col.separator()
            col.label("Texture")
            if gpcolor.fill_style not in ('TEXTURE', 'PATTERN'):
                col.prop(gpcolor, "texture_mix", text="Mix Texture")
            if gpcolor.fill_style in ('TEXTURE', 'PATTERN') or gpcolor.texture_mix is True:
                col.template_ID(gpcolor, "fill_image", open="image.open")
                split = col.split(percentage=0.5)
                subcol = split.column(align=True)
                subcol.prop(gpcolor, "texture_offset", text="Offset")
                subcol.prop(gpcolor, "texture_angle")
                subcol.prop(gpcolor, "texture_clamp", text="Clip Image")
                subcol = split.column(align=True)
                subcol.prop(gpcolor, "texture_scale", text="Scale")
                subcol.prop(gpcolor, "texture_opacity")


classes = (
    GPENCIL_UL_matslots,
    MATERIAL_PT_gpencil_slots,
    MATERIAL_PT_gpencil_strokecolor,
    MATERIAL_PT_gpencil_fillcolor,
)

if __name__ == "__main__":  # only for live edit.
    from bpy.utils import register_class
    for cls in classes:
        register_class(cls)
