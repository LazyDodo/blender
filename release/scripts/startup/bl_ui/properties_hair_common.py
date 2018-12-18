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

# <pep8-80 compliant>

import bpy


def draw_hair_display_settings(layout, settings):
    col = layout.column(align=True)
    col.label("Follicles:")
    col.prop(settings, "follicle_mode", expand=True)

    col = layout.column(align=True)
    col.label("Guide Curves:")
    col.prop(settings, "guide_mode", expand=True)

    layout.prop(settings, "shape")

    col = layout.column(align=True)
    col.prop(settings, "root_radius")
    col.prop(settings, "tip_radius")

    col = layout.column()
    col.prop(settings, "radius_scale")
    col.prop(settings, "use_close_tip")


class HAIR_PT_display_settings:
    # subclasses must define...
    # ~ bl_space_type = 'PROPERTIES'
    # ~ bl_region_type = 'WINDOW'
    bl_label = "Hair Display Settings"

    def draw(self, context):
        settings = context.draw_hair_display_settings
        draw_hair_display_settings(self.layout, hair_display_settings)


classes = (
)

if __name__ == "__main__":  # only for live edit.
    from bpy.utils import register_class
    for cls in classes:
        register_class(cls)