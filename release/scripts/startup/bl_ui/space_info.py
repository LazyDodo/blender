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
from bpy.types import Header, Menu


class INFO_HT_header(Header):
    bl_space_type = 'INFO'

    def draw(self, context):
        layout = self.layout
        layout.template_header()

        # Empty for now until info editor gets turned into log editor
        pass

"""
        if (bpy.app.assets_fail or bpy.app.assets_need_reload) and not bpy.app.assets_quiet:
            row.operator("script.assets_warn_clear", text="Ignore")
            if bpy.app.assets_need_reload is True and bpy.app.assets_quiet is False:
                row.operator("wm.assets_reload", icon='SCREEN_BACK', text="Reload Assets")
                row.label("Some assets have to be reloaded", icon='INFO')
            if bpy.app.assets_fail is True and bpy.app.assets_quiet is False:
                row.label("Some asset engine(s) failed to retrieve updated data about their assets...", icon='ERROR')
            return
"""


# Not really info, just add to re-usable location.
class INFO_MT_area(Menu):
    bl_label = "Area"

    def draw(self, context):
        layout = self.layout

        layout.operator("screen.area_dupli")
        if context.space_data.type == 'VIEW_3D':
            layout.operator("screen.region_quadview")
        layout.operator("screen.screen_full_area")
        layout.operator("screen.screen_full_area", text="Toggle Fullscreen Area").use_hide_panels = True


classes = (
    INFO_HT_header,
    INFO_MT_area,
)

if __name__ == "__main__":  # only for live edit.
    from bpy.utils import register_class
    for cls in classes:
        register_class(cls)
