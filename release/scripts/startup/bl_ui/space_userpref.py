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
from bpy.types import (
    Header,
    Menu,
    Panel,
)
from bpy.app.translations import pgettext_iface as iface_
from bpy.app.translations import contexts as i18n_contexts


class USERPREF_HT_header(Header):
    bl_space_type = 'USER_PREFERENCES'

    def draw(self, context):
        layout = self.layout

        layout.template_header()

        userpref = context.user_preferences


        if userpref.active_section == 'LIGHTS':
            layout.operator('wm.studiolight_install', text="Add MatCap").type = 'MATCAP'
            layout.operator('wm.studiolight_install', text="Add LookDev HDRI").type = 'WORLD'
            op = layout.operator('wm.studiolight_install', text="Add Studio Light")
            op.type = 'STUDIO'
            op.filter_glob = ".sl"


class USERPREF_PT_navigation(Panel):
    bl_label = ""
    bl_space_type = 'USER_PREFERENCES'
    bl_region_type = 'NAVIGATION_BAR'
    bl_options = {'HIDE_HEADER'}

    def draw(self, context):
        layout = self.layout
        layout.operator_context = 'EXEC_AREA'

        userpref = context.user_preferences

        col = layout.column()

        col.scale_x = 1.3
        col.scale_y = 1.3
        col.prop(userpref, "active_section", expand=True)

        col.separator_spacer()
        #Doesn't work vertically yet

        col.operator("wm.save_userpref")


class PreferencePanel(Panel):
    bl_space_type = 'USER_PREFERENCES'
    bl_region_type = 'WINDOW'

    def draw_props(self, context, layout):
        # Deriving classes should implement this.
        # TODO use abc module for abstract method support?
        pass

    def draw(self, context):
        layout = self.layout

        layout.use_property_split = True
        layout.use_property_decorate = False  # No animation.

        row = layout.row()
        row.label()

        col = row.column()
        col.ui_units_x = 50

        self.draw_props(context, col)

        row.label() # Needed so col above is centered.


class USERPREF_PT_interface_display(PreferencePanel):
    bl_label = "Display"

    @classmethod
    def poll(cls, context):
        userpref = context.user_preferences
        return (userpref.active_section == 'INTERFACE')

    def draw_props(self, context, layout):
        userpref = context.user_preferences
        view = userpref.view

        layout.prop(view, "ui_scale", text="Scale")
        layout.prop(view, "ui_line_width", text="Line Width")


class USERPREF_PT_interface_display_info(PreferencePanel):
    bl_label = "Information"
    bl_parent_id = "USERPREF_PT_interface_display"
    bl_options = {'DEFAULT_CLOSED'}

    def draw_props(self, context, layout):
        userpref = context.user_preferences
        view = userpref.view

        layout.prop(view, "show_tooltips")
        layout.prop(view, "show_object_info", text="Object Info")
        layout.prop(view, "show_large_cursors")
        layout.prop(view, "show_view_name", text="View Name")
        layout.prop(view, "show_playback_fps", text="Playback FPS")

class USERPREF_PT_interface_develop(PreferencePanel):
    bl_label = "Develop"
    bl_options = {'DEFAULT_CLOSED'}

    @classmethod
    def poll(cls, context):
        userpref = context.user_preferences
        return (userpref.active_section == 'INTERFACE')

    def draw_props(self, context, layout):
        userpref = context.user_preferences
        view = userpref.view

        layout.prop(view, "show_tooltips_python")
        layout.prop(view, "show_developer_ui")


class USERPREF_PT_interface_view_manipulation(PreferencePanel):
    bl_label = "View Manipulation"
    bl_options = {'DEFAULT_CLOSED'}

    @classmethod
    def poll(cls, context):
        userpref = context.user_preferences
        return (userpref.active_section == 'INTERFACE')

    def draw_props(self, context, layout):
        userpref = context.user_preferences
        view = userpref.view

        layout.prop(view, "smooth_view")
        layout.prop(view, "rotation_angle")

        layout.separator()

        layout.prop(view, "use_mouse_depth_cursor")
        layout.prop(view, "use_cursor_lock_adjust")

        layout.separator()

        layout.prop(view, "use_auto_perspective")
        layout.prop(view, "use_mouse_depth_navigate")

        layout.separator()

        layout.prop(view, "use_zoom_to_mouse")
        layout.prop(view, "use_rotate_around_active")

        layout.separator()

        layout.prop(view, "use_camera_lock_parent")


class USERPREF_PT_interface_viewports(PreferencePanel):
    bl_label = "Viewports"
    bl_options = {'DEFAULT_CLOSED'}

    @classmethod
    def poll(cls, context):
        userpref = context.user_preferences
        return (userpref.active_section == 'INTERFACE')


class USERPREF_PT_interface_viewports_3d(PreferencePanel):
    bl_label = "3D Viewports"
    bl_parent_id = "USERPREF_PT_interface_viewports"

    def draw_props(self, context, layout):
        userpref = context.user_preferences
        view = userpref.view

        layout.prop(view, "object_origin_size")

        layout.separator()

        layout.prop(view, "mini_axis_type", text="3D Viewport Axis")

        sub = layout.column()
        sub.active = view.mini_axis_type == 'MINIMAL'
        sub.prop(view, "mini_axis_size", text="Size")
        sub.prop(view, "mini_axis_brightness", text="Brightness")

        layout.separator()

        layout.prop(view, "gizmo_size", text="Gizmo Size")


class USERPREF_PT_interface_viewports_3d_weight_paint(PreferencePanel):
    bl_label = "Custom Weight Paint Range"
    bl_options = {'DEFAULT_CLOSED'}
    bl_parent_id = "USERPREF_PT_interface_viewports_3d"

    def draw_header(self, context):
        userpref = context.user_preferences
        system = userpref.system

        self.layout.prop(system, "use_weight_color_range", text="")

    def draw_props(self, context, layout):
        userpref = context.user_preferences
        system = userpref.system

        layout.active = system.use_weight_color_range
        layout.template_color_ramp(system, "weight_color_range", expand=True)


class USERPREF_PT_interface_viewports_2d(PreferencePanel):
    bl_label = "2D Viewports"
    bl_parent_id = "USERPREF_PT_interface_viewports"

    def draw_props(self, context, layout):
        userpref = context.user_preferences
        view = userpref.view

        layout.prop(view, "view2d_grid_spacing_min", text="Minimum Grid Spacing")
        layout.prop(view, "timecode_style")
        layout.prop(view, "view_frame_type")
        if view.view_frame_type == 'SECONDS':
            layout.prop(view, "view_frame_seconds")
        elif view.view_frame_type == 'KEYFRAMES':
            layout.prop(view, "view_frame_keyframes")


class USERPREF_PT_interface_menus(PreferencePanel):
    bl_label = "Menus"
    bl_options = {'DEFAULT_CLOSED'}

    @classmethod
    def poll(cls, context):
        userpref = context.user_preferences
        return (userpref.active_section == 'INTERFACE')

    def draw_props(self, context, layout):
        userpref = context.user_preferences
        view = userpref.view

        layout.prop(view, "show_splash")
        layout.prop(view, "use_quit_dialog")


class USERPREF_PT_interface_menus_mouse_over(PreferencePanel):
    bl_label = "Open on Mouse Over"
    bl_parent_id = "USERPREF_PT_interface_menus"
    bl_options = {'DEFAULT_CLOSED'}

    def draw_header(self, context):
        userpref = context.user_preferences
        view = userpref.view

        self.layout.prop(view, "use_mouse_over_open", text="")

    def draw_props(self, context, layout):
        userpref = context.user_preferences
        view = userpref.view

        layout.active = view.use_mouse_over_open

        layout.prop(view, "open_toplevel_delay", text="Top Level")
        layout.prop(view, "open_sublevel_delay", text="Sub Level")

class USERPREF_PT_interface_menus_pie(PreferencePanel):
    bl_label = "Pie Menus"
    bl_parent_id = "USERPREF_PT_interface_menus"
    bl_options = {'DEFAULT_CLOSED'}

    def draw_props(self, context, layout):
        userpref = context.user_preferences
        view = userpref.view

        layout.prop(view, "pie_animation_timeout")
        layout.prop(view, "pie_initial_timeout")
        layout.prop(view, "pie_menu_radius")
        layout.prop(view, "pie_menu_threshold")
        layout.prop(view, "pie_menu_confirm")


class USERPREF_PT_interface_templates(PreferencePanel):
    bl_label = "Templates"
    bl_options = {'DEFAULT_CLOSED'}

    @classmethod
    def poll(cls, context):
        userpref = context.user_preferences
        return (userpref.active_section == 'INTERFACE')

    def draw_props(self, context, layout):
        userpref = context.user_preferences
        view = userpref.view

        layout.label(text="Options intended for use with app-templates only.")
        layout.prop(view, "show_layout_ui")


class USERPREF_PT_edit_undo(PreferencePanel):
    bl_label = "Undo"

    @classmethod
    def poll(cls, context):
        userpref = context.user_preferences
        return (userpref.active_section == 'EDITING')

    def draw_props(self, context, layout):
        userpref = context.user_preferences
        edit = userpref.edit

        layout.prop(edit, "undo_steps", text="Steps")
        layout.prop(edit, "undo_memory_limit", text="Memory Limit")
        layout.prop(edit, "use_global_undo")


class USERPREF_PT_edit_objects(PreferencePanel):
    bl_label = "Objects"
    bl_options = {'DEFAULT_CLOSED'}

    @classmethod
    def poll(cls, context):
        userpref = context.user_preferences
        return (userpref.active_section == 'EDITING')

    def draw_props(self, context, layout):
        userpref = context.user_preferences
        edit = userpref.edit

        layout.prop(edit, "material_link", text="Link Materials to")
        layout.prop(edit, "object_align", text="Align New Objects to")
        layout.prop(edit, "use_enter_edit_mode", text="Enter Edit Mode for New Objects")


class USERPREF_PT_edit_gpencil(PreferencePanel):
    bl_label = "Grease Pencil"
    bl_options = {'DEFAULT_CLOSED'}

    @classmethod
    def poll(cls, context):
        userpref = context.user_preferences
        return (userpref.active_section == 'EDITING')

    def draw_props(self, context, layout):
        userpref = context.user_preferences
        edit = userpref.edit

        layout.prop(edit, "grease_pencil_manhattan_distance", text="Manhattan Distance")
        layout.prop(edit, "grease_pencil_euclidean_distance", text="Euclidean Distance")

class USERPREF_PT_edit_annotations(PreferencePanel):
    bl_label = "Annotations"
    bl_options = {'DEFAULT_CLOSED'}

    @classmethod
    def poll(cls, context):
        userpref = context.user_preferences
        return (userpref.active_section == 'EDITING')

    def draw_props(self, context, layout):
        userpref = context.user_preferences
        edit = userpref.edit

        layout.prop(edit, "grease_pencil_default_color", text="Default Color")
        layout.prop(edit, "grease_pencil_eraser_radius", text="Eraser Radius")
        layout.prop(edit, "use_grease_pencil_simplify_stroke", text="Simplify Stroke")


class USERPREF_PT_edit_animation(PreferencePanel):
    bl_label = "Animation"

    @classmethod
    def poll(cls, context):
        userpref = context.user_preferences
        return (userpref.active_section == 'EDITING')

    def draw_props(self, context, layout):
        userpref = context.user_preferences
        edit = userpref.edit

        layout.prop(edit, "use_negative_frames")

        layout.prop(edit, "use_visual_keying")
        layout.prop(edit, "use_keyframe_insert_needed", text="Only Insert Needed")

class USERPREF_PT_edit_animation_autokey(PreferencePanel):
    bl_label = "Auto-Keyframing"
    bl_options = {'DEFAULT_CLOSED'}
    bl_parent_id = "USERPREF_PT_edit_animation"

    def draw_header(self, context):
        userpref = context.user_preferences
        edit = userpref.edit

        self.layout.prop(edit, "use_auto_keying", text="")

    def draw_props(self, context, layout):
        userpref = context.user_preferences
        edit = userpref.edit

        layout.prop(edit, "use_auto_keying_warning")
        layout.prop(edit, "use_keyframe_insert_available", text="Only Insert Available")

class USERPREF_PT_edit_animation_fcurves(PreferencePanel):
    bl_label = "F-Curves"
    bl_options = {'DEFAULT_CLOSED'}
    bl_parent_id = "USERPREF_PT_edit_animation"

    def draw_props(self, context, layout):
        userpref = context.user_preferences
        edit = userpref.edit

        layout.prop(edit, "keyframe_new_interpolation_type", text="Default Interpolation")
        layout.prop(edit, "keyframe_new_handle_type", text="Default Handles")
        layout.prop(edit, "use_insertkey_xyz_to_rgb", text="XYZ to RGB")

        layout.separator()

        layout.prop(edit, "fcurve_unselected_alpha", text="F-Curve Visibility")


class USERPREF_PT_edit_transform(PreferencePanel):
    bl_label = "Transform"
    bl_options = {'DEFAULT_CLOSED'}

    @classmethod
    def poll(cls, context):
        userpref = context.user_preferences
        return (userpref.active_section == 'EDITING')

    def draw_props(self, context, layout):
        userpref = context.user_preferences
        edit = userpref.edit

        layout.prop(edit, "use_drag_immediately")
        layout.prop(edit, "use_numeric_input_advanced")


class USERPREF_PT_edit_duplicate_data(PreferencePanel):
    bl_label = "Duplicate Data"
    bl_options = {'DEFAULT_CLOSED'}

    @classmethod
    def poll(cls, context):
        userpref = context.user_preferences
        return (userpref.active_section == 'EDITING')

    def draw_props(self, context, layout):
        userpref = context.user_preferences
        edit = userpref.edit

        flow = layout.grid_flow(row_major=True, columns=0, even_columns=False, even_rows=False, align=True)

        col = flow.column()
        col.prop(edit, "use_duplicate_action", text="Action")
        col.prop(edit, "use_duplicate_armature", text="Armature")
        col.prop(edit, "use_duplicate_curve", text="Curve")
        #col.prop(edit, "use_duplicate_fcurve", text="F-Curve")
        col.prop(edit, "use_duplicate_light", text="Light")
        col = flow.column()
        col.prop(edit, "use_duplicate_material", text="Material")
        col.prop(edit, "use_duplicate_mesh", text="Mesh")
        col.prop(edit, "use_duplicate_metaball", text="Metaball")
        col.prop(edit, "use_duplicate_particle", text="Particle")
        col = flow.column()
        col.prop(edit, "use_duplicate_surface", text="Surface")
        col.prop(edit, "use_duplicate_text", text="Text")
        col.prop(edit, "use_duplicate_texture", text="Texture")


class USERPREF_PT_edit_misc(PreferencePanel):
    bl_label = "Miscellaneous"
    bl_options = {'DEFAULT_CLOSED'}

    @classmethod
    def poll(cls, context):
        userpref = context.user_preferences
        return (userpref.active_section == 'EDITING')

    def draw_props(self, context, layout):
        userpref = context.user_preferences
        edit = userpref.edit

        layout.prop(edit, "sculpt_paint_overlay_color", text="Sculpt Overlay Color")
        layout.prop(edit, "node_margin", text="Node Editor Auto-offset Margin")


class USERPREF_PT_interface_system_sound(PreferencePanel):
    bl_label = "Sound"
    bl_options = {'DEFAULT_CLOSED'}

    @classmethod
    def poll(cls, context):
        userpref = context.user_preferences
        return (userpref.active_section == 'SYSTEM_GENERAL')

    def draw_props(self, context, layout):
        userpref = context.user_preferences
        system = userpref.system

        layout.prop(system, "audio_device", expand=False)
        sub = layout.column()
        sub.active = system.audio_device not in {'NONE', 'Null'}
        #sub.prop(system, "use_preview_images")
        sub.prop(system, "audio_channels", text="Channels")
        sub.prop(system, "audio_mixing_buffer", text="Mixing Buffer")
        sub.prop(system, "audio_sample_rate", text="Sample Rate")
        sub.prop(system, "audio_sample_format", text="Sample Format")


class USERPREF_PT_interface_system_compute_device(PreferencePanel):
    bl_label = "Cycles Compute Device"
    bl_options = {'DEFAULT_CLOSED'}

    @classmethod
    def poll(cls, context):
        userpref = context.user_preferences
        return (userpref.active_section == 'SYSTEM_GENERAL')

    def draw_props(self, context, layout):
        userpref = context.user_preferences
        system = userpref.system

        col = layout.column()

        if bpy.app.build_options.cycles:
            addon = userpref.addons.get("cycles")
            if addon is not None:
                addon.preferences.draw_impl(col, context)
            del addon

        # NOTE: Disabled for until GPU side of OpenSubdiv is brought back.
        # if hasattr(system, "opensubdiv_compute_type"):
        #     col.label(text="OpenSubdiv compute:")
        #     col.row().prop(system, "opensubdiv_compute_type", text="")


class USERPREF_PT_interface_system_opengl(PreferencePanel):
    bl_label = "OpenGL"

    @classmethod
    def poll(cls, context):
        userpref = context.user_preferences
        return (userpref.active_section == 'SYSTEM_GENERAL')

    def draw_props(self, context, layout):
        import sys
        userpref = context.user_preferences
        system = userpref.system

        layout.prop(system, "gl_clip_alpha", slider=True)
        layout.prop(system, "gpu_viewport_quality")

        layout.prop(system, "anisotropic_filter")

        layout.prop(system, "multi_sample", text="Multisampling")
        layout.prop(system, "gpencil_multi_sample", text="Grease Pencil Multisampling")

        if sys.platform == "linux" and system.multi_sample != 'NONE':
            layout.label(text="Might fail for Mesh editing selection!")
            layout.separator()

        layout.prop(system, "use_region_overlap")
        layout.prop(system, "use_gpu_mipmap")
        layout.prop(system, "use_16bit_textures")

class USERPREF_PT_interface_system_opengl_selection(PreferencePanel):
    bl_label = "Selection"
    bl_parent_id = "USERPREF_PT_interface_system_opengl"
    bl_options = {'DEFAULT_CLOSED'}

    def draw_props(self, context, layout):
        userpref = context.user_preferences
        system = userpref.system

        layout.prop(system, "select_method", text="Selection Method")
        layout.prop(system, "use_select_pick_depth")


class USERPREF_PT_interface_system_textures(PreferencePanel):
    bl_label = "Textures"
    bl_options = {'DEFAULT_CLOSED'}

    @classmethod
    def poll(cls, context):
        userpref = context.user_preferences
        return (userpref.active_section == 'SYSTEM_GENERAL')

    def draw_props(self, context, layout):
        userpref = context.user_preferences
        system = userpref.system

        layout.prop(system, "gl_texture_limit", text="Limit Size")
        layout.prop(system, "texture_time_out", text="Time Out")
        layout.prop(system, "texture_collection_rate", text="Collection Rate")

        layout.separator()

        layout.prop(system, "image_draw_method", text="Image Display Method")


class USERPREF_PT_interface_system_text(PreferencePanel):
    bl_label = "Text"
    bl_options = {'DEFAULT_CLOSED'}

    @classmethod
    def poll(cls, context):
        userpref = context.user_preferences
        return (userpref.active_section == 'SYSTEM_GENERAL')

    def draw_props(self, context, layout):
        userpref = context.user_preferences
        system = userpref.system

        layout.prop(system, "use_text_antialiasing", text="Anti-aliasing")
        sub = layout.column()
        sub.active = system.use_text_antialiasing
        sub.prop(system, "text_hinting", text="Hinting")

        layout.prop(system, "font_path_ui")
        layout.prop(system, "font_path_ui_mono")


class USERPREF_PT_interface_system_text_translate(PreferencePanel):
    bl_label = "Translate UI"
    bl_options = {'DEFAULT_CLOSED'}
    bl_parent_id = "USERPREF_PT_interface_system_text"

    @classmethod
    def poll(cls, context):
        userpref = context.user_preferences
        if bpy.app.build_options.international:
            return (userpref.active_section == 'SYSTEM_GENERAL')

    def draw_header(self, context):
        userpref = context.user_preferences
        system = userpref.system

        self.layout.prop(system, "use_international_fonts", text="")

    def draw_props(self, context, layout):
        userpref = context.user_preferences
        system = userpref.system

        layout.active = system.use_international_fonts
        layout.prop(system, "language")

        layout.prop(system, "use_translate_tooltips", text="Translate Tooltips")
        layout.prop(system, "use_translate_interface", text="Translate Interface")
        layout.prop(system, "use_translate_new_dataname", text="Translate New Data")


class USERPREF_PT_interface_system_misc(PreferencePanel):
    bl_label = "Misc"
    bl_options = {'DEFAULT_CLOSED'}

    @classmethod
    def poll(cls, context):
        userpref = context.user_preferences
        return (userpref.active_section == 'SYSTEM_GENERAL')

    def draw_props(self, context, layout):
        userpref = context.user_preferences
        system = userpref.system

        layout.prop(system, "color_picker_type")

        layout.prop(system, "memory_cache_limit")

        layout.prop(system, "scrollback", text="Console Scrollback")


class USERPREF_MT_interface_theme_presets(Menu):
    bl_label = "Presets"
    preset_subdir = "interface_theme"
    preset_operator = "script.execute_preset"
    preset_type = 'XML'
    preset_xml_map = (
        ("user_preferences.themes[0]", "Theme"),
        ("user_preferences.ui_styles[0]", "ThemeStyle"),
    )
    draw = Menu.draw_preset

    def reset_cb(context):
        bpy.ops.ui.reset_default_theme()


class USERPREF_PT_theme(Panel):
    bl_space_type = 'USER_PREFERENCES'
    bl_label = "Themes"
    bl_region_type = 'WINDOW'
    bl_options = {'HIDE_HEADER'}

    # not essential, hard-coded UI delimiters for the theme layout
    ui_delimiters = {
        'VIEW_3D': {
            "text_grease_pencil",
            "text_keyframe",
            "speaker",
            "freestyle_face_mark",
            "split_normal",
            "bone_solid",
            "paint_curve_pivot",
        },
        'GRAPH_EDITOR': {
            "handle_vertex_select",
        },
        'IMAGE_EDITOR': {
            "paint_curve_pivot",
        },
        'NODE_EDITOR': {
            "layout_node",
        },
        'CLIP_EDITOR': {
            "handle_vertex_select",
        }
    }

    @staticmethod
    def _theme_generic(layout, themedata, theme_area):

        col = layout.column()

        def theme_generic_recurse(data):
            col.label(text=data.rna_type.name)
            row = col.row()
            subsplit = row.split(factor=0.95)

            padding1 = subsplit.split(factor=0.15)
            padding1.column()

            subsplit = row.split(factor=0.85)

            padding2 = subsplit.split(factor=0.15)
            padding2.column()

            colsub_pair = padding1.column(), padding2.column()

            props_type = {}

            for i, prop in enumerate(data.rna_type.properties):
                if prop.identifier == "rna_type":
                    continue

                props_type.setdefault((prop.type, prop.subtype), []).append(prop)

            th_delimiters = USERPREF_PT_theme.ui_delimiters.get(theme_area)
            for props_type, props_ls in sorted(props_type.items()):
                if props_type[0] == 'POINTER':
                    for i, prop in enumerate(props_ls):
                        theme_generic_recurse(getattr(data, prop.identifier))
                else:
                    if th_delimiters is None:
                        # simple, no delimiters
                        for i, prop in enumerate(props_ls):
                            colsub_pair[i % 2].row().prop(data, prop.identifier)
                    else:
                        # add hard coded delimiters
                        i = 0
                        for prop in props_ls:
                            colsub = colsub_pair[i]
                            colsub.row().prop(data, prop.identifier)
                            i = (i + 1) % 2
                            if prop.identifier in th_delimiters:
                                if i:
                                    colsub = colsub_pair[1]
                                    colsub.row().label(text="")
                                colsub_pair[0].row().label(text="")
                                colsub_pair[1].row().label(text="")
                                i = 0

        theme_generic_recurse(themedata)

    @staticmethod
    def _theme_widget_style(layout, widget_style):

        col = layout.column()

        row = col.row()

        row.column().prop(widget_style, "outline")
        row.column().prop(widget_style, "item", slider=True)
        row.column().prop(widget_style, "inner", slider=True)
        row.column().prop(widget_style, "inner_sel", slider=True)
        row.column().prop(widget_style, "text")
        row.column().prop(widget_style, "text_sel")

        col.separator()

        row = col.row()

        row.prop(widget_style, "roundness")
        row.prop(widget_style, "show_shaded")

        rowsub = row.row(align=True)
        rowsub.active = widget_style.show_shaded
        rowsub.prop(widget_style, "shadetop")
        rowsub.prop(widget_style, "shadedown")


    @staticmethod
    def _ui_font_style(layout, font_style):

        split = layout.split()

        col = split.column()
        col.label(text="Kerning Style:")
        col.row().prop(font_style, "font_kerning_style", expand=True)
        col.prop(font_style, "points")

        col = split.column()
        col.label(text="Shadow Offset:")
        col.prop(font_style, "shadow_offset_x", text="X")
        col.prop(font_style, "shadow_offset_y", text="Y")

        col = split.column()
        col.prop(font_style, "shadow")
        col.prop(font_style, "shadow_alpha")
        col.prop(font_style, "shadow_value")

        layout.separator()

    @classmethod
    def poll(cls, context):
        userpref = context.user_preferences
        return (userpref.active_section == 'THEMES')

    def draw(self, context):
        layout = self.layout

        theme = context.user_preferences.themes[0]

        row = layout.row()

        row.operator("wm.theme_install", text="Install", icon='FILEBROWSER')
        row.operator("ui.reset_default_theme", text="Reset", icon='LOOP_BACK')

        subrow = row.row(align=True)
        subrow.menu("USERPREF_MT_interface_theme_presets", text=USERPREF_MT_interface_theme_presets.bl_label)
        subrow.operator("wm.interface_theme_preset_add", text="", icon='ADD')
        subrow.operator("wm.interface_theme_preset_add", text="", icon='REMOVE').remove_active = True

        row.prop(theme, "theme_area", text="")


        if theme.theme_area == 'USER_INTERFACE':
            col = layout.column()
            ui = theme.user_interface

            col.label(text="Regular:")
            self._theme_widget_style(col, ui.wcol_regular)

            col.label(text="Tool:")
            self._theme_widget_style(col, ui.wcol_tool)

            col.label(text="Toolbar Item:")
            self._theme_widget_style(col, ui.wcol_toolbar_item)

            col.label(text="Radio Buttons:")
            self._theme_widget_style(col, ui.wcol_radio)

            col.label(text="Text:")
            self._theme_widget_style(col, ui.wcol_text)

            col.label(text="Option:")
            self._theme_widget_style(col, ui.wcol_option)

            col.label(text="Toggle:")
            self._theme_widget_style(col, ui.wcol_toggle)

            col.label(text="Number Field:")
            self._theme_widget_style(col, ui.wcol_num)

            col.label(text="Value Slider:")
            self._theme_widget_style(col, ui.wcol_numslider)

            col.label(text="Box:")
            self._theme_widget_style(col, ui.wcol_box)

            col.label(text="Menu:")
            self._theme_widget_style(col, ui.wcol_menu)

            col.label(text="Pie Menu:")
            self._theme_widget_style(col, ui.wcol_pie_menu)

            col.label(text="Pulldown:")
            self._theme_widget_style(col, ui.wcol_pulldown)

            col.label(text="Menu Back:")
            self._theme_widget_style(col, ui.wcol_menu_back)

            col.label(text="Tooltip:")
            self._theme_widget_style(col, ui.wcol_tooltip)

            col.label(text="Menu Item:")
            self._theme_widget_style(col, ui.wcol_menu_item)

            col.label(text="Scroll Bar:")
            self._theme_widget_style(col, ui.wcol_scroll)

            col.label(text="Progress Bar:")
            self._theme_widget_style(col, ui.wcol_progress)

            col.label(text="List Item:")
            self._theme_widget_style(col, ui.wcol_list_item)

            col.label(text="Tab:")
            self._theme_widget_style(col, ui.wcol_tab)

            ui_state = theme.user_interface.wcol_state
            col.label(text="State:")

            row = col.row()

            subsplit = row.split(factor=0.95)

            padding = subsplit.split(factor=0.15)
            colsub = padding.column()
            colsub = padding.column()
            colsub.row().prop(ui_state, "inner_anim")
            colsub.row().prop(ui_state, "inner_anim_sel")
            colsub.row().prop(ui_state, "inner_driven")
            colsub.row().prop(ui_state, "inner_driven_sel")
            colsub.row().prop(ui_state, "blend")

            subsplit = row.split(factor=0.85)

            padding = subsplit.split(factor=0.15)
            colsub = padding.column()
            colsub = padding.column()
            colsub.row().prop(ui_state, "inner_key")
            colsub.row().prop(ui_state, "inner_key_sel")
            colsub.row().prop(ui_state, "inner_overridden")
            colsub.row().prop(ui_state, "inner_overridden_sel")

            col.separator()
            col.separator()

            col.label(text="Styles:")

            row = col.row()

            subsplit = row.split(factor=0.95)

            padding = subsplit.split(factor=0.15)
            colsub = padding.column()
            colsub = padding.column()
            colsub.row().prop(ui, "menu_shadow_fac")
            colsub.row().prop(ui, "icon_alpha")
            colsub.row().prop(ui, "icon_saturation")
            colsub.row().prop(ui, "editor_outline")

            subsplit = row.split(factor=0.85)

            padding = subsplit.split(factor=0.15)
            colsub = padding.column()
            colsub = padding.column()
            colsub.row().prop(ui, "menu_shadow_width")
            colsub.row().prop(ui, "widget_emboss")

            col.separator()
            col.separator()

            col.label(text="Axis & Gizmo Colors:")

            row = col.row()

            subsplit = row.split(factor=0.95)

            padding = subsplit.split(factor=0.15)
            colsub = padding.column()
            colsub = padding.column()
            colsub.row().prop(ui, "axis_x")
            colsub.row().prop(ui, "axis_y")
            colsub.row().prop(ui, "axis_z")

            subsplit = row.split(factor=0.85)

            padding = subsplit.split(factor=0.15)
            colsub = padding.column()
            colsub = padding.column()
            colsub.row().prop(ui, "gizmo_primary")
            colsub.row().prop(ui, "gizmo_secondary")
            colsub.row().prop(ui, "gizmo_a")
            colsub.row().prop(ui, "gizmo_b")

            col.separator()
            col.separator()

            col.label(text="Icon Colors:")

            row = col.row()

            subsplit = row.split(factor=0.95)

            padding = subsplit.split(factor=0.15)
            colsub = padding.column()
            colsub = padding.column()
            colsub.row().prop(ui, "icon_collection")
            colsub.row().prop(ui, "icon_object")
            colsub.row().prop(ui, "icon_object_data")

            subsplit = row.split(factor=0.85)

            padding = subsplit.split(factor=0.15)
            colsub = padding.column()
            colsub = padding.column()
            colsub.row().prop(ui, "icon_modifier")
            colsub.row().prop(ui, "icon_shading")

            col.separator()
            col.separator()
        elif theme.theme_area == 'BONE_COLOR_SETS':
            col = layout.column()

            for i, ui in enumerate(theme.bone_color_sets, 1):
                col.label(text=iface_(f"Color Set {i:d}"), translate=False)

                row = col.row()

                subsplit = row.split(factor=0.95)

                padding = subsplit.split(factor=0.15)
                colsub = padding.column()
                colsub = padding.column()
                colsub.row().prop(ui, "normal")
                colsub.row().prop(ui, "select")
                colsub.row().prop(ui, "active")

                subsplit = row.split(factor=0.85)

                padding = subsplit.split(factor=0.15)
                colsub = padding.column()
                colsub = padding.column()
                colsub.row().prop(ui, "show_colored_constraints")
        elif theme.theme_area == 'STYLE':
            col = layout.column()

            style = context.user_preferences.ui_styles[0]

            col.label(text="Panel Title:")
            self._ui_font_style(col, style.panel_title)

            col.separator()

            col.label(text="Widget:")
            self._ui_font_style(col, style.widget)

            col.separator()

            col.label(text="Widget Label:")
            self._ui_font_style(col, style.widget_label)
        else:
            self._theme_generic(layout, getattr(theme, theme.theme_area.lower()), theme.theme_area)


class USERPREF_PT_file_paths(PreferencePanel):
    bl_label = "File Paths"

    @classmethod
    def poll(cls, context):
        userpref = context.user_preferences
        return (userpref.active_section == 'SYSTEM_FILES')

    def draw_props(self, context, layout):
        userpref = context.user_preferences
        paths = userpref.filepaths
        system = userpref.system

        layout.prop(paths, "render_output_directory", text="Render Output")
        layout.prop(paths, "render_cache_directory", text="Render Cache")
        layout.prop(paths, "font_directory", text="Fonts")
        layout.prop(paths, "texture_directory", text="Textures")
        layout.prop(paths, "script_directory", text="Scripts")
        layout.prop(paths, "sound_directory", text="Sounds")
        layout.prop(paths, "temporary_directory", text="Temp")
        layout.prop(paths, "i18n_branches_directory", text="I18n Branches")
        layout.prop(paths, "image_editor", text="Image Editor")
        layout.prop(paths, "animation_player_preset", text="Playback Preset")
        layout.prop(paths, "animation_player", text="Animation Player")


class USERPREF_PT_file_autorun(PreferencePanel):
    bl_label = "Auto Run Python Scripts"
    bl_options = {'DEFAULT_CLOSED'}

    @classmethod
    def poll(cls, context):
        userpref = context.user_preferences
        return (userpref.active_section == 'SYSTEM_FILES')

    def draw_header(self, context):
        userpref = context.user_preferences
        paths = userpref.filepaths
        system = userpref.system

        self.layout.prop(system, "use_scripts_auto_execute", text="")

    def draw_props(self, context, layout):
        userpref = context.user_preferences
        paths = userpref.filepaths
        system = userpref.system

        if system.use_scripts_auto_execute:
            box = layout.box()
            row = box.row()
            row.label(text="Excluded Paths:")
            row.operator("wm.userpref_autoexec_path_add", text="", icon='ADD', emboss=False)
            for i, path_cmp in enumerate(userpref.autoexec_paths):
                row = box.row()
                row.prop(path_cmp, "path", text="")
                row.prop(path_cmp, "use_glob", text="", icon='FILTER')
                row.operator("wm.userpref_autoexec_path_remove", text="", icon='X', emboss=False).index = i


class USERPREF_PT_file_saveload(PreferencePanel):
    bl_label = "Save & Load"
    bl_options = {'DEFAULT_CLOSED'}

    @classmethod
    def poll(cls, context):
        userpref = context.user_preferences
        return (userpref.active_section == 'SYSTEM_FILES')

    def draw_props(self, context, layout):
        userpref = context.user_preferences
        paths = userpref.filepaths
        system = userpref.system


        layout.prop(paths, "use_relative_paths")
        layout.prop(paths, "use_file_compression")
        layout.prop(paths, "use_load_ui")
        layout.prop(paths, "use_filter_files")
        layout.prop(paths, "show_hidden_files_datablocks")
        layout.prop(paths, "hide_recent_locations")
        layout.prop(paths, "hide_system_bookmarks")
        layout.prop(paths, "show_thumbnails")
        layout.prop(paths, "use_save_preview_images")

        layout.separator()

        layout.prop(paths, "save_version")
        layout.prop(paths, "recent_files")


class USERPREF_PT_file_saveload_autosave(PreferencePanel):
    bl_label = "Auto Save"
    bl_parent_id = "USERPREF_PT_file_saveload"
    bl_options = {'DEFAULT_CLOSED'}

    def draw_props(self, context, layout):
        userpref = context.user_preferences
        paths = userpref.filepaths
        system = userpref.system

        layout.prop(paths, "use_keep_session")
        layout.prop(paths, "use_auto_save_temporary_files")
        sub = layout.column()
        sub.active = paths.use_auto_save_temporary_files
        sub.prop(paths, "auto_save_time", text="Timer (mins)")


class USERPREF_PT_file_saveload_texteditor(PreferencePanel):
    bl_label = "Text Editor"
    bl_parent_id = "USERPREF_PT_file_saveload"
    bl_options = {'DEFAULT_CLOSED'}

    def draw_props(self, context, layout):
        userpref = context.user_preferences
        paths = userpref.filepaths
        system = userpref.system

        layout.prop(system, "use_tabs_as_spaces")
        layout.prop(system, "author", text="Author")


class USERPREF_MT_ndof_settings(Menu):
    # accessed from the window key-bindings in C (only)
    bl_label = "3D Mouse Settings"

    def draw(self, context):
        layout = self.layout

        input_prefs = context.user_preferences.inputs

        is_view3d = context.space_data.type == 'VIEW_3D'

        layout.prop(input_prefs, "ndof_sensitivity")
        layout.prop(input_prefs, "ndof_orbit_sensitivity")
        layout.prop(input_prefs, "ndof_deadzone")

        if is_view3d:
            layout.separator()
            layout.prop(input_prefs, "ndof_show_guide")

            layout.separator()
            layout.label(text="Orbit Style")
            layout.row().prop(input_prefs, "ndof_view_navigate_method", text="")
            layout.row().prop(input_prefs, "ndof_view_rotate_method", text="")
            layout.separator()
            layout.label(text="Orbit Options")
            layout.prop(input_prefs, "ndof_rotx_invert_axis")
            layout.prop(input_prefs, "ndof_roty_invert_axis")
            layout.prop(input_prefs, "ndof_rotz_invert_axis")

        # view2d use pan/zoom
        layout.separator()
        layout.label(text="Pan Options")
        layout.prop(input_prefs, "ndof_panx_invert_axis")
        layout.prop(input_prefs, "ndof_pany_invert_axis")
        layout.prop(input_prefs, "ndof_panz_invert_axis")
        layout.prop(input_prefs, "ndof_pan_yz_swap_axis")

        layout.label(text="Zoom Options")
        layout.prop(input_prefs, "ndof_zoom_invert")

        if is_view3d:
            layout.separator()
            layout.label(text="Fly/Walk Options")
            layout.prop(input_prefs, "ndof_fly_helicopter", icon='NDOF_FLY')
            layout.prop(input_prefs, "ndof_lock_horizon", icon='NDOF_DOM')


class USERPREF_PT_input_mouse(PreferencePanel):
    bl_label = "Mouse"

    @classmethod
    def poll(cls, context):
        userpref = context.user_preferences
        return (userpref.active_section == 'INPUT')

    def draw_props(self, context, layout):
        userpref = context.user_preferences
        inputs = userpref.inputs

        layout.prop(inputs, "drag_threshold")
        layout.prop(inputs, "tweak_threshold")
        layout.prop(inputs, "mouse_double_click_time", text="Double Click Speed")
        layout.prop(inputs, "use_mouse_emulate_3_button")
        layout.prop(inputs, "use_mouse_continuous")


class USERPREF_PT_input_view(PreferencePanel):
    bl_label = "View"

    @classmethod
    def poll(cls, context):
        userpref = context.user_preferences
        return (userpref.active_section == 'INPUT')

    def draw_props(self, context, layout):
        import sys
        userpref = context.user_preferences
        inputs = userpref.inputs

        layout.row().prop(inputs, "view_rotate_method", expand=True)

        layout.row().prop(inputs, "view_zoom_method", text="Zoom Method")
        if inputs.view_zoom_method in {'DOLLY', 'CONTINUE'}:
            layout.row().prop(inputs, "view_zoom_axis", expand=True)
            layout.prop(inputs, "invert_mouse_zoom", text="Invert Mouse Zoom Direction")

        layout.prop(inputs, "invert_zoom_wheel", text="Invert Wheel Zoom Direction")
        #sub.prop(view, "wheel_scroll_lines", text="Scroll Lines")

        if sys.platform == "darwin":
            sub = layout.column()
            sub.prop(inputs, "use_trackpad_natural", text="Natural Trackpad Direction")

class USERPREF_PT_input_view_fly_walk(PreferencePanel):
    bl_label = "Fly & Walk"
    bl_parent_id = "USERPREF_PT_input_view"
    bl_options = {'DEFAULT_CLOSED'}

    def draw_props(self, context, layout):
        userpref = context.user_preferences
        inputs = userpref.inputs

        layout.row().prop(inputs, "navigation_mode", expand=True)

        layout.label(text="Walk Navigation:")

        walk = inputs.walk_navigation

        layout.prop(walk, "use_mouse_reverse")
        layout.prop(walk, "mouse_speed")
        layout.prop(walk, "teleport_time")

        sub = layout.column(align=True)
        sub.prop(walk, "walk_speed")
        sub.prop(walk, "walk_speed_factor")


class USERPREF_PT_input_view_fly_walk_gravity(PreferencePanel):
    bl_label = "Gravity"
    bl_parent_id = "USERPREF_PT_input_view_fly_walk"

    def draw_header(self, context):
        userpref = context.user_preferences
        inputs = userpref.inputs
        walk = inputs.walk_navigation

        self.layout.prop(walk, "use_gravity", text="")

    def draw_props(self, context, layout):
        userpref = context.user_preferences
        inputs = userpref.inputs
        walk = inputs.walk_navigation

        layout.active = walk.use_gravity
        layout.prop(walk, "view_height")
        layout.prop(walk, "jump_height")


class USERPREF_PT_input_devices(PreferencePanel):
    bl_label = "Devices"

    @classmethod
    def poll(cls, context):
        userpref = context.user_preferences
        return (userpref.active_section == 'INPUT')

    def draw_props(self, context, layout):
        userpref = context.user_preferences
        inputs = userpref.inputs

        layout.prop(inputs, "use_emulate_numpad")

class USERPREF_PT_input_devices_tablet(PreferencePanel):
    bl_label = "Tablet"
    bl_parent_id = "USERPREF_PT_input_devices"
    bl_options = {'DEFAULT_CLOSED'}

    def draw_props(self, context, layout):
        userpref = context.user_preferences
        inputs = userpref.inputs

        layout.prop(inputs, "pressure_threshold_max")
        layout.prop(inputs, "pressure_softness")


class USERPREF_PT_input_devices_ndof(PreferencePanel):
    bl_label = "NDOF"
    bl_parent_id = "USERPREF_PT_input_devices"
    bl_options = {'DEFAULT_CLOSED'}

    @classmethod
    def poll(cls, context):
        userpref = context.user_preferences
        inputs = userpref.inputs
        if inputs.use_ndof:
            return (userpref.active_section == 'INPUT')

    def draw_props(self, context, layout):
        userpref = context.user_preferences
        inputs = userpref.inputs

        layout.prop(inputs, "ndof_sensitivity", text="Pan Sensitivity")
        layout.prop(inputs, "ndof_orbit_sensitivity", text="Orbit Sensitivity")
        layout.prop(inputs, "ndof_deadzone", text="Deadzone")

        layout.separator()

        layout.row().prop(inputs, "ndof_view_navigate_method", expand=True)
        layout.row().prop(inputs, "ndof_view_rotate_method", expand=True)


class USERPREF_MT_keyconfigs(Menu):
    bl_label = "KeyPresets"
    preset_subdir = "keyconfig"
    preset_operator = "wm.keyconfig_activate"

    def draw(self, context):
        Menu.draw_preset(self, context)


class USERPREF_PT_keymap(Panel):
    bl_space_type = 'USER_PREFERENCES'
    bl_label = "Keymap"
    bl_region_type = 'WINDOW'
    bl_options = {'HIDE_HEADER'}

    @classmethod
    def poll(cls, context):
        userpref = context.user_preferences
        return (userpref.active_section == 'KEYMAP')

    @staticmethod
    def draw_input_prefs(inputs, layout):
        import sys


    def draw(self, context):
        from rna_keymap_ui import draw_keymaps

        layout = self.layout

        #import time

        #start = time.time()

        userpref = context.user_preferences
        inputs = userpref.inputs

        col = layout.column()

        # Input settings
        self.draw_input_prefs(inputs, col)

        # Keymap Settings
        draw_keymaps(context, col)

        #print("runtime", time.time() - start)


class USERPREF_MT_addons_online_resources(Menu):
    bl_label = "Online Resources"

    # menu to open web-pages with addons development guides
    def draw(self, context):
        layout = self.layout

        layout.operator(
            "wm.url_open", text="Add-ons Catalog", icon='URL',
        ).url = "http://wiki.blender.org/index.php/Extensions:2.6/Py/Scripts"

        layout.separator()

        layout.operator(
            "wm.url_open", text="How to share your add-on", icon='URL',
        ).url = "http://wiki.blender.org/index.php/Dev:Py/Sharing"
        layout.operator(
            "wm.url_open", text="Add-on Guidelines", icon='URL',
        ).url = "http://wiki.blender.org/index.php/Dev:2.5/Py/Scripts/Guidelines/Addons"
        layout.operator(
            "wm.url_open", text="API Concepts", icon='URL',
        ).url = bpy.types.WM_OT_doc_view._prefix + "/info_quickstart.html"
        layout.operator(
            "wm.url_open", text="Add-on Tutorial", icon='URL',
        ).url = bpy.types.WM_OT_doc_view._prefix + "/info_tutorial_addon.html"


class USERPREF_PT_addons(Panel):
    bl_space_type = 'USER_PREFERENCES'
    bl_label = "Add-ons"
    bl_region_type = 'WINDOW'
    bl_options = {'HIDE_HEADER'}

    _support_icon_mapping = {
        'OFFICIAL': 'FILE_BLEND',
        'COMMUNITY': 'COMMUNITY',
        'TESTING': 'EXPERIMENTAL',
    }

    @classmethod
    def poll(cls, context):
        userpref = context.user_preferences
        return (userpref.active_section == 'ADDONS')

    @staticmethod
    def is_user_addon(mod, user_addon_paths):
        import os

        if not user_addon_paths:
            for path in (
                    bpy.utils.script_path_user(),
                    bpy.utils.script_path_pref(),
            ):
                if path is not None:
                    user_addon_paths.append(os.path.join(path, "addons"))

        for path in user_addon_paths:
            if bpy.path.is_subdir(mod.__file__, path):
                return True
        return False

    @staticmethod
    def draw_error(layout, message):
        lines = message.split("\n")
        box = layout.box()
        sub = box.row()
        sub.label(text=lines[0])
        sub.label(icon='ERROR')
        for l in lines[1:]:
            box.label(text=l)

    def draw(self, context):
        import os
        import addon_utils

        layout = self.layout

        userpref = context.user_preferences
        used_ext = {ext.module for ext in userpref.addons}

        addon_user_dirs = tuple(
            p for p in (
                os.path.join(userpref.filepaths.script_directory, "addons"),
                bpy.utils.user_resource('SCRIPTS', "addons"),
            )
            if p
        )

        # Development option for 2.8x, don't show users bundled addons
        # unless they have been updated for 2.8x.
        # Developers can turn them on with '--debug'
        show_official_27x_addons = bpy.app.debug

        # collect the categories that can be filtered on
        addons = [
            (mod, addon_utils.module_bl_info(mod))
            for mod in addon_utils.modules(refresh=False)
        ]

        row = layout.row()
        row.operator("wm.addon_install", icon='ADD', text="Install")
        row.operator("wm.addon_refresh", icon='FILE_REFRESH', text="Refresh")
        row.menu("USERPREF_MT_addons_online_resources", text="Resources")
        row = layout.row()
        row.prop(context.window_manager, "addon_support", expand=True)
        row.prop(context.window_manager, "addon_filter", text="")
        row.prop(context.window_manager, "addon_search", text="", icon='VIEWZOOM')

        col = layout.column()

        # set in addon_utils.modules_refresh()
        if addon_utils.error_duplicates:
            box = col.box()
            row = box.row()
            row.label(text="Multiple add-ons with the same name found!")
            row.label(icon='ERROR')
            box.label(text="Delete one of each pair to resolve:")
            for (addon_name, addon_file, addon_path) in addon_utils.error_duplicates:
                box.separator()
                sub_col = box.column(align=True)
                sub_col.label(text=addon_name + ":")
                sub_col.label(text="    " + addon_file)
                sub_col.label(text="    " + addon_path)

        if addon_utils.error_encoding:
            self.draw_error(
                col,
                "One or more addons do not have UTF-8 encoding\n"
                "(see console for details)",
            )

        filter = context.window_manager.addon_filter
        search = context.window_manager.addon_search.lower()
        support = context.window_manager.addon_support

        # initialized on demand
        user_addon_paths = []

        for mod, info in addons:
            module_name = mod.__name__

            is_enabled = module_name in used_ext

            if info["support"] not in support:
                continue

            # check if addon should be visible with current filters
            if (
                    (filter == "All") or
                    (filter == info["category"]) or
                    (filter == "Enabled" and is_enabled) or
                    (filter == "Disabled" and not is_enabled) or
                    (filter == "User" and (mod.__file__.startswith(addon_user_dirs)))
            ):
                if search and search not in info["name"].lower():
                    if info["author"]:
                        if search not in info["author"].lower():
                            continue
                    else:
                        continue

                # Skip 2.7x add-ons included with Blender, unless in debug mode.
                is_addon_27x = info.get("blender", (0,)) < (2, 80)
                if (
                        is_addon_27x and
                        (not show_official_27x_addons) and
                        (not mod.__file__.startswith(addon_user_dirs))
                ):
                    continue

                # Addon UI Code
                col_box = col.column()
                box = col_box.box()
                colsub = box.column()
                row = colsub.row(align=True)

                row.operator(
                    "wm.addon_expand",
                    icon='DISCLOSURE_TRI_DOWN' if info["show_expanded"] else 'DISCLOSURE_TRI_RIGHT',
                    emboss=False,
                ).module = module_name

                row.operator(
                    "wm.addon_disable" if is_enabled else "wm.addon_enable",
                    icon='CHECKBOX_HLT' if is_enabled else 'CHECKBOX_DEHLT', text="",
                    emboss=False,
                ).module = module_name

                sub = row.row()
                sub.active = is_enabled
                sub.label(text="%s: %s" % (info["category"], info["name"]))

                # WARNING: 2.8x exception, may be removed
                # use disabled state for old add-ons, chances are they are broken.
                if is_addon_27x:
                    sub.label(text="upgrade to 2.8x required")
                    sub.label(icon='ERROR')
                # Remove code above after 2.8x migration is complete.
                elif info["warning"]:
                    sub.label(icon='ERROR')

                # icon showing support level.
                sub.label(icon=self._support_icon_mapping.get(info["support"], 'QUESTION'))

                # Expanded UI (only if additional info is available)
                if info["show_expanded"]:
                    if info["description"]:
                        split = colsub.row().split(factor=0.15)
                        split.label(text="Description:")
                        split.label(text=info["description"])
                    if info["location"]:
                        split = colsub.row().split(factor=0.15)
                        split.label(text="Location:")
                        split.label(text=info["location"])
                    if mod:
                        split = colsub.row().split(factor=0.15)
                        split.label(text="File:")
                        split.label(text=mod.__file__, translate=False)
                    if info["author"]:
                        split = colsub.row().split(factor=0.15)
                        split.label(text="Author:")
                        split.label(text=info["author"], translate=False)
                    if info["version"]:
                        split = colsub.row().split(factor=0.15)
                        split.label(text="Version:")
                        split.label(text=".".join(str(x) for x in info["version"]), translate=False)
                    if info["warning"]:
                        split = colsub.row().split(factor=0.15)
                        split.label(text="Warning:")
                        split.label(text="  " + info["warning"], icon='ERROR')

                    user_addon = USERPREF_PT_addons.is_user_addon(mod, user_addon_paths)
                    tot_row = bool(info["wiki_url"]) + bool(user_addon)

                    if tot_row:
                        split = colsub.row().split(factor=0.15)
                        split.label(text="Internet:")
                        if info["wiki_url"]:
                            split.operator(
                                "wm.url_open", text="Documentation", icon='HELP',
                            ).url = info["wiki_url"]
                        # Only add "Report a Bug" button if tracker_url is set
                        # or the add-on is bundled (use official tracker then).
                        if info.get("tracker_url") or not user_addon:
                            split.operator(
                                "wm.url_open", text="Report a Bug", icon='URL',
                            ).url = info.get(
                                "tracker_url",
                                "https://developer.blender.org/maniphest/task/edit/form/2",
                            )
                        if user_addon:
                            split.operator(
                                "wm.addon_remove", text="Remove", icon='CANCEL',
                            ).module = mod.__name__

                        for _ in range(4 - tot_row):
                            split.separator()

                    # Show addon user preferences
                    if is_enabled:
                        addon_preferences = userpref.addons[module_name].preferences
                        if addon_preferences is not None:
                            draw = getattr(addon_preferences, "draw", None)
                            if draw is not None:
                                addon_preferences_class = type(addon_preferences)
                                box_prefs = col_box.box()
                                box_prefs.label(text="Preferences:")
                                addon_preferences_class.layout = box_prefs
                                try:
                                    draw(context)
                                except:
                                    import traceback
                                    traceback.print_exc()
                                    box_prefs.label(text="Error (see console)", icon='ERROR')
                                del addon_preferences_class.layout

        # Append missing scripts
        # First collect scripts that are used but have no script file.
        module_names = {mod.__name__ for mod, info in addons}
        missing_modules = {ext for ext in used_ext if ext not in module_names}

        if missing_modules and filter in {"All", "Enabled"}:
            col.column().separator()
            col.column().label(text="Missing script files")

            module_names = {mod.__name__ for mod, info in addons}
            for module_name in sorted(missing_modules):
                is_enabled = module_name in used_ext
                # Addon UI Code
                box = col.column().box()
                colsub = box.column()
                row = colsub.row(align=True)

                row.label(text="", icon='ERROR')

                if is_enabled:
                    row.operator(
                        "wm.addon_disable", icon='CHECKBOX_HLT', text="", emboss=False,
                    ).module = module_name

                row.label(text=module_name, translate=False)


class StudioLightPanelMixin():
    bl_space_type = 'USER_PREFERENCES'
    bl_region_type = 'WINDOW'

    @classmethod
    def poll(cls, context):
        userpref = context.user_preferences
        return (userpref.active_section == 'LIGHTS')

    def _get_lights(self, userpref):
        return [light for light in userpref.studio_lights if light.is_user_defined and light.type == self.sl_type]

    def draw(self, context):
        layout = self.layout
        userpref = context.user_preferences
        lights = self._get_lights(userpref)

        self.draw_light_list(layout, lights)

    def draw_light_list(self, layout, lights):
        if lights:
            flow = layout.column_flow(columns=4)
            for studio_light in lights:
                self.draw_studio_light(flow, studio_light)
        else:
            layout.label(text="No custom {} configured".format(self.bl_label))

    def draw_studio_light(self, layout, studio_light):
        box = layout.box()
        row = box.row()

        row.template_icon(layout.icon(studio_light), scale=6.0)
        col = row.column()
        op = col.operator('wm.studiolight_uninstall', text="", icon='REMOVE')
        op.index = studio_light.index

        if studio_light.type == 'STUDIO':
            op = col.operator('wm.studiolight_copy_settings', text="", icon='IMPORT')
            op.index = studio_light.index

        box.label(text=studio_light.name)


class USERPREF_PT_studiolight_matcaps(Panel, StudioLightPanelMixin):
    bl_label = "MatCaps"
    sl_type = 'MATCAP'


class USERPREF_PT_studiolight_world(Panel, StudioLightPanelMixin):
    bl_label = "LookDev HDRIs"
    sl_type = 'WORLD'


class USERPREF_PT_studiolight_lights(Panel, StudioLightPanelMixin):
    bl_label = "Studio Lights"
    sl_type = 'STUDIO'


class USERPREF_PT_studiolight_light_editor(Panel):
    bl_label = "Studio Light Editor"
    bl_parent_id = "USERPREF_PT_studiolight_lights"
    bl_space_type = 'USER_PREFERENCES'
    bl_region_type = 'WINDOW'

    def opengl_light_buttons(self, layout, light):

        col = layout.column()
        col.active = light.use

        col.prop(light, "use", text="Use Light")
        col.prop(light, "diffuse_color", text="Diffuse")
        col.prop(light, "specular_color", text="Specular")
        col.prop(light, "smooth")
        col.prop(light, "direction")

    def draw(self, context):
        layout = self.layout

        userpref = context.user_preferences
        system = userpref.system

        row = layout.row()
        row.prop(system, "edit_studio_light", toggle=True)
        row.operator('wm.studiolight_new', text="Save as Studio light", icon="FILE_TICK")

        layout.separator()

        layout.use_property_split = True
        column = layout.split()
        column.active = system.edit_studio_light

        light = system.solid_lights[0]
        colsplit = column.split(factor=0.85)
        self.opengl_light_buttons(colsplit, light)

        light = system.solid_lights[1]
        colsplit = column.split(factor=0.85)
        self.opengl_light_buttons(colsplit, light)

        light = system.solid_lights[2]
        colsplit = column.split(factor=0.85)
        self.opengl_light_buttons(colsplit, light)

        light = system.solid_lights[3]
        self.opengl_light_buttons(column, light)

        layout.separator()

        layout.prop(system, "light_ambient")


classes = (
    USERPREF_HT_header,
    USERPREF_PT_navigation,

    USERPREF_PT_interface_display,
    USERPREF_PT_interface_display_info,
    USERPREF_PT_interface_view_manipulation,
    USERPREF_PT_interface_viewports,
    USERPREF_PT_interface_viewports_3d,
    USERPREF_PT_interface_viewports_3d_weight_paint,
    USERPREF_PT_interface_viewports_2d,
    USERPREF_PT_interface_menus,
    USERPREF_PT_interface_menus_mouse_over,
    USERPREF_PT_interface_menus_pie,
    USERPREF_PT_interface_develop,
    USERPREF_PT_interface_templates,

    USERPREF_PT_edit_undo,
    USERPREF_PT_edit_objects,
    USERPREF_PT_edit_animation,
    USERPREF_PT_edit_animation_autokey,
    USERPREF_PT_edit_animation_fcurves,
    USERPREF_PT_edit_transform,
    USERPREF_PT_edit_duplicate_data,
    USERPREF_PT_edit_gpencil,
    USERPREF_PT_edit_annotations,
    USERPREF_PT_edit_misc,

    USERPREF_PT_interface_system_opengl,
    USERPREF_PT_interface_system_opengl_selection,
    USERPREF_PT_interface_system_sound,
    USERPREF_PT_interface_system_compute_device,
    USERPREF_PT_interface_system_textures,
    USERPREF_PT_interface_system_text,
    USERPREF_PT_interface_system_text_translate,

    USERPREF_MT_interface_theme_presets,

    USERPREF_PT_theme,

    USERPREF_PT_file_paths,
    USERPREF_PT_file_autorun,
    USERPREF_PT_file_saveload,
    USERPREF_PT_file_saveload_autosave,
    USERPREF_PT_file_saveload_texteditor,

    USERPREF_MT_ndof_settings,
    USERPREF_MT_keyconfigs,

    USERPREF_PT_input_mouse,
    USERPREF_PT_input_view,
    USERPREF_PT_input_view_fly_walk,
    USERPREF_PT_input_view_fly_walk_gravity,
    USERPREF_PT_input_devices,
    USERPREF_PT_input_devices_tablet,
    USERPREF_PT_input_devices_ndof,

    USERPREF_PT_keymap,
    USERPREF_MT_addons_online_resources,
    USERPREF_PT_addons,

    USERPREF_PT_studiolight_lights,
    USERPREF_PT_studiolight_light_editor,
    USERPREF_PT_studiolight_matcaps,
    USERPREF_PT_studiolight_world,
)

if __name__ == "__main__":  # only for live edit.
    from bpy.utils import register_class
    for cls in classes:
        register_class(cls)
