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
    Panel,
    UIList,
)

from rna_prop_ui import PropertyPanel
from bl_operators.presets import PresetMenu

from .properties_physics_common import (
    point_cache_ui,
    effector_weights_ui,
)


class SCENE_UL_keying_set_paths(UIList):
    def draw_item(self, context, layout, data, item, icon, active_data, active_propname, index):
        # assert(isinstance(item, bpy.types.KeyingSetPath)
        kspath = item
        icon = layout.enum_item_icon(kspath, "id_type", kspath.id_type)
        if self.layout_type in {'DEFAULT', 'COMPACT'}:
            # Do not make this one editable in uiList for now...
            layout.label(text=kspath.data_path, translate=False, icon_value=icon)
        elif self.layout_type == 'GRID':
            layout.alignment = 'CENTER'
            layout.label(text="", icon_value=icon)


class SceneButtonsPanel:
    bl_space_type = 'PROPERTIES'
    bl_region_type = 'WINDOW'
    bl_context = "scene"

    @classmethod
    def poll(cls, context):
        return (context.engine in cls.COMPAT_ENGINES)


class SCENE_PT_scene(SceneButtonsPanel, Panel):
    bl_label = "Scene"
    COMPAT_ENGINES = {'BLENDER_RENDER', 'BLENDER_EEVEE', 'BLENDER_WORKBENCH'}

    def draw(self, context):
        layout = self.layout
        layout.use_property_split = True
        layout.use_property_decorate = False

        scene = context.scene

        layout.prop(scene, "camera")
        layout.prop(scene, "background_set")
        layout.prop(scene, "active_clip")


class SCENE_PT_unit(SceneButtonsPanel, Panel):
    bl_label = "Units"
    bl_options = {'DEFAULT_CLOSED'}
    COMPAT_ENGINES = {'BLENDER_RENDER', 'BLENDER_EEVEE', 'BLENDER_WORKBENCH'}

    def draw(self, context):
        layout = self.layout

        unit = context.scene.unit_settings

        layout.use_property_split = True
        layout.use_property_decorate = False

        layout.prop(unit, "system")

        col = layout.column()
        col.enabled = unit.system != 'NONE'
        col.prop(unit, "scale_length")
        col.prop(unit, "use_separate")

        col = layout.column()
        col.prop(unit, "system_rotation", text="Rotation")
        subcol = col.column()
        subcol.enabled = unit.system != 'NONE'
        subcol.prop(unit, "length_unit", text="Length")
        subcol.prop(unit, "mass_unit", text="Mass")
        subcol.prop(unit, "time_unit", text="Time")


class SceneKeyingSetsPanel:

    @staticmethod
    def draw_keyframing_settings(context, layout, ks, ksp):
        SceneKeyingSetsPanel._draw_keyframing_setting(
            context, layout, ks, ksp, "Needed",
            "use_insertkey_override_needed", "use_insertkey_needed",
            userpref_fallback="use_keyframe_insert_needed",
        )
        SceneKeyingSetsPanel._draw_keyframing_setting(
            context, layout, ks, ksp, "Visual",
            "use_insertkey_override_visual", "use_insertkey_visual",
            userpref_fallback="use_visual_keying",
        )
        SceneKeyingSetsPanel._draw_keyframing_setting(
            context, layout, ks, ksp, "XYZ to RGB",
            "use_insertkey_override_xyz_to_rgb", "use_insertkey_xyz_to_rgb",
        )

    @staticmethod
    def _draw_keyframing_setting(context, layout, ks, ksp, label, toggle_prop, prop, userpref_fallback=None):
        if ksp:
            item = ksp

            if getattr(ks, toggle_prop):
                owner = ks
                propname = prop
            else:
                owner = context.user_preferences.edit
                if userpref_fallback:
                    propname = userpref_fallback
                else:
                    propname = prop
        else:
            item = ks

            owner = context.user_preferences.edit
            if userpref_fallback:
                propname = userpref_fallback
            else:
                propname = prop

        row = layout.row(align=True)

        subrow = row.row(align=True)
        subrow.active = getattr(item, toggle_prop)

        if subrow.active:
            subrow.prop(item, prop, text=label)
        else:
            subrow.prop(owner, propname, text=label)

        row.prop(item, toggle_prop, text="", icon='STYLUS_PRESSURE', toggle=True)  # XXX: needs dedicated icon


class SCENE_PT_keying_sets(SceneButtonsPanel, SceneKeyingSetsPanel, Panel):
    bl_label = "Keying Sets"
    bl_options = {'DEFAULT_CLOSED'}
    COMPAT_ENGINES = {'BLENDER_RENDER', 'BLENDER_EEVEE', 'BLENDER_WORKBENCH'}

    def draw(self, context):
        layout = self.layout

        scene = context.scene

        row = layout.row()

        col = row.column()
        col.template_list("UI_UL_list", "keying_sets", scene, "keying_sets", scene.keying_sets, "active_index", rows=1)

        col = row.column(align=True)
        col.operator("anim.keying_set_add", icon='ZOOMIN', text="")
        col.operator("anim.keying_set_remove", icon='ZOOMOUT', text="")

        layout.use_property_split = True
        layout.use_property_decorate = False  # No animation.

        flow = layout.grid_flow(row_major=False, columns=0, even_columns=False, even_rows=False, align=False)

        ks = scene.keying_sets.active
        if ks and ks.is_path_absolute:
            col = flow.column()
            col.prop(ks, "bl_description")

            subcol = flow.column()
            subcol.operator_context = 'INVOKE_DEFAULT'
            subcol.operator("anim.keying_set_export", text="Export to File").filepath = "keyingset.py"


class SCENE_PT_keyframing_settings(SceneButtonsPanel, SceneKeyingSetsPanel, Panel):
    bl_label = "Keyframing Settings"
    bl_parent_id = "SCENE_PT_keying_sets"
    COMPAT_ENGINES = {'BLENDER_RENDER', 'BLENDER_EEVEE', "BLENDER_LANPR"}

    @classmethod
    def poll(cls, context):
        ks = context.scene.keying_sets.active
        return (ks and ks.is_path_absolute)

    def draw(self, context):
        layout = self.layout
        layout.use_property_split = True
        layout.use_property_decorate = False  # No animation.

        scene = context.scene
        ks = scene.keying_sets.active

        flow = layout.grid_flow(row_major=True, columns=0, even_columns=False, even_rows=False, align=True)

        col = flow.column(align=True)
        col.alignment = 'RIGHT'
        col.label(text="General Override")

        self.draw_keyframing_settings(context, col, ks, None)

        ksp = ks.paths.active
        if ksp:
            col.separator()

            col = flow.column(align=True)
            col.alignment = 'RIGHT'
            col.label(text="Active Set Override")

            self.draw_keyframing_settings(context, col, ks, ksp)


class SCENE_PT_keying_set_paths(SceneButtonsPanel, SceneKeyingSetsPanel, Panel):
    bl_label = "Active Keying Set"
    bl_parent_id = "SCENE_PT_keying_sets"
    COMPAT_ENGINES = {'BLENDER_RENDER', 'BLENDER_EEVEE', 'BLENDER_WORKBENCH'}

    @classmethod
    def poll(cls, context):
        ks = context.scene.keying_sets.active
        return (ks and ks.is_path_absolute)

    def draw(self, context):
        layout = self.layout

        scene = context.scene
        ks = scene.keying_sets.active

        row = layout.row()
        row.label(text="Paths:")

        row = layout.row()

        col = row.column()
        col.template_list("SCENE_UL_keying_set_paths", "", ks, "paths", ks.paths, "active_index", rows=1)

        col = row.column(align=True)
        col.operator("anim.keying_set_path_add", icon='ZOOMIN', text="")
        col.operator("anim.keying_set_path_remove", icon='ZOOMOUT', text="")

        # TODO: 1) the template_any_ID needs to be fixed for the text alignment.
        #       2) use_property_decorate has to properly skip the non animatable properties.
        #          Properties affected with needless draw:
        #          group_method, template_any_ID dropdown, use_entire_array

        layout.use_property_split = True
        layout.use_property_decorate = False  # No animation (remove this later on).

        flow = layout.grid_flow(row_major=False, columns=0, even_columns=False, even_rows=False, align=True)

        ksp = ks.paths.active
        if ksp:
            col = flow.column(align=True)
            col.alignment = 'RIGHT'

            col.template_any_ID(ksp, "id", "id_type", text="Target ID-Block")

            col.separator()

            col.template_path_builder(ksp, "data_path", ksp.id, text="Data Path")

            col = flow.column()

            col.prop(ksp, "use_entire_array", text="Array All Items")

            if not ksp.use_entire_array:
                col.prop(ksp, "array_index", text="Index")

            col.separator()

            col.prop(ksp, "group_method", text="F-Curve Grouping")
            if ksp.group_method == 'NAMED':
                col.prop(ksp, "group")


class SCENE_PT_color_management(SceneButtonsPanel, Panel):
    bl_label = "Color Management"
    bl_options = {'DEFAULT_CLOSED'}
    COMPAT_ENGINES = {'BLENDER_RENDER', 'BLENDER_EEVEE', 'BLENDER_OPENGL', 'BLENDER_LANPR'}

    def draw(self, context):
        layout = self.layout
        layout.use_property_split = True

        scene = context.scene
        view = scene.view_settings

        flow = layout.grid_flow(row_major=True, columns=0, even_columns=False, even_rows=False, align=True)

        col = flow.column()
        col.prop(scene.display_settings, "display_device")

        col.separator()

        col.prop(view, "view_transform")
        col.prop(view, "look")

        col = flow.column()
        col.prop(view, "exposure")
        col.prop(view, "gamma")

        col.separator()

        col.prop(scene.sequencer_colorspace_settings, "name", text="Sequencer")


class SCENE_PT_color_management_curves(SceneButtonsPanel, Panel):
    bl_label = "Use Curves"
    bl_parent_id = "SCENE_PT_color_management"
    bl_options = {'DEFAULT_CLOSED'}
    COMPAT_ENGINES = {'BLENDER_RENDER', 'BLENDER_EEVEE', 'BLENDER_OPENGL', 'BLENDER_LANPR'}

    def draw_header(self, context):

        scene = context.scene
        view = scene.view_settings

        self.layout.prop(view, "use_curve_mapping", text="")

    def draw(self, context):
        layout = self.layout

        scene = context.scene
        view = scene.view_settings

        layout.use_property_split = False
        layout.enabled = view.use_curve_mapping

        layout.template_curve_mapping(view, "curve_mapping", levels=True)


class SCENE_PT_audio(SceneButtonsPanel, Panel):
    bl_label = "Audio"
    bl_options = {'DEFAULT_CLOSED'}
    COMPAT_ENGINES = {'BLENDER_RENDER', 'BLENDER_EEVEE', 'BLENDER_WORKBENCH'}

    def draw(self, context):
        layout = self.layout
        layout.use_property_split = True

        scene = context.scene
        rd = context.scene.render
        ffmpeg = rd.ffmpeg

        flow = layout.grid_flow(row_major=True, columns=0, even_columns=True, even_rows=False, align=True)

        col = flow.column()
        col.prop(scene, "audio_volume")

        col.separator()

        col.prop(scene, "audio_distance_model")
        col.prop(ffmpeg, "audio_channels")

        col.separator()

        col = flow.column()
        col.prop(ffmpeg, "audio_mixrate", text="Sample Rate")

        col.separator()

        col = col.column(align=True)
        col.prop(scene, "audio_doppler_speed", text="Doppler Speed")
        col.prop(scene, "audio_doppler_factor", text="Doppler Factor")

        col.separator()

        layout.operator("sound.bake_animation")


class SCENE_PT_physics(SceneButtonsPanel, Panel):
    bl_label = "Gravity"
    bl_options = {'DEFAULT_CLOSED'}
    COMPAT_ENGINES = {'BLENDER_RENDER', 'BLENDER_EEVEE', 'BLENDER_WORKBENCH'}

    def draw_header(self, context):
        self.layout.prop(context.scene, "use_gravity", text="")

    def draw(self, context):
        layout = self.layout
        layout.use_property_split = True

        scene = context.scene

        layout.active = scene.use_gravity

        layout.prop(scene, "gravity")


class SCENE_PT_rigid_body_world(SceneButtonsPanel, Panel):
    bl_label = "Rigid Body World"
    bl_options = {'DEFAULT_CLOSED'}
    COMPAT_ENGINES = {'BLENDER_RENDER', 'BLENDER_EEVEE', 'BLENDER_WORKBENCH'}

    @classmethod
    def poll(cls, context):
        return (context.engine in cls.COMPAT_ENGINES)

    def draw_header(self, context):
        scene = context.scene
        rbw = scene.rigidbody_world
        if rbw is not None:
            self.layout.prop(rbw, "enabled", text="")

    def draw(self, context):
        layout = self.layout
        layout.use_property_split = True

        scene = context.scene
        rbw = scene.rigidbody_world

        if rbw is None:
            layout.operator("rigidbody.world_add")
        else:
            layout.operator("rigidbody.world_remove")


class SCENE_PT_rigid_body_world_settings(SceneButtonsPanel, Panel):
    bl_label = "Settings"
    bl_parent_id = "SCENE_PT_rigid_body_world"
    COMPAT_ENGINES = {'BLENDER_RENDER', 'BLENDER_EEVEE'}

    @classmethod
    def poll(cls, context):
        scene = context.scene
        return scene and scene.rigidbody_world and (context.engine in cls.COMPAT_ENGINES)

    def draw(self, context):
        layout = self.layout
        layout.use_property_split = True

        scene = context.scene
        rbw = scene.rigidbody_world

        if rbw:
            flow = layout.grid_flow(row_major=True, columns=0, even_columns=True, even_rows=False, align=True)

            col = flow.column()
            col.active = rbw.enabled

            col = col.column()
            col.prop(rbw, "collection")
            col.prop(rbw, "constraints")

            col = col.column()
            col.prop(rbw, "time_scale", text="Speed")

            col = flow.column()
            col.active = rbw.enabled
            col.prop(rbw, "use_split_impulse")

            col = col.column()
            col.prop(rbw, "steps_per_second", text="Steps Per Second")
            col.prop(rbw, "solver_iterations", text="Solver Iterations")


class SCENE_PT_rigid_body_cache(SceneButtonsPanel, Panel):
    bl_label = "Cache"
    bl_parent_id = "SCENE_PT_rigid_body_world"
    bl_options = {'DEFAULT_CLOSED'}
    COMPAT_ENGINES = {'BLENDER_RENDER', 'BLENDER_EEVEE', 'BLENDER_WORKBENCH'}

    @classmethod
    def poll(cls, context):
        scene = context.scene
        return scene and scene.rigidbody_world and (context.engine in cls.COMPAT_ENGINES)

    def draw(self, context):
        scene = context.scene
        rbw = scene.rigidbody_world

        point_cache_ui(self, context, rbw.point_cache, rbw.point_cache.is_baked is False and rbw.enabled, 'RIGID_BODY')


class SCENE_PT_rigid_body_field_weights(SceneButtonsPanel, Panel):
    bl_label = "Field Weights"
    bl_parent_id = "SCENE_PT_rigid_body_world"
    bl_options = {'DEFAULT_CLOSED'}
    COMPAT_ENGINES = {'BLENDER_RENDER', 'BLENDER_EEVEE', 'BLENDER_WORKBENCH'}

    @classmethod
    def poll(cls, context):
        scene = context.scene
        return scene and scene.rigidbody_world and (context.engine in cls.COMPAT_ENGINES)

    def draw(self, context):
        scene = context.scene
        rbw = scene.rigidbody_world

        effector_weights_ui(self, context, rbw.effector_weights, 'RIGID_BODY')


class SCENE_PT_simplify(SceneButtonsPanel, Panel):
    bl_label = "Simplify"
    bl_options = {'DEFAULT_CLOSED'}
    COMPAT_ENGINES = {'BLENDER_RENDER', 'BLENDER_EEVEE', 'BLENDER_OPENGL', 'BLENDER_LANPR'}

    def draw_header(self, context):
        rd = context.scene.render
        self.layout.prop(rd, "use_simplify", text="")

    def draw(self, context):
        pass


class SCENE_PT_simplify_viewport(SceneButtonsPanel, Panel):
    bl_label = "Viewport"
    bl_parent_id = "SCENE_PT_simplify"
    COMPAT_ENGINES = {'BLENDER_RENDER', 'BLENDER_EEVEE', 'BLENDER_OPENGL', 'BLENDER_LANPR'}

    def draw(self, context):
        layout = self.layout
        layout.use_property_split = True

        rd = context.scene.render

        layout.active = rd.use_simplify

        flow = layout.grid_flow(row_major=True, columns=0, even_columns=False, even_rows=False, align=True)

        col = flow.column()
        col.prop(rd, "simplify_subdivision", text="Max Subdivision")

        col = flow.column()
        col.prop(rd, "simplify_child_particles", text="Max Child Particles")


class SCENE_PT_simplify_render(SceneButtonsPanel, Panel):
    bl_label = "Render"
    bl_parent_id = "SCENE_PT_simplify"
    COMPAT_ENGINES = {'BLENDER_RENDER', 'BLENDER_EEVEE', 'BLENDER_OPENGL', 'BLENDER_LANPR'}

    def draw(self, context):
        layout = self.layout
        layout.use_property_split = True

        rd = context.scene.render

        layout.active = rd.use_simplify

        flow = layout.grid_flow(row_major=True, columns=0, even_columns=False, even_rows=False, align=True)

        col = flow.column()
        col.prop(rd, "simplify_subdivision_render", text="Max Subdivision")

        col = flow.column()
        col.prop(rd, "simplify_child_particles_render", text="Max Child Particles")


class SCENE_PT_simplify_greasepencil(SceneButtonsPanel, Panel):
    bl_label = "Grease Pencil"
    bl_parent_id = "SCENE_PT_simplify"
    COMPAT_ENGINES = {'BLENDER_RENDER', 'BLENDER_GAME', 'BLENDER_CLAY', 'BLENDER_EEVEE'}
    bl_options = {'DEFAULT_CLOSED'}

    def draw_header(self, context):
        rd = context.scene.render
        self.layout.prop(rd, "simplify_gpencil", text="")

    def draw(self, context):
        layout = self.layout
        layout.use_property_split = True

        rd = context.scene.render

        layout.active = rd.simplify_gpencil

        col = layout.column()
        col.prop(rd, "simplify_gpencil_onplay", text="Playback Only")
        col.prop(rd, "simplify_gpencil_view_modifier", text="Modifiers")
        col.prop(rd, "simplify_gpencil_shader_fx", text="ShaderFX")

        col = layout.column(align=True)
        col.prop(rd, "simplify_gpencil_view_fill")
        sub = col.column()
        sub.active = rd.simplify_gpencil_view_fill
        sub.prop(rd, "simplify_gpencil_remove_lines", text="Lines")

class SCENE_PT_viewport_display(SceneButtonsPanel, Panel):
    bl_label = "Viewport Display"
    bl_options = {'DEFAULT_CLOSED'}

    @classmethod
    def poll(cls, context):
        return True

    def draw(self, context):
        layout = self.layout
        layout.use_property_split = True
        scene = context.scene
        col = layout.column()
        col.prop(scene.display, "light_direction")
        col.prop(scene.display, "shadow_shift")

class LANPR_linesets(UIList):
    def draw_item(self, context, layout, data, item, icon, active_data, active_propname, index):
        lineset = item
        if self.layout_type in {'DEFAULT', 'COMPACT'}:
            split = layout.split(factor=0.6)
            split.label(text="Layer")
            row = split.row(align=True)
            row.prop(lineset, "color", text="", icon_value=icon)
            row.prop(lineset, "thickness", text="", icon_value=icon)
        elif self.layout_type == 'GRID':
            layout.alignment = 'CENTER'
            layout.label("", icon_value=icon)

def lanpr_get_composition_scene(scene):
    n = scene.name+'_lanpr_comp'
    for s in bpy.data.scenes:
        if s.name == n: return s
    return None

def lanpr_is_composition_scene(scene):
    return scene.name.endswith('_lanpr_comp')

class SCENE_PT_lanpr(SceneButtonsPanel, Panel):
    COMPAT_ENGINES = {'BLENDER_RENDER', 'BLENDER_LANPR', 'BLENDER_OPENGL', 'BLENDER_EEVEE'}
    bl_label = "LANPR"
    bl_options = {'DEFAULT_CLOSED'}
    
    @classmethod
    def poll(cls, context):
        return True

    def draw(self, context):
        layout = self.layout
        scene = context.scene
        lanpr = scene.lanpr
        active_layer = lanpr.layers.active_layer 
        
        sc = lanpr_get_composition_scene(scene)
        
        if lanpr_is_composition_scene(scene):
            row = layout.row()
            row.scale_y=1.5
            row.operator("lanpr.goto_original_scene")

        if sc is not None:
            layout.label(text = 'You are adjusting values for LANPR compostion scene.')
            row = layout.row()
            row.scale_y=1.5
            row.operator("lanpr.goto_composition_scene")
            layout.operator("lanpr.remove_composition_scene")
            scene = sc
            lanpr = scene.lanpr
            active_layer = lanpr.layers.active_layer
            return
        elif scene.render.engine!='BLENDER_LANPR':
            layout.label(text = 'Select LANPR engine or use composition scene.')
            layout.operator("lanpr.make_composition_scene")
            return
            

        layout.prop(lanpr, "master_mode", expand=True) 

        if lanpr.master_mode == "DPIX" or lanpr.master_mode == "SOFTWARE":
            
            layout.prop(lanpr, "background_color")

            if lanpr.master_mode == "SOFTWARE":
                layout.label(text="Enable On Demand:")
                row = layout.row()
                row.prop(lanpr,"enable_intersections", text = "Intersection Lines")
                row.prop(lanpr,"enable_chaining", text = "Chaining (SLOW!)")
                layout.label(text="RUN:")
                layout.operator("scene.lanpr_calculate", icon='RENDER_STILL')

                split = layout.split(factor=0.7)
                col = split.column()
                col.label(text="Layer Composition:")
                col = split.column()
                col.operator("scene.lanpr_auto_create_line_layer", text = "Default")#, icon = "ZOOMIN")
                layout.template_list("LANPR_linesets", "", lanpr, "layers", lanpr.layers, "active_layer_index", rows=4)
                if active_layer:
                    split = layout.split()
                    col = split.column()
                    col.operator("scene.lanpr_add_line_layer")#icon="ZOOMIN")
                    col.operator("scene.lanpr_delete_line_layer")#, icon="ZOOMOUT")
                    col = split.column()
                    col.operator("scene.lanpr_move_line_layer").direction = "UP"
                    col.operator("scene.lanpr_move_line_layer").direction = "DOWN"
                    layout.operator("scene.lanpr_rebuild_all_commands")
                else:
                    layout.operator("scene.lanpr_add_line_layer")
            
            if active_layer:
                layout.label(text= "Normal:")
                layout.prop(active_layer,"normal_mode", expand = True)
                if active_layer.normal_mode != "DISABLED":
                    layout.prop(active_layer,"normal_control_object")
                    layout.prop(active_layer,"normal_effect_inverse", toggle = True)
                    layout.prop(active_layer,"normal_ramp_begin")
                    layout.prop(active_layer,"normal_ramp_end")
                    layout.prop(active_layer,"normal_thickness_begin", slider=True)
                    layout.prop(active_layer,"normal_thickness_end", slider=True)
            
        else:
            layout.label(text="Vectorization:")
            layout.prop(lanpr, "enable_vector_trace", expand = True)


class SCENE_PT_lanpr_line_types(SceneButtonsPanel, Panel):
    bl_label = "Types"
    bl_parent_id = "SCENE_PT_lanpr"
    COMPAT_ENGINES = {'BLENDER_RENDER', 'BLENDER_LANPR', 'BLENDER_OPENGL'}

    @classmethod
    def poll(cls, context):
        scene = context.scene
        lanpr = scene.lanpr
        active_layer = lanpr.layers.active_layer
        return active_layer and lanpr.master_mode != "SNAKE"

    def draw(self, context):
        layout = self.layout
        scene = context.scene
        lanpr = scene.lanpr
        active_layer = lanpr.layers.active_layer
        if active_layer and lanpr.master_mode == "DPIX":
            active_layer = lanpr.layers[0]

        layout.operator("scene.lanpr_enable_all_line_types")

        split = layout.split(factor=0.3)
        col = split.column()
        col.prop(active_layer, "enable_contour", text="Contour", toggle=True)
        col.prop(active_layer, "enable_crease", text="Crease", toggle=True)
        col.prop(active_layer, "enable_edge_mark", text="Mark", toggle=True)
        col.prop(active_layer, "enable_material_seperate", text="Material", toggle=True)
        col.prop(active_layer, "enable_intersection", text="Intersection", toggle=True)
        col = split.column()
        row = col.row(align = True)
        #row.enabled = active_layer.enable_contour this is always enabled now
        row.prop(active_layer, "color", text="")
        row.prop(active_layer, "thickness", text="")
        row = col.row(align = True)
        row.enabled = active_layer.enable_crease
        row.prop(active_layer, "crease_color", text="")
        row.prop(active_layer, "thickness_crease", text="")
        row = col.row(align = True)
        row.enabled = active_layer.enable_edge_mark
        row.prop(active_layer, "edge_mark_color", text="")
        row.prop(active_layer, "thickness_edge_mark", text="")
        row = col.row(align = True)
        row.enabled = active_layer.enable_material_seperate
        row.prop(active_layer, "material_color", text="")
        row.prop(active_layer, "thickness_material", text="")
        row = col.row(align = True)
        if lanpr.enable_intersections:
            row.enabled = active_layer.enable_intersection
            row.prop(active_layer, "intersection_color", text="")
            row.prop(active_layer, "thickness_intersection", text="")
        else:
            row.label(text= "Intersection Calculation Disabled")

        if lanpr.master_mode == "DPIX" and active_layer.enable_intersection:
            row = col.row(align = True)
            row.prop(lanpr,"enable_intersections", toggle = True, text = "Enable")
            if lanpr.enable_intersections:
                row.operator("scene.lanpr_calculate", text= "Recalculate")

        if lanpr.master_mode == "SOFTWARE":
            row = layout.row(align=True)
            row.prop(active_layer, "qi_begin")
            row.prop(active_layer, "qi_end")

class SCENE_PT_lanpr_line_components(SceneButtonsPanel, Panel):
    bl_label = "Including"
    bl_parent_id = "SCENE_PT_lanpr"
    COMPAT_ENGINES = {'BLENDER_RENDER', 'BLENDER_LANPR', 'BLENDER_OPENGL'}

    @classmethod
    def poll(cls, context):
        scene = context.scene
        lanpr = scene.lanpr
        active_layer = lanpr.layers.active_layer
        return active_layer and lanpr.master_mode == "SOFTWARE" and not lanpr.enable_chaining

    def draw(self, context):
        layout = self.layout
        scene = context.scene
        lanpr = scene.lanpr
        active_layer = lanpr.layers.active_layer

        layout.operator("scene.lanpr_add_line_component")#, icon = "ZOOMIN")
        
        i=0
        for c in active_layer.components:
            split = layout.split(factor=0.85)
            col = split.column()
            sp2 = col.split(factor=0.4)
            cl = sp2.column()
            cl.prop(c,"component_mode", text = "")
            cl = sp2.column()
            if c.component_mode == "OBJECT":
                cl.prop(c,"object_select", text = "")
            elif c.component_mode == "MATERIAL":
                cl.prop(c,"material_select", text = "")
            elif c.component_mode == "COLLECTION":
                cl.prop(c,"collection_select", text = "")
            col = split.column()
            col.operator("scene.lanpr_delete_line_component", text="").index=i
            i=i+1


class SCENE_PT_lanpr_line_effects(SceneButtonsPanel, Panel):
    bl_label = "Effects"
    bl_parent_id = "SCENE_PT_lanpr"
    COMPAT_ENGINES = {'BLENDER_RENDER', 'BLENDER_LANPR', 'BLENDER_OPENGL'}

    @classmethod
    def poll(cls, context):
        scene = context.scene
        lanpr = scene.lanpr
        return lanpr.master_mode == "DPIX"

    def draw(self, context):
        layout = self.layout
        scene = context.scene
        lanpr = scene.lanpr
        active_layer = lanpr.layers.active_layer

        row = layout.row(align = True)
        row.prop(lanpr, "crease_threshold")
        row.prop(lanpr, "crease_fade_threshold")
        row = layout.row(align = True)
        row.prop(lanpr, "depth_width_influence")
        row.prop(lanpr, "depth_width_curve")
        row = layout.row(align = True)
        row.prop(lanpr, "depth_alpha_influence")
        row.prop(lanpr, "depth_alpha_curve")

class SCENE_PT_lanpr_snake_sobel_parameters(SceneButtonsPanel, Panel):
    bl_label = "Sobel Parameters"
    bl_parent_id = "SCENE_PT_lanpr"
    COMPAT_ENGINES = {'BLENDER_RENDER', 'BLENDER_LANPR', 'BLENDER_OPENGL'}

    @classmethod
    def poll(cls, context):
        scene = context.scene
        lanpr = scene.lanpr
        return lanpr.master_mode == "SNAKE"

    def draw(self, context):
        layout = self.layout
        scene = context.scene
        lanpr = scene.lanpr
        layout.prop(lanpr, "depth_clamp")
        layout.prop(lanpr, "depth_strength")
        layout.prop(lanpr, "normal_clamp")
        layout.prop(lanpr, "normal_strength")
        if lanpr.enable_vector_trace == "DISABLED":
            layout.prop(lanpr, "display_thinning_result")

class SCENE_PT_lanpr_snake_settings(SceneButtonsPanel, Panel):
    bl_label = "Snake Settings"
    bl_parent_id = "SCENE_PT_lanpr"
    COMPAT_ENGINES = {'BLENDER_RENDER', 'BLENDER_LANPR', 'BLENDER_OPENGL'}

    @classmethod
    def poll(cls, context):
        scene = context.scene
        lanpr = scene.lanpr
        return lanpr.master_mode == "SNAKE" and lanpr.enable_vector_trace == "ENABLED"

    def draw(self, context):
        layout = self.layout
        scene = context.scene
        lanpr = scene.lanpr

        split = layout.split()
        col = split.column()
        col.prop(lanpr, "background_color")
        col = split.column()
        col.prop(lanpr, "line_color")
        
        layout.prop(lanpr, "line_thickness")

        split = layout.split()
        col = split.column()
        col.prop(lanpr, "depth_width_influence")
        col.prop(lanpr, "depth_alpha_influence")
        col = split.column()
        col.prop(lanpr, "depth_width_curve")
        col.prop(lanpr, "depth_alpha_curve")
        
        layout.label(text="Taper:")
        layout.prop(lanpr, "use_same_taper", expand = True)
        if lanpr.use_same_taper == "DISABLED":
            split = layout.split()
            col = split.column(align = True)
            col.label(text="Left:")
            col.prop(lanpr,"taper_left_distance")
            col.prop(lanpr,"taper_left_strength")
            col = split.column(align = True)
            col.label(text="Right:")
            col.prop(lanpr,"taper_right_distance")
            col.prop(lanpr,"taper_right_strength")
        else:
            split = layout.split()
            col = split.column(align = True)
            col.prop(lanpr,"taper_left_distance")
            col.prop(lanpr,"taper_left_strength") 

        layout.label(text="Tip Extend:")
        layout.prop(lanpr, "enable_tip_extend",  expand = True)
        if lanpr.enable_tip_extend == "ENABLED":
            layout.label(text="---INOP---")
            layout.prop(lanpr,"extend_length")

class SCENE_PT_lanpr_software_chain_styles(SceneButtonsPanel, Panel):
    bl_label = "Chain Styles"
    bl_parent_id = "SCENE_PT_lanpr"
    COMPAT_ENGINES = {'BLENDER_RENDER', 'BLENDER_LANPR', 'BLENDER_OPENGL'}

    @classmethod
    def poll(cls, context):
        scene = context.scene
        lanpr = scene.lanpr
        return lanpr.master_mode == "SOFTWARE" and lanpr.enable_chaining

    def draw(self, context):
        layout = self.layout
        scene = context.scene
        lanpr = scene.lanpr
        layout.label(text="Taper:")
        layout.prop(lanpr, "use_same_taper", expand = True)
        if lanpr.use_same_taper == "DISABLED":
            split = layout.split()
            col = split.column(align = True)
            col.label(text="Left:")
            col.prop(lanpr,"taper_left_distance")
            col.prop(lanpr,"taper_left_strength")
            col = split.column(align = True)
            col.label(text="Right:")
            col.prop(lanpr,"taper_right_distance")
            col.prop(lanpr,"taper_right_strength")
        else:
            split = layout.split()
            col = split.column(align = True)
            col.prop(lanpr,"taper_left_distance")
            col.prop(lanpr,"taper_left_strength") 

        layout.label(text="Tip Extend:")
        layout.prop(lanpr, "enable_tip_extend",  expand = True)
        if lanpr.enable_tip_extend == "ENABLED":
            layout.label(text="---INOP---")
            layout.prop(lanpr,"extend_length")


class SCENE_PT_viewport_display_ssao(SceneButtonsPanel, Panel):
    bl_label = "Screen Space Ambient Occlusion"
    bl_parent_id = "SCENE_PT_viewport_display"

    @classmethod
    def poll(cls, context):
        return True

    def draw(self, context):
        layout = self.layout
        layout.use_property_split = True
        scene = context.scene
        col = layout.column()
        col.prop(scene.display, "matcap_ssao_samples")
        col.prop(scene.display, "matcap_ssao_distance")
        col.prop(scene.display, "matcap_ssao_attenuation")


class SCENE_PT_custom_props(SceneButtonsPanel, PropertyPanel, Panel):
    COMPAT_ENGINES = {'BLENDER_RENDER', 'BLENDER_EEVEE', 'BLENDER_WORKBENCH'}
    _context_path = "scene"
    _property_type = bpy.types.Scene


classes = (
    SCENE_UL_keying_set_paths,
    SCENE_PT_scene,
    SCENE_PT_unit,
    SCENE_PT_keying_sets,
    SCENE_PT_keying_set_paths,
    SCENE_PT_keyframing_settings,
    SCENE_PT_color_management,
    SCENE_PT_color_management_curves,
    SCENE_PT_viewport_display,
    SCENE_PT_viewport_display_ssao,
    SCENE_PT_audio,
    SCENE_PT_physics,
    SCENE_PT_rigid_body_world,
    SCENE_PT_rigid_body_world_settings,
    SCENE_PT_rigid_body_cache,
    SCENE_PT_rigid_body_field_weights,
    SCENE_PT_simplify,
    SCENE_PT_simplify_viewport,
    SCENE_PT_simplify_render,
    SCENE_PT_simplify_greasepencil,
    SCENE_PT_custom_props,

    SCENE_PT_lanpr,
    SCENE_PT_lanpr_line_types,
    SCENE_PT_lanpr_line_components,
    SCENE_PT_lanpr_line_effects,
    SCENE_PT_lanpr_snake_sobel_parameters,
    SCENE_PT_lanpr_snake_settings,
    SCENE_PT_lanpr_software_chain_styles,

    LANPR_linesets,
)

if __name__ == "__main__":  # only for live edit.
    from bpy.utils import register_class
    for cls in classes:
        register_class(cls)
