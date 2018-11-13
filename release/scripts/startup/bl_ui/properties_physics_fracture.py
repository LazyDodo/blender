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
from bpy.types import Panel, Menu, UIList
from bpy.app.translations import pgettext_iface as iface_
from bl_operators.presets import PresetMenu


class FRACTURE_PT_presets(PresetMenu):
    bl_label = "Fracture Presets"
    preset_subdir = "fracture"
    preset_operator = "script.execute_preset"
    preset_add_operator = "fracture.preset_add"

class PhysicButtonsPanel():
    bl_space_type = 'PROPERTIES'
    bl_region_type = 'WINDOW'
    bl_context = "physics"

    @classmethod
    def poll(cls, context):
        ob = context.object
        return ob and (ob.type == 'MESH') and context.fracture

class PHYSICS_PT_fracture_anim_mesh(PhysicButtonsPanel, Panel):
    bl_label = "Animated Mesh"
    bl_options = {'DEFAULT_CLOSED'}
    bl_parent_id = 'PHYSICS_PT_fracture'

    def draw_header(self, context):
        md = context.fracture
        self.layout.prop(md, "use_animated_mesh", text="")

    def draw(self, context):
        layout = self.layout
        md = context.fracture
        layout.context_pointer_set("modifier", md)
        layout.active = md.use_animated_mesh
        row = layout.row()
        row.prop(md, "animated_mesh_limit")
        row.prop(md, "use_animated_mesh_rotation")
        row = layout.row()
        row.prop(md, "animated_mesh_input")
        row = layout.row()
        row.operator("object.fracture_anim_bind", text="Bind", icon="UV_VERTEXSEL")

class PHYSICS_PT_fracture_advanced(PhysicButtonsPanel, Panel):
    bl_label = "Advanced"
    bl_options = {'DEFAULT_CLOSED'}
    bl_parent_id = 'PHYSICS_PT_fracture'

    def draw(self, context):
        layout = self.layout
        md = context.fracture
        ob = context.object

        layout.label(text="Fracture Point Source")
        col = layout.column()
        col.context_pointer_set("modifier", md)
        col.prop(md, "point_source")
        if 'GRID' in md.point_source:
            sub = col.split(factor=0.33)
            sub.prop(md, "grid_resolution")
            sub.prop(md, "grid_offset")
            sub.prop(md, "grid_spacing")
        #if 'GREASE_PENCIL' in md.point_source:
        #    col.prop(md, "use_greasepencil_edges")
        #    col.prop(md, "grease_offset")
        #    col.prop(md, "grease_decimate")
        #    col.prop(md, "cutter_axis")
        if 'EXTRA_PARTICLES' in md.point_source or 'EXTRA_VERTS' in md.point_source:
            col.prop(md, "extra_group", text="Helpers")
        if 'CUSTOM' in md.point_source:
            col.prop(md, "cutter_group")
            if (md.cutter_group):
                col.prop(md, "keep_cutter_shards")
                col.label(text="Material Index Offset")
                row = col.row(align=True)
                row.prop(md, "material_offset_intersect", text="Intersect")
                row.prop(md, "material_offset_difference", text="Difference")
        row = col.row()
        row.prop(md, "dm_group", text="Pack Group")
        row.prop(md, "use_constraint_group", text="Constraints Only")
        col.operator("object.fracture_pack", text="Pack", icon="PACKAGE")
        col.prop(md, "use_particle_birth_coordinates")
        col.prop(md, "percentage")
        sub = col.column(align=True)
        sub.prop_search(md, "thresh_vertex_group", ob, "vertex_groups", text = "Threshold")
        sub.prop_search(md, "passive_vertex_group", ob, "vertex_groups", text = "Passive")
        sub.prop_search(md, "inner_vertex_group", ob, "vertex_groups", text = "Inner")
        sub.prop(md, "inner_crease")
        if (md.frac_algorithm in {'BISECT_FAST', 'BISECT_FAST_FILL', 'BOOLEAN_FRACTAL'}):
            col.prop(md, "orthogonality_factor", text="Rectangular Alignment")

class PHYSICS_PT_fracture_dynamic(PhysicButtonsPanel, Panel):
      bl_label = "Dynamic"
      bl_parent_id = 'PHYSICS_PT_fracture'
      bl_options = {'DEFAULT_CLOSED'}

      def draw_header(self, context):
          md = context.fracture
          self.layout.prop(md, "use_dynamic", text="")

      def draw(self, context):
          md = context.fracture
          layout = self.layout
          layout.active = md.use_dynamic
          row = layout.row(align=True)
          row.prop(md, "dynamic_force")
          row.prop(md, "dynamic_percentage")
          col = layout.column(align=True)
          col.prop(md, "dynamic_new_constraints")
          row = col.row(align=True)
          row.prop(md, "limit_impact")
          row.prop(md, "dynamic_min_size")

class PHYSICS_PT_fracture(PhysicButtonsPanel, Panel):
    bl_label = "Fracture"

    def draw_header_preset(self, context):
        FRACTURE_PT_presets.draw_panel_header(self.layout)

    def draw(self, context):
        md = context.fracture
        layout = self.layout

        layout.context_pointer_set("modifier", md)
        row = layout.row()
        row.operator("object.fracture_refresh", text="Execute Fracture", icon='MOD_EXPLODE').reset = True

class PHYSICS_PT_fracture_basic(PhysicButtonsPanel, Panel):
    bl_label = "Basic"
    bl_parent_id = 'PHYSICS_PT_fracture'
    bl_options = {'DEFAULT_CLOSED'}

    def icon(self, bool):
        if bool:
            return 'TRIA_DOWN'
        else:
            return 'TRIA_RIGHT'

    def draw(self, context):
        layout = self.layout

        md = context.fracture
        ob = context.object

        layout.prop(md, "frac_algorithm")
        if md.frac_algorithm in {'BOOLEAN', 'BOOLEAN_FRACTAL'}:
            col = layout.column(align=True)
            col.prop(md, "boolean_double_threshold")
        col = layout.column(align=True)
        col.prop(md, "shard_count")
        col.prop(md, "point_seed")

        if md.frac_algorithm in {'BOOLEAN', 'BISECT_FILL', 'BISECT_FAST_FILL', 'BOOLEAN_FRACTAL'}:
            col = layout.column()
            col.prop(md, "inner_material")
            col.prop_search(md, "uv_layer", ob.data, "uv_layers", icon="GROUP_UVS")
        if md.frac_algorithm == 'BOOLEAN_FRACTAL':
            col = layout.column(align=True)
            row = col.row(align=True)
            row.prop(md, "fractal_cuts")
            row.prop(md, "fractal_iterations")
            row = col.row(align=True)
            row.prop(md, "fractal_amount")
            row.prop(md, "physics_mesh_scale")
        row = layout.row(align=True)
        row.prop(md, "splinter_axis")
        row = layout.row(align=True)
        row.prop(md, "splinter_length")
        row = layout.row()
        row.prop(md, "split_islands")
        row.prop(md, "use_smooth")
        row = layout.row()
        row.prop(md, "auto_execute")
        row.prop(md, "execute_threaded", text="Threaded (WIP)")

class PHYSICS_PT_fracture_simulation(PhysicButtonsPanel, Panel):
    bl_label = "Constraints"
    bl_options = {'DEFAULT_CLOSED'}
    bl_parent_id = 'PHYSICS_PT_fracture'

    @classmethod
    def poll(cls, context):
        md = context.fracture
        return PhysicButtonsPanel.poll(context)

    def draw_header(self, context):
        md = context.fracture
        self.layout.prop(md, "use_constraints", text="")

    def draw(self, context):
        layout = self.layout
        md = context.fracture
        ob = context.object

        layout.active = md.use_constraints
        layout.label(text="Constraint Building Settings")
        row = layout.row()
        row.prop(md, "use_breaking")
        row.prop(md, "activate_broken")
        row = layout.row()
        row.prop(md, "use_constraint_collision")
        row.prop(md, "use_self_collision")

        col = layout.column(align=True)
        col.prop(md, "constraint_target")
        col.prop(md, "constraint_type")
        col = layout.column(align=True)
        col.prop(md, "constraint_limit", text="Constraint limit, per MeshIsland")
        col.prop(md, "contact_dist")

        layout.label(text="Constraint Cluster Settings")
        layout.prop(md, "cluster_count")
        col = layout.column(align=True)
        col.prop(md, "cluster_group")
        col.prop(md, "cluster_constraint_type", text="Cluster Type")

        if md.use_compounds:
            layout.label(text="Compound Breaking Settings")
        else:
            layout.label(text="Constraint Breaking Settings")

        col = layout.column(align=True)
        col.prop(md, "breaking_threshold", text="Threshold")
        col.prop(md, "cluster_breaking_threshold")

        if md.use_compounds:
            #layout.label("Compound Damage Propagation Settings")
            col = layout.column(align=True)
            col.prop(md, "minimum_impulse")
            #col.prop(md, "impulse_dampening")
            #col.prop(md, "directional_factor")
            col.prop(md, "mass_threshold_factor")
        else:
            layout.label(text="Constraint Special Breaking Settings")
            col = layout.column(align=True)
            row = col.row(align=True)
            row.prop(md, "breaking_percentage", text="Percentage")
            row.prop(md, "cluster_breaking_percentage", text="Cluster Percentage")

            row = col.row(align=True)
            row.prop(md, "breaking_angle", text="Angle")
            row.prop(md, "cluster_breaking_angle", text="Cluster Angle")

            row = col.row(align=True)
            row.prop(md, "breaking_distance", text="Distance")
            row.prop(md, "cluster_breaking_distance", text="Cluster Distance")

            col = layout.column(align=True)
            col.prop(md, "solver_iterations_override")
            col.prop(md, "cluster_solver_iterations_override")

            row = layout.row(align=True)
            row.prop(md, "breaking_angle_weighted")
            row.prop(md, "breaking_distance_weighted")

            row = layout.row(align=True)
            row.prop(md, "breaking_percentage_weighted")
            row.prop(md, "use_mass_dependent_thresholds", text="Mass Dependent Thresholds")

        if not md.use_compounds:
            layout.label(text="Constraint Deform Settings")
            col = layout.column(align=True)
            row = col.row(align=True)
            row.prop(md, "deform_angle", text="Deforming Angle")
            row.prop(md, "cluster_deform_angle", text="Cluster Deforming Angle")

            row = col.row(align=True)
            row.prop(md, "deform_distance", text="Deforming Distance")
            row.prop(md, "cluster_deform_distance", text="Cluster Deforming Distance")

            col.prop(md, "deform_weakening")
            row = layout.row(align=True)
            row.prop(md, "deform_angle_weighted")
            row.prop(md, "deform_distance_weighted")


class PHYSICS_PT_fracture_utilities(PhysicButtonsPanel, Panel):
    bl_label = "Utilities"
    bl_options = {'DEFAULT_CLOSED'}
    bl_parent_id = 'PHYSICS_PT_fracture'

    @classmethod
    def poll(cls, context):
        md = context.fracture
        return PhysicButtonsPanel.poll(context) # and md.fracture_mode != 'EXTERNAL'

    def draw(self, context):
        layout = self.layout
        md = context.fracture
        layout.prop(md, "autohide_filter_group", text = "Filter Group")
        col = layout.column(align=True)
        col.prop(md, "autohide_dist")
        col.prop(md, "automerge_dist")
        row = layout.row()
        row.prop(md, "keep_distort")
        row.prop(md, "do_merge")
        row = layout.row()
        row.prop(md, "use_centroids")
        row.prop(md, "use_vertices")
        row = layout.row()
        row.prop(md, "fix_normals")
        row.prop(md, "nor_range")

        col = layout.column(align=True)
        col.context_pointer_set("modifier", md)
        col.operator("object.rigidbody_convert_to_objects", text = "Convert To Objects", icon="UGLYPACKAGE")
        col.operator("object.rigidbody_convert_to_keyframes", text = "Convert To Keyframed Objects", icon="KEY_HLT")

classes = (
    FRACTURE_PT_presets,
    PHYSICS_PT_fracture,
    PHYSICS_PT_fracture_basic,
    PHYSICS_PT_fracture_advanced,
    PHYSICS_PT_fracture_dynamic,
    PHYSICS_PT_fracture_simulation,
    PHYSICS_PT_fracture_utilities,
    PHYSICS_PT_fracture_anim_mesh,
)

if __name__ == "__main__":  # only for live edit.
    from bpy.utils import register_class
    for cls in classes:
        register_class(cls)
