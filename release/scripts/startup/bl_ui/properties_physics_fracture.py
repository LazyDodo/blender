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
        layout.use_property_split = True
        flow = layout.grid_flow(row_major=True, columns=0, even_columns=True, even_rows=False, align=True);
        col = layout.column()
        col.prop(md, "animated_mesh_limit")
        col.prop(md, "use_animated_mesh_rotation")
        col.prop(md, "animated_mesh_input")
        col.operator("object.fracture_anim_bind", text="Bind", icon="UV_VERTEXSEL")

class PHYSICS_PT_fracture_advanced(PhysicButtonsPanel, Panel):
    bl_label = "Advanced"
    bl_options = {'DEFAULT_CLOSED'}
    bl_parent_id = 'PHYSICS_PT_fracture'

    def draw(self, context):
        layout = self.layout
        md = context.fracture
        ob = context.object

        layout.context_pointer_set("modifier", md)
        layout.use_property_split = True
        flow = layout.grid_flow(row_major=True, columns=0, even_columns=True, even_rows=False, align=True);
        col = flow.column()
        col.prop(md, "point_source")
        if 'GRID' in md.point_source:
            col.prop(md, "grid_resolution")
            col.prop(md, "grid_offset")
            col.prop(md, "grid_spacing")
        if 'EXTRA_PARTICLES' in md.point_source or 'EXTRA_VERTS' in md.point_source:
            col.prop(md, "extra_group", text="Helpers")
        if 'CUSTOM' in md.point_source:
            col.prop(md, "cutter_group")
            if (md.cutter_group):
                col.prop(md, "keep_cutter_shards")
                col.prop(md, "material_offset_intersect", text="Intersect Material")
                col.prop(md, "material_offset_difference", text="Difference Material")
        col.prop(md, "use_particle_birth_coordinates")
        col.prop(md, "percentage")

        col.prop_search(md, "thresh_vertex_group", ob, "vertex_groups", text = "Threshold")
        col.prop_search(md, "passive_vertex_group", ob, "vertex_groups", text = "Passive")
        col.prop_search(md, "inner_vertex_group", ob, "vertex_groups", text = "Inner")
        col.prop(md, "inner_crease")

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

          layout.use_property_split = True
          flow = layout.grid_flow(row_major=True, columns=0, even_columns=True, even_rows=False, align=True);
          col = flow.column()
          col.prop(md, "dynamic_shard_count")
          col.prop(md, "dynamic_force")
          col.prop(md, "dynamic_percentage")
          col.prop(md, "dynamic_new_constraints")
          col.prop(md, "dynamic_activation_size")
          col.prop(md, "dynamic_min_size")
          col.prop(md, "limit_impact")

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

        layout.use_property_split = True
        flow = layout.grid_flow(row_major=True, columns=0, even_columns=True, even_rows=False, align=True);
        col = flow.column()

        col.prop(md, "frac_algorithm")
        if md.frac_algorithm in {'BOOLEAN', 'BOOLEAN_FRACTAL'}:
            col.prop(md, "boolean_double_threshold")
        col.prop(md, "shard_count")
        col.prop(md, "point_seed")

        if (md.frac_algorithm in {'BISECT_FAST', 'BISECT_FAST_FILL', 'BOOLEAN_FRACTAL'}):
            col.prop(md, "rectangular_alignment", text="Rectangular Alignment")
        if md.frac_algorithm in {'BOOLEAN', 'BISECT_FILL', 'BISECT_FAST_FILL', 'BOOLEAN_FRACTAL'}:
            col.prop(md, "inner_material")
            col.prop_search(md, "uv_layer", ob.data, "uv_layers", icon="GROUP_UVS")
        if md.frac_algorithm == 'BOOLEAN_FRACTAL':
            col.prop(md, "fractal_cuts")
            col.prop(md, "fractal_iterations")
            col.prop(md, "fractal_amount")
        col.prop(md, "splinter_axis")
        col.prop(md, "splinter_length")
        col.prop(md, "split_islands")
        col.prop(md, "use_smooth")
        col.prop(md, "auto_execute")

class PHYSICS_PT_fracture_constraints(PhysicButtonsPanel, Panel):
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
        layout.use_property_split = True
        flow = layout.grid_flow(row_major=True, columns=0, even_columns=True, even_rows=False, align=True);
        col = flow.column()

        col.prop(md, "use_breaking")
        col.prop(md, "activate_broken")
        col.prop(md, "use_constraint_collision")
        col.prop(md, "use_self_collision")

        col.prop(md, "constraint_type")
        col.prop(md, "constraint_target")
        col.prop(md, "constraint_limit")
        col.prop(md, "contact_dist")
        col.prop(md, "contact_size")

        col.prop(md, "cluster_count")
        col.prop(md, "cluster_group")
        col.prop(md, "cluster_constraint_type", text="Cluster Type")


class PHYSICS_PT_fracture_utilities(PhysicButtonsPanel, Panel):
    bl_label = "Utilities"
    bl_options = {'DEFAULT_CLOSED'}
    bl_parent_id = 'PHYSICS_PT_fracture'

    @classmethod
    def poll(cls, context):
        md = context.fracture
        return PhysicButtonsPanel.poll(context)

    def draw(self, context):
        layout = self.layout
        md = context.fracture

        layout.use_property_split = True
        flow = layout.grid_flow(row_major=True, columns=0, even_columns=True, even_rows=False, align=True);
        col = flow.column()

        col.prop(md, "autohide_filter_group", text = "Filter Group")
        col.prop(md, "autohide_dist")
        col.prop(md, "automerge_dist")
        col.prop(md, "perform_merge")
        col.prop(md, "keep_distort")
        col.prop(md, "use_centroids")
        col.prop(md, "use_vertices")
        col.prop(md, "fix_normals")
        col.prop(md, "normal_search_radius")

        col.context_pointer_set("modifier", md)
        col.operator("object.rigidbody_convert_to_objects", text = "Convert To Objects", icon="UGLYPACKAGE")
        col.operator("object.rigidbody_convert_to_keyframes", text = "Convert To Keyframed Objects", icon="KEY_HLT")

class PHYSICS_PT_fracture_constraints_breaking(PhysicButtonsPanel, Panel):
    bl_label = "Breaking"
    bl_options = {'DEFAULT_CLOSED'}
    bl_parent_id = 'PHYSICS_PT_fracture_constraints'

    @classmethod
    def poll(cls, context):
        md = context.fracture
        return PhysicButtonsPanel.poll(context)

    def draw(self, context):
        layout = self.layout
        md = context.fracture

        layout.active = md.use_constraints
        layout.use_property_split = True
        flow = layout.grid_flow(row_major=True, columns=0, even_columns=True, even_rows=False, align=True);
        col = flow.column()
        col.prop(md, "breaking_threshold", text="Threshold")
        col.prop(md, "cluster_breaking_threshold")
        col.prop(md, "breaking_percentage", text="Percentage")
        col.prop(md, "cluster_breaking_percentage", text="Cluster Percentage")

        col.prop(md, "breaking_angle", text="Angle")
        col.prop(md, "cluster_breaking_angle", text="Cluster Angle")

        col.prop(md, "breaking_distance", text="Distance")
        col.prop(md, "cluster_breaking_distance", text="Cluster Distance")
        col.prop(md, "solver_iterations_override")
        col.prop(md, "cluster_solver_iterations_override")

        col.prop(md, "breaking_angle_weighted")
        col.prop(md, "breaking_distance_weighted")
        col.prop(md, "breaking_percentage_weighted")
        col.prop(md, "use_mass_dependent_thresholds", text="Mass Dependent Thresholds")

class PHYSICS_PT_fracture_constraints_deforming(PhysicButtonsPanel, Panel):
    bl_label = "Deforming"
    bl_options = {'DEFAULT_CLOSED'}
    bl_parent_id = 'PHYSICS_PT_fracture_constraints'

    @classmethod
    def poll(cls, context):
        md = context.fracture
        return PhysicButtonsPanel.poll(context)

    def draw(self, context):
        layout = self.layout
        md = context.fracture

        layout.active = md.use_constraints
        layout.use_property_split = True
        flow = layout.grid_flow(row_major=True, columns=0, even_columns=True, even_rows=False, align=True);
        col = flow.column()

        col.prop(md, "deform_angle", text="Deforming Angle")
        col.prop(md, "cluster_deform_angle", text="Cluster Deforming Angle")

        col.prop(md, "deform_distance", text="Deforming Distance")
        col.prop(md, "cluster_deform_distance", text="Cluster Deforming Distance")

        col.prop(md, "deform_weakening")
        col.prop(md, "deform_angle_weighted")
        col.prop(md, "deform_distance_weighted")

class PHYSICS_PT_fracture_basic_packing(PhysicButtonsPanel, Panel):
    bl_label = "Packing"
    bl_options = {'DEFAULT_CLOSED'}
    bl_parent_id = 'PHYSICS_PT_fracture_basic'

    @classmethod
    def poll(cls, context):
        md = context.fracture
        return PhysicButtonsPanel.poll(context)

    def draw(self, context):
        layout = self.layout
        md = context.fracture

        layout.context_pointer_set("modifier", md)
        layout.use_property_split = True
        flow = layout.grid_flow(row_major=True, columns=0, even_columns=True, even_rows=False, align=True);
        col = flow.column()
        col.prop(md, "pack_group", text="Pack Group")
        col.prop(md, "use_constraint_group", text="Constraints Only")


classes = (
    FRACTURE_PT_presets,
    PHYSICS_PT_fracture,
    PHYSICS_PT_fracture_basic,
    PHYSICS_PT_fracture_basic_packing,
    PHYSICS_PT_fracture_advanced,
    PHYSICS_PT_fracture_utilities,
    PHYSICS_PT_fracture_dynamic,
    PHYSICS_PT_fracture_constraints,
    PHYSICS_PT_fracture_anim_mesh,
    PHYSICS_PT_fracture_constraints_breaking,
    PHYSICS_PT_fracture_constraints_deforming,
)

if __name__ == "__main__":  # only for live edit.
    from bpy.utils import register_class
    for cls in classes:
        register_class(cls)
