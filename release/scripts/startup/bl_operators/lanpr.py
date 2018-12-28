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
import string

def lanpr_get_composition_scene(scene):
    n = scene.name+'_lanpr_comp'
    for s in bpy.data.scenes:
        if s.name == n: return s
    return None

def lanpr_make_composition_scene(scene):
    name = scene.name;
    new_name = scene.name+'_lanpr_comp'
    scene.name = new_name
    bpy.ops.scene.new(type='LINK_OBJECTS')
    for s in bpy.data.scenes:
        if s.name == new_name+'.001':
            new_scene = s
            break
    scene.name = name
    new_scene.name = new_name
    
    s = new_scene
    
    s.render.engine = 'BLENDER_LANPR'
    s.use_nodes = True
    
    comp = s.node_tree
    
    comp.nodes.clear()
    n1 = comp.nodes.new("CompositorNodeRLayers")
    n1.scene = scene
    n1.location = (0,0)
    n2 = comp.nodes.new("CompositorNodeRLayers")
    n2.scene = s
    n2.location = (0,-300)
    
    mix = comp.nodes.new("CompositorNodeAlphaOver")
    mix.location = (300,-150)
    comp.links.new(n1.outputs['Image'],mix.inputs[1])
    comp.links.new(n2.outputs['Image'],mix.inputs[2])
    
    out = comp.nodes.new("CompositorNodeComposite")
    out.location = (500,-150)
    comp.links.new(mix.outputs['Image'],out.inputs['Image'])
    
    
class LANPR_make_composition_scene(bpy.types.Operator):
    """Make Composition Scene"""
    bl_idname = "lanpr.make_composition_scene"
    bl_label = "Make Composition Scene"

    @classmethod
    def poll(cls, context):
        return (lanpr_get_composition_scene(context.scene) is None)

    def execute(self, context):
        lanpr_make_composition_scene(context.scene)
        return {'FINISHED'}

def lanpr_remove_composition_scene(scene):
    bpy.data.scenes.remove(lanpr_get_composition_scene(scene))

class LANPR_remove_composition_scene(bpy.types.Operator):
    """Remove Composition Scene"""
    bl_idname = "lanpr.remove_composition_scene"
    bl_label = "Remove Composition Scene"

    @classmethod
    def poll(cls, context):
        return (lanpr_get_composition_scene(context.scene) is not None)

    def execute(self, context):
        lanpr_remove_composition_scene(context.scene)
        return {'FINISHED'}
    
def lanpr_is_composition_scene(scene):
    return scene.name.endswith('_lanpr_comp')

def lanpr_goto_original_scene(scene):
    name = scene.name[:-len('_lanpr_comp')]
    for s in bpy.data.scenes:
        if s.name == name:
            bpy.context.window.scene = s
            break
            
class LANPR_goto_original_scene(bpy.types.Operator):
    """Goto Original Scene"""
    bl_idname = "lanpr.goto_original_scene"
    bl_label = "Goto Original Scene"

    @classmethod
    def poll(cls, context):
        return lanpr_is_composition_scene(context.scene)

    def execute(self, context):
        lanpr_goto_original_scene(context.scene)
        return {'FINISHED'}
    
def lanpr_goto_composition_scene(scene):
    name = scene.name+'_lanpr_comp'
    for s in bpy.data.scenes:
        if s.name == name:
            bpy.context.window.scene = s
            break
            
class LANPR_goto_composition_scene(bpy.types.Operator):
    """Goto Composition Scene"""
    bl_idname = "lanpr.goto_composition_scene"
    bl_label = "Goto Composition Scene"

    @classmethod
    def poll(cls, context):
        return lanpr_get_composition_scene(context.scene) is not None

    def execute(self, context):
        lanpr_goto_composition_scene(context.scene)
        return {'FINISHED'}

@persistent
def lanpr_render_composited_still(scene):

@persistent
def lanpr_render_composited_still(scene):
    
            
class LANPR_render_composited_still(bpy.types.Operator):
    """Render Composited Still"""
    bl_idname = "lanpr.goto_composition_scene"
    bl_label = "Render Composited Still"

    @classmethod
    def poll(cls, context):
        return lanpr_get_composition_scene(context.scene) is not None

    def execute(self, context):
        lanpr_goto_composition_scene(context.scene)
        return {'FINISHED'}

classes=(
    LANPR_make_composition_scene,
    LANPR_remove_composition_scene,
    LANPR_goto_original_scene,
    LANPR_goto_composition_scene,
)