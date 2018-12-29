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

GLOBAL_OUTPUT_PATH = ''

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

def lanpr_get_original_scene(scene):
    name = scene.name[:-len('_lanpr_comp')]
    for s in bpy.data.scenes:
        if s.name == name: return s
    return None

def lanpr_goto_original_scene(scene):
    s = lanpr_get_original_scene(scene)
    if s: bpy.context.window.scene = s
    
            
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


#callbacks

GC = None

def lanpr_render_next_frame(sc):
    global GLOBAL_OUTPUT_PATH
    sc.frame_current = sc.frame_current+1
    if sc.frame_current>sc.frame_end: 
        bpy.app.handlers.render_complete.remove(lanpr_render_next_frame)
        bpy.context.scene.render.filepath = GLOBAL_OUTPUT_PATH
        return
    
    bpy.app.handlers.render_cancel.append(lanpr_render_canceled)
    bpy.app.handlers.render_complete.remove(lanpr_render_next_frame)
    
    lanpr_render_backdrop_first(sc)

def lanpr_render_this_scene_next(scene):
    
    bpy.app.handlers.render_complete.remove(lanpr_render_this_scene_next)
    bpy.app.handlers.render_cancel.remove(lanpr_render_canceled)
    
    bpy.app.handlers.render_cancel.append(lanpr_render_canceled)
    
    sc = lanpr_get_composition_scene(scene)
    write = sc.lanpr.composite_render_animation
    
    bpy.context.scene.render.filepath = GLOBAL_OUTPUT_PATH + '/%04d'%sc.frame_current + bpy.context.scene.render.file_extension
    
    if sc.lanpr.composite_render_animation:
        bpy.app.handlers.render_complete.append(lanpr_render_next_frame)
        global GC
        bpy.ops.render.render(scene=sc.name, write_still = write)
    else:
        #'INVOKE_DEAFULT' still cause trouble on windows.
        #bpy.ops.render.render(GC,'INVOKE_DEFAULT',scene=sc.name)
        bpy.ops.render.render(scene=sc.name)
    
def lanpr_render_canceled(scene):
    
    bpy.app.handlers.render_complete.remove(lanpr_render_this_scene_next)
    
    bpy.app.handlers.render_cancel.remove(lanpr_render_canceled)

def lanpr_render_backdrop_first(this_scene):
    s = lanpr_get_original_scene(this_scene)
    if not s: return

    s.frame_current = this_scene.frame_current
    bpy.app.handlers.render_complete.append(lanpr_render_this_scene_next)
    bpy.ops.render.render(scene=s.name)
            
class LANPR_render_composited(bpy.types.Operator):
    """Render Composited"""
    bl_idname = "lanpr.render_composited"
    bl_label = "Render Composited"

    @classmethod
    def poll(cls, context):
        return True
    
    def execute(self, context):
        if bpy.context.scene.lanpr.composite_render_animation:
            s = lanpr_get_original_scene(bpy.context.scene)
            bpy.context.scene.frame_current = bpy.context.scene.frame_start
            s.frame_current = bpy.context.scene.frame_start
            bpy.context.scene.frame_end = s.frame_end
            bpy.context.scene.render.filepath = s.render.filepath
            
        global GLOBAL_OUTPUT_PATH
        GLOBAL_OUTPUT_PATH = bpy.context.scene.render.filepath

        bpy.app.handlers.render_cancel.append(lanpr_render_canceled)
        global GC
        GC = bpy.context.copy()
        
        lanpr_render_backdrop_first(bpy.context.scene)
        
        return {'FINISHED'}

classes=(
    LANPR_make_composition_scene,
    LANPR_remove_composition_scene,
    LANPR_goto_original_scene,
    LANPR_goto_composition_scene,
    LANPR_render_composited,
)