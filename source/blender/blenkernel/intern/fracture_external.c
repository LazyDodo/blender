#include <stdio.h>

#include "MEM_guardedalloc.h"

#include "BKE_fracture.h"
#include "BKE_mesh.h"
#include "BKE_deform.h"
#include "BKE_material.h"
#include "BKE_main.h"
#include "BKE_modifier.h"
#include "BKE_rigidbody.h"
#include "BKE_pointcache.h"
#include "BKE_object.h"

#include "BLI_listbase.h"
#include "BLI_ghash.h"
#include "BLI_math.h"
#include "BLI_math_vector.h"
#include "BLI_string.h"

#include "DNA_object_types.h"
#include "DNA_modifier_types.h"
#include "DNA_mesh_types.h"
#include "DNA_meshdata_types.h"
#include "DNA_fracture_types.h"
#include "DNA_material_types.h"
#include "DNA_rigidbody_types.h"
#include "DNA_scene_types.h"

static MeshIsland* fracture_shard_to_island(Scene *scene, FractureModifierData *fmd, Shard *s, int vertstart, float quat[4]);
static Shard* fracture_object_to_shard( Object *own, Object* target);
static void fracture_update_shards(FractureModifierData *fmd, Shard *s);

int BKE_fracture_update_visual_mesh(FractureModifierData *fmd, Object *ob, bool do_custom_data)
{
    MeshIsland *mi;
    Mesh *dm = fmd->visible_mesh_cached;
    int vertstart = 0, totvert = 0, totpoly = 0, polystart = 0, matstart = 1, defstart = 0, loopstart = 0;
    MVert *mv = NULL;
    MPoly *mp = NULL, *mpoly = NULL, *ppoly = NULL, *pp = NULL, *spoly = NULL, *sp = NULL, *tpoly = NULL, *tp = NULL;
    int i = 0, j = 0;
    MDeformVert *dvert = NULL;
    Mesh *me = (Mesh*)ob->data;
    Shard *s, *t;

    if (dm)
    {
        vertstart = dm->totvert;
        BKE_mesh_free(dm);
        dm = fmd->visible_mesh_cached = NULL;
    }

    fmd->visible_mesh_cached = BKE_fracture_assemble_mesh_from_shards(fmd, do_custom_data, fmd->pack_storage.first != NULL);

    if (!fmd->visible_mesh_cached)
        return 0;

    //store start mesh in order to be able to change autohide dist based on it later in sim too !
    if (fmd->dm) {
        BKE_mesh_free(fmd->dm);
        fmd->dm = NULL;
    }

    BKE_mesh_nomain_to_mesh(fmd->visible_mesh_cached, fmd->dm, ob, CD_MASK_MESH, false);

    dm = fmd->visible_mesh_cached;
    mv = dm->mvert;
    totvert = dm->totvert;
    mpoly = dm->mpoly;
    dvert = CustomData_get_layer(&me->vdata, CD_MDEFORMVERT);

    CustomData_merge(&dm->ldata, &me->ldata, CD_MASK_MLOOPUV, CD_CALLOC, dm->totloop);

    s = fmd->frac_mesh->shard_map.first;
    t = fmd->pack_storage.first;

    //update existing island's vert refs, if any...should have used indexes instead :S
    for (mi = fmd->meshIslands.first; mi; mi = mi->next)
    {
        MVert *pvert = mi->physics_mesh->mvert;
        float inv_size[3] = {1.0f, 1.0f, 1.0f};

        inv_size[0] = 1.0f / s->impact_size[0];
        inv_size[1] = 1.0f / s->impact_size[1];
        inv_size[2] = 1.0f / s->impact_size[2];

        for (i = 0; i < mi->vertex_count; i++)
        {
            //just update pointers, dont need to reallocate something
            MVert *v = NULL, *pv = NULL;
            int index;

            //also correct indexes
            if (mi->vertex_indices[i] >= totvert)
            {
                index = mi->vertex_indices[i];
                mi->vertex_indices[i] -= (vertstart - totvert);
                printf("I: %d, O: %d, N: %d\n", i, index, mi->vertex_indices[i]);
            }

            index = mi->vertex_indices[i];
            v = mv + index;
            mi->vertices_cached[i] = v;

            pv = pvert + i;
            mul_v3_v3(pv->co, inv_size);

            //transform vertex properly ? compensate for shrunken shard ?
            //sub_v3_v3v3(loc, mi->centroid, s->raw_centroid);
            //loc_quat_size_to_mat4(mat, loc , rot, s->impact_size);

            //invert_m4_m4(imat, mat);
#if 0
            pv = pvert + i;
            mul_m4_v3(imat, pv->co);
#endif

            //eliminate shrink but take also difference in centroids into account here
            //sub_v3_v3(v->co, s->centroid);
            //mul_m4_v3(imat, v->co);
            //add_v3_v3(v->co, s->centroid);

//			sub_v3_v3(v->co, s->centroid);
//			mul_v3_v3(v->co, inv_size);
//			add_v3_v3(v->co, s->centroid);

            //printf("%d %d\n", index, dm->getNumVerts(dm));

            //hrm perhaps we need to update rest coordinates, too...
            mi->vertco[3*i] = v->co[0];
            mi->vertco[3*i+1] = v->co[1];
            mi->vertco[3*i+2] = v->co[2];

            mi->vertno[3*i] = v->no[0];
            mi->vertno[3*i+1] = v->no[1];
            mi->vertno[3*i+2] = v->no[2];

            if (ob && dvert)
            {
                int l;
                MDeformVert *dv = dvert + mi->vertex_indices[i];
                if (dv && dv->dw)
                {
                    for (l = 0; l < dv->totweight; l++)
                    {
                        MDeformWeight *dw = dv->dw;
                        //refill mapping data, to make it accessible for each vert (for dumb mapping function)
                        int key = defstart + l;
                        int index = GET_INT_FROM_POINTER(BLI_ghash_lookup(fmd->defgrp_index_map, SET_INT_IN_POINTER(key)));
                        //printf("Got: %d %d\n", key, index);
                        if (dw->def_nr == l)
                            dw->def_nr = index;
                    }

                    //XXX TODO store this on physics mesh too ? to be able to reload it from blend
                }
            }
        }

        defstart += mi->totdef;

        totpoly = mi->physics_mesh->totpoly;
        ppoly = mi->physics_mesh->mpoly;
        spoly = s->mpoly;
        if (t) {
            tpoly = t->mpoly;
        }

        for (j = 0, mp = mpoly + polystart, pp = ppoly, sp = spoly; j < totpoly; j++, mp++, pp++, sp++)
        {
            /* material index lookup and correction, avoid having the same material in different slots */
            int index = GET_INT_FROM_POINTER(BLI_ghash_lookup(fmd->material_index_map,
                                             SET_INT_IN_POINTER(mp->mat_nr + matstart)));

            if (index > 0)
                index--;

            mp->mat_nr = index;
            //store this on physics mesh as well, and on shard too so for being able to reload it from blend later (without
            // having a materialmap then)
            pp->mat_nr = index;
            sp->mat_nr = index;

            //also dont forget pack storage
            if (tpoly) {
                tp = tpoly + j;
                tp->mat_nr = index;
            }
        }

        /* fortunately we know how many faces "belong" to this meshisland, too */
        polystart += totpoly;
        matstart += mi->totcol;
        loopstart += s->totloop;

        if (s) {
            s = s->next;
        }

        if (t) {
            t = t->next;
        }
    }

    return vertstart;
}

int fracture_collect_defgrp(Object* o, Object* ob, int defstart, GHash** def_index_map)
{
    bDeformGroup *vgroup, *ngroup;
    int k = 0;

    /* create vertexgroups on new object, if they dont exist already there*/
    for (vgroup = o->defbase.first; vgroup; vgroup = vgroup->next) {
        int index = defgroup_name_index(ob, vgroup->name);
        int key = defstart + k;

        if (index == -1) {
            // old group index + defstart to make it somehow linearized
            ngroup = MEM_callocN(sizeof(bDeformGroup), "collect deformGroup");
            memcpy(ngroup, vgroup, sizeof(bDeformGroup));
            BLI_addtail(&ob->defbase, ngroup);
            index = BLI_listbase_count(&ob->defbase)-1;
        }

        if (!BLI_ghash_haskey(*def_index_map, SET_INT_IN_POINTER(key)))
            BLI_ghash_insert(*def_index_map, SET_INT_IN_POINTER(key), SET_INT_IN_POINTER(index));

        k++;
    }

    if (ob->defbase.first && ob->actdef == 0)
        ob->actdef = 1;

    return k;
}

short BKE_fracture_collect_materials(Main* bmain, Object* o, Object* ob, int matstart, GHash** mat_index_map)
{
    short *totcolp = NULL;
    Material ***matarar = NULL;
    int j;

    /* append materials to target object, if not existing yet */
    totcolp = give_totcolp(o);
    matarar = give_matarar(o);

    for (j = 0; j < (*totcolp); j++)
    {
        void *key;
        int index = BKE_object_material_slot_find_index(ob, (*matarar)[j]);
        if (index == 0)
        {
            index = ob->totcol+1;
            assign_material(bmain, ob, (*matarar)[j], index, BKE_MAT_ASSIGN_USERPREF);
        }

        key = SET_INT_IN_POINTER(matstart+j);
        if (!BLI_ghash_haskey(*mat_index_map, key))
            BLI_ghash_insert(*mat_index_map, key, SET_INT_IN_POINTER(index));
    }

    return (*totcolp);
}

static void pack_storage_add(FractureModifierData *fmd, Shard* s)
{
    Shard *t = BKE_fracture_shard_copy(s);
    BLI_addtail(&fmd->pack_storage, t);
}

static void fracture_collect_layer(CustomData* src, CustomData *dst, int totelem, int cd_type, int dst_offset, int count)
{
    int layerstart = CustomData_get_layer_index(src, cd_type);
    int totlayer = CustomData_number_of_layers(src, cd_type);
    int j;

    for (j = 0; j < totlayer; j++)
    {
        const char *name = CustomData_get_layer_name(src, cd_type, j);

        //find index of named layer in dst mesh
        int index = CustomData_get_named_layer_index(dst, cd_type, name);
        if (index == -1)
        {
            //add layer if not there
            //void *layer = CustomData_get_layer_named(src, cd_type, name);
            CustomData_add_layer_named(dst, cd_type, CD_CALLOC, NULL, totelem, name);
        }

        index = CustomData_get_named_layer_index(dst, cd_type, name);
        CustomData_copy_data_layer(src, dst, j+layerstart, index, 0, dst_offset, count);
    }
}

void BKE_fracture_collect_layers(Shard* s, Mesh *dm, int vertstart, int polystart, int loopstart, int edgestart)
{
    int totloop = dm->totloop;
    int totvert = dm->totvert;
    int totedge = dm->totedge;

    fracture_collect_layer(&s->vertData, &dm->vdata, totvert, CD_MDEFORMVERT, vertstart, s->totvert);
    fracture_collect_layer(&s->loopData, &dm->ldata, totloop, CD_MLOOPUV, loopstart, s->totloop);
    fracture_collect_layer(&s->edgeData, &dm->edata, totedge, CD_CREASE, edgestart, s->totedge);
    fracture_collect_layer(&s->edgeData, &dm->edata, totedge, CD_BWEIGHT, edgestart, s->totedge);
    fracture_collect_layer(&s->vertData, &dm->vdata, totvert, CD_PROP_FLT, vertstart, s->totvert);
}

MeshIsland* BKE_fracture_mesh_island_add(Main* bmain, FractureModifierData *fmd, Object* own, Object *target, Scene *scene)
{
    MeshIsland *mi;
    Shard *s;
    int vertstart = 0;
    short totcol = 0, totdef = 0;
    float loc[3], quat[4], iquat[4];

    if (fmd->fracture_mode != MOD_FRACTURE_EXTERNAL || own->type != OB_MESH || !own->data)
        return NULL;

    if (target->type != OB_MESH || !target->data)
        return NULL;

    //lets see whether we need to add loc here too XXX TODO
    mat4_to_loc_quat(loc, quat, target->obmat);

    s = fracture_object_to_shard(own, target);
    copy_v3_v3(s->centroid, loc);

    fracture_update_shards(fmd, s);

    vertstart = fmd->frac_mesh->progress_counter;
    fmd->frac_mesh->progress_counter += s->totvert;

    //hrm need to rebuild ALL islands since vertex refs are bonkers now after mesh has changed
    invert_qt_qt(iquat, quat);
    mi = fracture_shard_to_island(scene, fmd, s, vertstart, iquat);

    copy_qt_qt(mi->rot, quat);
    copy_v3_v3(mi->centroid, loc);

    mi->rigidbody = BKE_rigidbody_create_shard(bmain, scene, own, target, mi);
    if (mi->rigidbody)
    {
        mi->rigidbody->meshisland_index = mi->id;
    }

    BLI_strncpy(mi->name, target->id.name + 2, MAX_ID_NAME - 2);

    //handle materials
    if (!fmd->material_index_map)
    {
        fmd->material_index_map = BLI_ghash_int_new("mat_index_map");
        fmd->matstart = 1;
    }

    totcol = BKE_fracture_collect_materials(bmain, target, own, fmd->matstart, &fmd->material_index_map);
    if (totcol < 0)
        totcol = 0;
    fmd->matstart += totcol;
    mi->totcol = totcol;

    /*XXXXX TODO deal with material deletion, and reorder (in material code) */

    //handle vertexgroups, too
    if (!fmd->defgrp_index_map)
    {
        fmd->defgrp_index_map = BLI_ghash_int_new("defgrp_index_map");
        fmd->defstart = 0;
    }

    totdef = fracture_collect_defgrp(target, own, fmd->defstart, &fmd->defgrp_index_map);
    if (totdef < 0)
        totdef = 0;
    fmd->defstart += totdef;
    mi->totdef = totdef;

    //XXX TODO handle UVs, shapekeys and more ?
//	fracture_collect_uv_tex(target, own);

    //add shard to pack storage
    pack_storage_add(fmd, s);

    return mi;
}

void BKE_fracture_mesh_island_free(FractureModifierData *rmd, MeshIsland *mi, bool remove_rigidbody, Scene* scene)
{
    if (mi->physics_mesh) {
        BKE_mesh_free(mi->physics_mesh);
        mi->physics_mesh = NULL;
    }

    if (mi->rigidbody) {
        if (remove_rigidbody)
            BKE_rigidbody_remove_shard(scene, mi);
        MEM_freeN(mi->rigidbody);
        mi->rigidbody = NULL;
    }

    if (mi->vertco) {
        MEM_freeN(mi->vertco);
        mi->vertco = NULL;
    }

    if (mi->vertno) {
        MEM_freeN(mi->vertno);
        mi->vertno = NULL;
    }

    if (mi->vertices) {
        //MEM_freeN(mi->vertices);
        mi->vertices = NULL; /*borrowed only !!!*/
    }

    if (mi->vertices_cached) {
        MEM_freeN(mi->vertices_cached);
        mi->vertices_cached = NULL;
    }

    if (mi->bb != NULL) {
        MEM_freeN(mi->bb);
        mi->bb = NULL;
    }

    if (mi->participating_constraints != NULL) {
        int i = 0;
        for (i = 0; i < mi->participating_constraint_count; i++)
        {
            RigidBodyShardCon *con = mi->participating_constraints[i];
            if (con) {
                con->mi1 = NULL;
                con->mi2 = NULL;
            }
        }

        MEM_freeN(mi->participating_constraints);
        mi->participating_constraints = NULL;
        mi->participating_constraint_count = 0;
    }

    if (mi->vertex_indices) {
        MEM_freeN(mi->vertex_indices);
        mi->vertex_indices = NULL;
    }

    if (mi->rots) {
        MEM_freeN(mi->rots);
        mi->rots = NULL;
    }

    if (mi->locs) {
        MEM_freeN(mi->locs);
        mi->locs = NULL;
    }

    if (mi->acc_sequence)
    {
        MEM_freeN(mi->acc_sequence);
        mi->acc_sequence = NULL;
    }

    mi->frame_count = 0;

    MEM_freeN(mi);
    mi = NULL;
}

void pack_storage_remove(FractureModifierData *fmd, Shard *s)
{
    Shard *t = fmd->pack_storage.first;
    while(t)
    {
        if (t->shard_id == s->shard_id)
        {
            BLI_remlink(&fmd->pack_storage, t);
            BKE_fracture_shard_free(t, true);
            break;
        }

        t = t->next;
    }
}

void BKE_fracture_mesh_island_remove(FractureModifierData *fmd, MeshIsland *mi, Scene* scene)
{
    if (BLI_listbase_is_single(&fmd->meshIslands))
    {
        BKE_fracture_mesh_island_remove_all(fmd, scene);
        return;
    }

    if (fmd->frac_mesh && mi)
    {
        int index = BLI_findindex(&fmd->meshIslands, mi);
        Shard *s = BLI_findlink(&fmd->frac_mesh->shard_map, index);
        if (s)
        {
            int i;
            BLI_remlink(&fmd->frac_mesh->shard_map, s);
            pack_storage_remove(fmd, s);
            BKE_fracture_shard_free(s, true);
            fmd->frac_mesh->shard_count--;

            BLI_remlink(&fmd->meshIslands, mi);
            for (i = 0; i < mi->participating_constraint_count; i++)
            {
                RigidBodyShardCon *con = mi->participating_constraints[i];
                BLI_remlink(&fmd->meshConstraints, con);
                BKE_rigidbody_remove_shard_con(scene, con);
            }

            BKE_fracture_mesh_island_free(fmd, mi, true, scene);
        }
    }
}

static void pack_storage_remove_all(FractureModifierData*fmd)
{
    Shard *s;
    while (fmd->pack_storage.first) {
        s = fmd->pack_storage.first;
        BLI_remlink(&fmd->pack_storage, s);
        BKE_fracture_shard_free(s, true);
    }
}

void BKE_fracture_mesh_island_remove_all(FractureModifierData *fmd, Scene* scene)
{
    MeshIsland *mi;

    //free all shards
    BKE_fracture_fracmesh_free(fmd->frac_mesh, true);
    MEM_freeN(fmd->frac_mesh);
    fmd->frac_mesh = NULL;

    //free pack storage
    pack_storage_remove_all(fmd);

    //free all constraints first
    BKE_fracture_constraints_free(fmd, scene);

    //free all meshislands
    while (fmd->meshIslands.first) {
        mi = fmd->meshIslands.first;
        BLI_remlink(&fmd->meshIslands, mi);
        BKE_fracture_mesh_island_free(fmd, mi, true, scene);
    }

    fmd->meshIslands.first = NULL;
    fmd->meshIslands.last = NULL;

    //free visual_mesh
    if (fmd->visible_mesh_cached)
    {
        BKE_mesh_free(fmd->visible_mesh_cached);
        fmd->visible_mesh_cached = NULL;
    }

    if (fmd->material_index_map)
    {
        BLI_ghash_free(fmd->material_index_map, NULL, NULL);
        fmd->material_index_map = NULL;
        fmd->matstart = 1;
    }

    if (fmd->defgrp_index_map)
    {
        BLI_ghash_free(fmd->defgrp_index_map, NULL, NULL);
        fmd->defgrp_index_map = NULL;
        fmd->defstart = 0;
    }
}

Mesh* BKE_fracture_external_apply(FractureModifierData *fmd, Object* ob, Mesh* pack_dm, Mesh* derivedData, Scene *scene)
{
    Mesh *final_dm;

    if (ob->type != OB_MESH)
    {	//sanity check
        if (pack_dm != derivedData)
        {
            BKE_mesh_free(pack_dm);
            pack_dm = NULL;
        }

        return derivedData;
    }

    fmd->refresh = false;
    fmd->shards_to_islands = false;

    if (!fmd->visible_mesh_cached)
    {
        BKE_fracture_update_visual_mesh(fmd, ob, true);
        fmd->refresh_autohide = true;

        if (fmd->face_pairs != NULL) {
            BLI_ghash_free(fmd->face_pairs, NULL, NULL);
            fmd->face_pairs = NULL;
        }
    }

    if (fmd->visible_mesh_cached && fmd->dm) {

        if (fmd->refresh_autohide) {

            if (fmd->autohide_dist > 0) {

                if (!fmd->face_pairs)
                {
                    fmd->face_pairs = BLI_ghash_int_new("face_pairs");
                }

                BKE_fracture_face_pairs(fmd, fmd->dm, ob);

                if (!fmd->distortion_cached)
                {
                    BKE_fracture_shared_verts_free(&fmd->shared_verts);
                    BKE_fracture_shared_vert_groups(fmd, fmd->dm, &fmd->shared_verts);
                }
            }

            fmd->refresh_autohide = false;
        }

        if (fmd->autohide_dist > 0 || fmd->automerge_dist > 0)
        {
            final_dm = BKE_fracture_autohide_do(fmd, fmd->visible_mesh_cached, ob, scene);
        }
        else {
            BKE_mesh_nomain_to_mesh(fmd->visible_mesh_cached, final_dm, ob, CD_MASK_MESH, false);
            if (!fmd->fix_normals) {
                BKE_mesh_calc_normals(final_dm);
            }
        }
    }

    return final_dm;
}

static Shard* fracture_object_to_shard( Object *own, Object* target)
{
    Mesh *dm;
    Shard *s = NULL;
    Object *ob_eval = NULL; //get eval object from target

    MVert* mvert, *mv;
    MPoly* mpoly;
    MLoop* mloop;
    MEdge* medge;
    SpaceTransform trans;
    float mat[4][4], size[3] = {1.0f, 1.0f, 1.0f};

    int totvert, totpoly, totloop, totedge, v;
    bool do_free = false;

    //dm = target->derivedFinal; //eval ob, modifier mesh (TODO)
    dm = BKE_modifier_get_evaluated_mesh_from_evaluated_object(ob_eval, &do_free);

    if (!dm)
    {	//fallback if no derivedFinal available
        dm = target->data;
    }

    unit_m4(mat);
    BLI_space_transform_from_matrices(&trans, target->obmat, mat);
    //BLI_SPACE_TRANSFORM_SETUP(&trans, target, own);
    mat4_to_size(size, target->obmat);

    mvert = dm->mvert;
    mpoly = dm->mpoly;
    mloop = dm->mloop;
    medge = dm->medge;
    totvert = dm->totvert;
    totpoly = dm->totpoly;
    totloop = dm->totloop;
    totedge = dm->totedge;

    // create temp shard -> that necessary at all ?
    s = BKE_fracture_shard_create(mvert, mpoly, mloop, medge, totvert, totpoly, totloop, totedge, true);

    //use this as size holder, and rawcentroid is the old ob location
    copy_v3_v3(s->impact_size, size);

    //compare centroid in worldspace with location
    mul_v3_m4v3(s->raw_centroid, target->obmat, s->centroid);

    for (v = 0, mv = s->mvert; v < s->totvert; v++, mv++)
    {
        //shrink the shard ? (and take centroid diff into account here, too)
        BLI_space_transform_apply(&trans, mv->co);
    }

    BLI_space_transform_apply(&trans, s->centroid);

    BKE_fracture_custom_data_mesh_to_shard(s, dm);
    BKE_shard_calc_minmax(s);

    if (do_free && dm)
    {
        BKE_mesh_free(dm);
    }

    return s;
}

static void fracture_update_shards(FractureModifierData *fmd, Shard *s)
{
    FracMesh* fm;

    if (!fmd->frac_mesh)
    {
        fmd->frac_mesh = BKE_fracture_fracmesh_create();
        fmd->frac_mesh->progress_counter = 0; //XXXX ABUSE this for vertstart now, threading doesnt work anyway yet
        fmd->matstart = 1; //TODO, is this 1-based ?
    }

    fm = fmd->frac_mesh;
    BLI_addtail(&fm->shard_map, s);
    s->shard_id = fm->shard_count;
    fm->shard_count++;
}

static MeshIsland* fracture_shard_to_island(Scene* scene, FractureModifierData *fmd, Shard *s, int vertstart, float quat[4])
{
    MeshIsland *mi;
    int k = 0, j = 0, totvert;
    MVert *mverts = NULL, *verts, *mv;

    //create mesh island and intialize
    mi = MEM_callocN(sizeof(MeshIsland), "meshIsland");
    BLI_addtail(&fmd->meshIslands, mi);
    mi->participating_constraints = NULL;
    mi->participating_constraint_count = 0;
    mi->thresh_weight = 0.0f;
    mi->ground_weight = 0.0f;
    mi->vertex_count = s->totvert;
    mi->totcol = 0;
    mi->totdef = 0;

    //link up the visual mesh verts
    mi->vertices_cached = MEM_mallocN(sizeof(MVert *) * s->totvert, "vert_cache");
    if (fmd->visible_mesh_cached) /*ensure to be NULL in "pack, unpack" methods */
        mverts = fmd->visible_mesh_cached->mvert;
    mi->vertex_indices = MEM_mallocN(sizeof(int) * mi->vertex_count, "mi->vertex_indices");

    for (k = 0; k < s->totvert; k++) {
        if (mverts)
        {
            mi->vertices_cached[k] = mverts + vertstart + k;
        }
        else
        {
            mi->vertices_cached[k] = NULL;
        }
        mi->vertex_indices[k] = vertstart + k;
    }

    //some dummy setup, necessary here ?
    mi->locs = MEM_mallocN(sizeof(float)*3, "mi->locs");
    mi->rots = MEM_mallocN(sizeof(float)*4, "mi->rots");
    mi->frame_count = 0;
    if (scene && scene->rigidbody_world)
    {
        /*modifier might have no linked scene yet after creation on an inactive layer */
        /*so just try a fallback here */
        mi->start_frame = scene->rigidbody_world->shared->pointcache->startframe;
    }
    else
    {
        mi->start_frame = 1;
    }

    mi->physics_mesh = BKE_fracture_shard_to_mesh(s, true);
    totvert = mi->physics_mesh->totvert;
    verts = mi->physics_mesh->mvert;

    mi->vertco = MEM_mallocN(sizeof(float) * 3 * totvert, "vertco");
    mi->vertno = MEM_mallocN(sizeof(short) * 3 * totvert, "vertno");

    for (mv = verts, j = 0; j < totvert; mv++, j++) {
        short no[3];

        mi->vertco[j * 3] = mv->co[0];
        mi->vertco[j * 3 + 1] = mv->co[1];
        mi->vertco[j * 3 + 2] = mv->co[2];

        copy_v3_v3_short(no, mv->no);

        mi->vertno[j * 3] = no[0];
        mi->vertno[j * 3 + 1] = no[1];
        mi->vertno[j * 3 + 2] = no[2];

        /* then eliminate centroid in vertex coords*/
        sub_v3_v3(mv->co, s->centroid);

        mul_qt_v3(quat, mv->co);
    }

    copy_v3_v3(mi->centroid, s->centroid);
    mi->id = s->shard_id;
    mi->bb = BKE_boundbox_alloc_unit();
    BKE_boundbox_init_from_minmax(mi->bb, s->min, s->max);
    mi->particle_index = -1;

    //this info isnt necessary here... constraints will be provided too !
    mi->neighbor_ids = s->neighbor_ids;
    mi->neighbor_count = s->neighbor_count;

    return mi;
}

