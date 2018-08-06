#include "MEM_guardedalloc.h"
#include "BKE_fracture.h"

#include "DNA_fracture_types.h"
#include "DNA_object_types.h"
#include "DNA_mesh_types.h"
#include "DNA_meshdata_types.h"
#include "DNA_modifier_types.h"
#include "DNA_scene_types.h"

Mesh *BKE_fracture_dynamic_apply(FractureModifierData *fmd, Object *ob, Mesh *derivedData)
{
#if 0
    char (*names)[66] = NULL;
    bool valid_fractured_mesh = false;

    if (fmd->refresh_constraints)
    {
        BKE_free_constraints(fmd);
    }

    /* group_dm, clean_dm not necessary here as we dont support non-mesh objects and subobject_groups here */
    // this should be done from an operator, really...
    if (fmd->refresh)
    {
        BKE_fracture_dynamic_initialize(fmd, ob, derivedData, names);

        fmd->current_mi_entry->is_new = false;
        fmd->refresh = false;
        fmd->distortion_cached = false;

        //constraint + autohide refresh
       // BKE_fracture_refresh_all_constraints(scene, ob, fmd);

        BKE_fracture_autohide_refresh(fmd, ob);

        if (fmd->autohide_dist > 0 && !fmd->distortion_cached) {
            BKE_fracture_automerge_refresh(fmd);
        }

        if (names) {
            MEM_freeN(names);
        }
    }
    else {
        //here we should just read out our mi->rigidbodies loc / rot and redraw with BKE_update_cell
        //BKE_fracture_dynamic_update(fmd, scene);
    }

    /*XXX better rename this, it checks whether we have a valid fractured mesh */
    valid_fractured_mesh = !fmd->explo_shared || (fmd->explo_shared && fmd->dm && fmd->frac_mesh);

   /* if ((!valid_fractured_mesh) || (fmd->visible_mesh == NULL && fmd->visible_mesh_cached == NULL)) {
        do_clear(fmd);
    }*/

    return BKE_fracture_mesh_result(fmd, derivedData, ob, valid_fractured_mesh);
#endif
    return derivedData;
}

#if 0
int BKE_fracture_dynamic_update(FractureModifierData *fmd, Object* ob, Scene* scene, Mesh* dm)
{
    int frame = (int)BKE_scene_frame_get(scene);
    /*HERE we must know which shard(s) to fracture... hmm shards... we should "merge" states which happen in the same frame automatically !*/
    /* TODO_2 dynamic, this is called from rigidbody system only !!! so move out of the md loop as well, to BKE */
    if (!(BKE_lookup_mesh_state(fmd, frame, false)))
    {
        /*simulation mode*/
        /* bullet callbacks may happen multiple times per frame, in next frame we can evaluate them all,
         * so we need some array of shardIDs or shards to fracture each *
         * we need to loop over those shard IDs here, but lookup of shard ids might be slow, but fracturing of many shards is slower...
         * should not have a visible effect in general */

        int count = 0;
        int totpoint = 0;
        FractureID *fid = fmd->fracture_ids.first;

        while(fid) {
            FracPointCloud points;
            points = BKE_fracture_points(fmd, ob, dm, fid->shardID);
            totpoint += points.totpoints;
            MEM_freeN(points.points);
            points.totpoints = 0;
            fid = fid->next;
        }

        if (totpoint == 0)
        {
            fmd->update_dynamic = false;
        }
        else
        {
            if (fmd->update_dynamic)
            {
                if (!fmd->is_dynamic_external)
                {
                    BKE_free_constraints(fmd);
                }

                printf("ADD NEW 2: %s \n", ob->id.name);
                fmd->update_dynamic = false;
                BKE_fracture_dynamic_new_entries_add(fmd, dm, ob);
            }

            while(fmd->fracture_ids.first){
                fid = (FractureID*)fmd->fracture_ids.first;
                do_fracture(fmd, fid->shardID, ob, dm);
                BLI_remlink(&fmd->fracture_ids, fid);
                MEM_freeN(fid);
                count++;
            }

            if (count > 0)
            {
                //BKE_free_constraints(fmd);
                printf("REFRESH: %s \n", ob->id.name);
                scene->rigidbody_world->flag |= RBW_FLAG_OBJECT_CHANGED;
                scene->rigidbody_world->flag |= RBW_FLAG_REFRESH_MODIFIERS;
            }
        }
    }

    fmd->last_frame = frame;
}

void BKE_fracture_dynamic_initialize(FractureModifierData *fmd, Object *ob, Mesh *dm, char (**names)[66])
{
    /*TODO_1 refresh, move to BKE and just call from operator for prefractured case*/
    int mi_count = 0;

    printf("ADD NEW 1: %s \n", ob->id.name);
    if ((fmd->last_frame == INT_MAX)) //TODO event if we restart at cache startframe
    {
        if (fmd->reset_shards)
        {
            free_simulation(fmd, true, true);
            free_modifier(fmd, true, true);
            fmd->last_frame = 1;
        }
        else
        {
            free_simulation(fmd, false, true);
            fmd->last_frame = 1;
        }

        //try to exec handlers only ONCE
        if (fmd->frac_mesh == NULL) {
            // fire a callback which can then load external data right NOW
            BLI_callback_exec(G.main, &ob->id, BLI_CB_EVT_FRACTURE_REFRESH);
            if (fmd->frac_mesh != NULL) {
                bool tmp = fmd->shards_to_islands;

                fmd->shards_to_islands = false;
                BKE_fracture_create_dm(fmd, true, false);
                fmd->shards_to_islands = tmp;

                BKE_fracture_update_visual_mesh(fmd, ob, true);

                //store names here... gahhhh that is so clumsy
                if (names) {
                    MeshIsland *mi;
                    int i = 0, count = 0;
                    count = BLI_listbase_count(&fmd->meshIslands);

                    (*names) = MEM_callocN(sizeof(char*) * 66 * count, "names");
                    for (mi = fmd->meshIslands.first; mi; mi = mi->next) {
                        //names[i] = MEM_callocN(sizeof(char*) * 66, "names");
                        BLI_snprintf((*names)[i], sizeof(mi->name), "%s", mi->name);
                        i++;
                    }
                    mi_count = i;
                }

                add_new_entries(fmd, dm, ob);
            }
        }
    }

    /* here we just create the fracmesh, in dynamic case we add the first sequence entry as well */
    if (fmd->frac_mesh == NULL) {
         add_new_entries(fmd, dm, ob);
    }
}

bool BKE_fracture_dynamic_lookup_mesh_state(FractureModifierData *fmd, int frame, int do_lookup)
{
    bool changed = false;
    bool forward = false;
    bool backward = false;

    backward = ((fmd->last_frame > frame) && fmd->current_mi_entry && fmd->current_mi_entry->prev);
    forward = ((fmd->last_frame < frame) && (fmd->current_mi_entry) && (fmd->current_mi_entry->next != NULL) &&
               (fmd->current_mi_entry->is_new == false));

    if (backward)
    {
        if (do_lookup)
        {
            while (fmd->current_mi_entry && fmd->current_mi_entry->prev &&
                   frame <= fmd->current_mi_entry->prev->frame)
            {
                printf("Jumping backward because %d is smaller than %d\n", frame, fmd->current_mi_entry->prev->frame);
                changed = true;
                //BKE_free_constraints(fmd);
                BKE_get_prev_entries(fmd);
            }
        }
    }
    else if (forward)
    {
        if (do_lookup)
        {
            while ((fmd->current_mi_entry) && (fmd->current_mi_entry->next != NULL) &&
                   (fmd->current_mi_entry->is_new == false) &&
                   frame > fmd->current_mi_entry->frame)
            {
                printf("Jumping forward because %d is greater than %d\n", frame, fmd->current_mi_entry->frame);
                changed = true;
                //BKE_free_constraints(fmd);
                BKE_get_next_entries(fmd);
            }
        }
    }

    if (do_lookup)
    {
        return changed;
    }
    else
    {
        if (forward || backward)
        {
            fmd->modifier.scene->rigidbody_world->flag |= RBW_FLAG_REFRESH_MODIFIERS;
            fmd->modifier.scene->rigidbody_world->flag |= RBW_FLAG_OBJECT_CHANGED;
        }

        return forward || backward;
    }
}

void BKE_fracture_dynamic_islands_create(FractureModifierData *fmd, Object *ob, DerivedMesh* dm,
                                         DerivedMesh *orig_dm, DerivedMesh *old_cached)
{
    double start = 0.0;
    MDeformVert *ivert = NULL;

    copy_m4_m4(fmd->origmat, ob->obmat);

    /* refracture, convert the fracture shards to new meshislands here *
    * shards = fracture datastructure
    * meshisland = simulation datastructure */

    if (fmd->fix_normals)
    {
        start = PIL_check_seconds_timer();
    }

    ivert = do_islands_from_shards(fmd, ob, orig_dm);

    if (fmd->fix_normals) {
        printf("Fixing normals done, %g\n", PIL_check_seconds_timer() - start);
    }

    fill_vgroup(fmd, fmd->visible_mesh_cached, ivert, ob, old_cached);

    printf("Islands: %d\n", BLI_listbase_count(&fmd->meshIslands));
    /* Grrr, due to stupid design of mine (listbase as value in struct instead of pointer)
     * we have to synchronize the lists here again */

    /* need to ensure(!) old pointers keep valid, else the whole meshisland concept is broken */
    fmd->current_mi_entry->visible_dm = fmd->visible_mesh_cached;
    fmd->current_mi_entry->meshIslands = fmd->meshIslands;
}

ShardSequence* BKE_fracture_dynamic_shard_sequence_add(FractureModifierData* fmd, float frame, DerivedMesh* dm, Object *ob)
{
    ShardSequence *ssq = MEM_mallocN(sizeof(ShardSequence), "shard_sequence_add");

    /*copy last state, to be modified now */
    if (fmd->frac_mesh == NULL) {
        Shard *s = NULL;
        bool temp = fmd->shards_to_islands;
        fmd->frac_mesh = BKE_create_fracture_container();

        /*fill with initial shards*/
        if (fmd->shards_to_islands) {
            do_halving(fmd, ob, dm, dm, true, -1);
        }
        else {
            /* create first shard covering the entire mesh */
            s = BKE_create_fracture_shard(dm->getVertArray(dm),
                                          dm->getPolyArray(dm),
                                          dm->getLoopArray(dm),
                                          dm->getEdgeArray(dm),
                                          dm->numVertData,
                                          dm->numPolyData,
                                          dm->numLoopData,
                                          dm->numEdgeData,
                                          true);

            s = BKE_custom_data_to_shard(s, dm);
            s->flag = SHARD_INTACT;
            s->shard_id = 0;
            BLI_addtail(&fmd->frac_mesh->shard_map, s);
            fmd->frac_mesh->shard_count = 1;
        }

        //build fmd->dm here !
        fmd->shards_to_islands = false;
        BKE_fracture_create_dm(fmd, true, false);
        fmd->shards_to_islands = temp;

        ssq->frac_mesh = fmd->frac_mesh;
    }
    else {
        ssq->frac_mesh = copy_fracmesh(fmd->frac_mesh);
    }

    ssq->is_new = true;
    ssq->frame = frame;
    BLI_addtail(&fmd->shard_sequence, ssq);

    return ssq;
}

MeshIslandSequence* BKE_fracture_dynamic_meshisland_sequence_add(FractureModifierData* fmd,
                                                                 float frame, Object *ob, DerivedMesh *dm)
{
    MeshIslandSequence *msq = MEM_mallocN(sizeof(MeshIslandSequence), "meshisland_sequence_add");
    msq->frame = frame;

    if (BLI_listbase_is_empty(&fmd->meshIslands)) {
        msq->meshIslands.first = NULL;
        msq->meshIslands.last = NULL;
        fmd->visible_mesh_cached = CDDM_copy(fmd->dm);
        do_islands_from_shards(fmd, ob, dm);
        msq->meshIslands = fmd->meshIslands;
        msq->visible_dm = fmd->visible_mesh_cached;
        fmd->refresh = false;
        msq->is_new = false;
    }
    else {
        msq->meshIslands.first = NULL;
        msq->meshIslands.last = NULL;
        msq->visible_dm = NULL;
        msq->is_new = true;
    }

    BLI_addtail(&fmd->meshIsland_sequence, msq);

    return msq;
}

void BKE_fracture_dynamic_new_entries_add(FractureModifierData* fmd, DerivedMesh *dm, Object* ob)
{
    int frame = (int)BKE_scene_frame_get(fmd->modifier.scene);
    int end = 250; //TODO might be problematic with longer sequences, take proper end value ?

    if (fmd->modifier.scene->rigidbody_world)
    {
        end = fmd->modifier.scene->rigidbody_world->pointcache->endframe;
    }

    if (fmd->current_shard_entry)
    {
        fmd->current_shard_entry->is_new = false;
        fmd->current_shard_entry->frame = frame;
    }
    fmd->current_shard_entry = shard_sequence_add(fmd, end, dm, ob);
    fmd->frac_mesh = fmd->current_shard_entry->frac_mesh;

    if (fmd->current_mi_entry) {
        fmd->current_mi_entry->frame = frame;
    }

    fmd->current_mi_entry = meshisland_sequence_add(fmd, end, ob, dm);
    fmd->meshIslands = fmd->current_mi_entry->meshIslands;
}

void BKE_fracture_dynamic_sequences_free(FractureModifierData *fmd)
{
    /* in dynamic mode we have to get rid of the entire Meshisland sequence */
    MeshIslandSequence *msq;
    ShardSequence *ssq;

    while (fmd->meshIsland_sequence.first) {
        msq = fmd->meshIsland_sequence.first;
        BLI_remlink(&fmd->meshIsland_sequence, msq);
        free_meshislands(fmd, &msq->meshIslands, do_free_rigidbody);
        MEM_freeN(msq);
        msq = NULL;
    }

    fmd->meshIsland_sequence.first = NULL;
    fmd->meshIsland_sequence.last = NULL;

    fmd->meshIslands.first = NULL;
    fmd->meshIslands.last = NULL;

    fmd->current_mi_entry = NULL;

    while (fmd->shard_sequence.first)
    {
        ssq = fmd->shard_sequence.first;
        BLI_remlink(&fmd->shard_sequence, ssq);
        BKE_fracmesh_free(ssq->frac_mesh, true);
        MEM_freeN(ssq->frac_mesh);
        MEM_freeN(ssq);
    }

    fmd->shard_sequence.first = NULL;
    fmd->shard_sequence.last = NULL;
    fmd->current_shard_entry = NULL;
    fmd->frac_mesh = NULL;
}

void BKE_fracture_dynamic_free(FractureModifierData *fmd, bool do_free_seq, bool do_free_rigidbody)
{
    if ((!fmd->refresh && !fmd->refresh_constraints)) {
        /* free entire modifier */
        BKE_fracture_modifier_free(fmd, do_free_seq, do_free_rigidbody);
    }
    else if (fmd->refresh_constraints && !fmd->is_dynamic_external) {
        /* refresh constraints only */
        BKE_free_constraints(fmd);
    }
}
#endif
