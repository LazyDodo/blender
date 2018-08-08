#include "MEM_guardedalloc.h"
#include "BKE_fracture.h"
#include "BKE_main.h"
#include "BKE_mesh.h"
#include "BKE_pointcache.h"

#include "BLI_listbase.h"
#include "BLI_string.h"

#include "DNA_fracture_types.h"
#include "DNA_object_types.h"
#include "DNA_mesh_types.h"
#include "DNA_meshdata_types.h"
#include "DNA_modifier_types.h"
#include "DNA_scene_types.h"
#include "DNA_rigidbody_types.h"

#include "DEG_depsgraph_query.h"

#include "limits.h"
#include "PIL_time.h"

static void fracture_dynamic_new_entries_add(FractureModifierData* fmd, Mesh *dm, Object* ob, Scene *scene);
static void fracture_dynamic_initialize(FractureModifierData *fmd, Object *ob, Mesh *dm, char (**names)[66], Scene *scene);
static void get_next_entries(FractureModifierData *fmd);
static void get_prev_entries(FractureModifierData *fmd);

Mesh *BKE_fracture_dynamic_apply(FractureModifierData *fmd, Object *ob, Mesh *derivedData, Scene* scene)
{
    char (*names)[66] = NULL;
    bool valid_fractured_mesh = false;

    if (fmd->refresh_constraints)
    {
        BKE_fracture_constraints_free(fmd, scene);
    }

    /* group_dm, clean_dm not necessary here as we dont support non-mesh objects and subobject_groups here */
    // this should be done from an operator, really...
    if (fmd->refresh)
    {
        fracture_dynamic_initialize(fmd, ob, derivedData, &names, scene);

        fmd->shared->current_mi_entry->is_new = false;
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
    valid_fractured_mesh = !fmd->explo_shared || (fmd->explo_shared && fmd->shared->dm && fmd->shared->frac_mesh);

   /* if ((!valid_fractured_mesh) || (fmd->visible_mesh == NULL && fmd->visible_mesh_cached == NULL)) {
        do_clear(fmd);
    }*/

    return BKE_fracture_result_mesh(fmd, derivedData, ob, valid_fractured_mesh, scene);
}

static void fracture_dynamic_update(FractureModifierData *fmd, Object* ob, Scene* scene, Mesh* dm, Depsgraph* depsgraph, Main* bmain)
{
    int frame = (int)BKE_scene_frame_get(scene);
    /*HERE we must know which shard(s) to fracture... hmm shards... we should "merge" states which happen in the same frame automatically !*/
    /* TODO_2 dynamic, this is called from rigidbody system only !!! so move out of the md loop as well, to BKE */
    if (!(BKE_fracture_dynamic_lookup_mesh_state(fmd, frame, false, scene)))
    {
        /*simulation mode*/
        /* bullet callbacks may happen multiple times per frame, in next frame we can evaluate them all,
         * so we need some array of shardIDs or shards to fracture each *
         * we need to loop over those shard IDs here, but lookup of shard ids might be slow, but fracturing of many shards is slower...
         * should not have a visible effect in general */

        int count = 0;
        int totpoint = 0;
        FractureID *fid = fmd->shared->fracture_ids.first;

        while(fid) {
            FracPointCloud points;
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
                    BKE_fracture_constraints_free(fmd, scene);
                }

                printf("ADD NEW 2: %s \n", ob->id.name);
                fmd->update_dynamic = false;
                fracture_dynamic_new_entries_add(fmd, dm, ob, scene);
            }

            while(fmd->shared->fracture_ids.first){
                fid = (FractureID*)fmd->shared->fracture_ids.first;
                BKE_fracture_do(fmd, fid->shardID, ob, dm, depsgraph, bmain);
                BLI_remlink(&fmd->shared->fracture_ids, fid);
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

static void fracture_dynamic_initialize(FractureModifierData *fmd, Object *ob, Mesh *dm, char (**names)[66], Scene* scene)
{
    /*TODO_1 refresh, move to BKE and just call from operator for prefractured case*/

    printf("ADD NEW 1: %s \n", ob->id.name);
    if ((fmd->last_frame == INT_MAX)) //TODO event if we restart at cache startframe
    {
        if (fmd->reset_shards)
        {
            BKE_fracture_simulation_free(fmd, true, true, scene);
            BKE_fracture_free(fmd, true, true, scene);
            fmd->last_frame = 1;
        }
        else
        {
            BKE_fracture_simulation_free(fmd, false, true, scene);
            fmd->last_frame = 1;
        }

        //try to exec handlers only ONCE
        if (fmd->shared->frac_mesh == NULL) {
            // fire a callback which can then load external data right NOW
            //BLI_callback_exec(G.main, &ob->id, BLI_CB_EVT_FRACTURE_REFRESH);
            if (fmd->shared->frac_mesh != NULL) {
                bool tmp = fmd->shards_to_islands;

                fmd->shards_to_islands = false;
                BKE_fracture_assemble_mesh_from_shards(fmd, true, false);
                fmd->shards_to_islands = tmp;

                BKE_fracture_update_visual_mesh(fmd, ob, true);

                //store names here... gahhhh that is so clumsy
                if (names) {
                    MeshIsland *mi;
                    int i = 0, count = 0;
                    count = BLI_listbase_count(&fmd->shared->meshIslands);

                    (*names) = MEM_callocN(sizeof(char*) * 66 * count, "names");
                    for (mi = fmd->shared->meshIslands.first; mi; mi = mi->next) {
                        //names[i] = MEM_callocN(sizeof(char*) * 66, "names");
                        BLI_snprintf((*names)[i], sizeof(mi->name), "%s", mi->name);
                        i++;
                    }
                    //mi_count = i;
                }

                fracture_dynamic_new_entries_add(fmd, dm, ob, scene);
            }
        }
    }

    /* here we just create the fracmesh, in dynamic case we add the first sequence entry as well */
    if (fmd->shared->frac_mesh == NULL) {
         fracture_dynamic_new_entries_add(fmd, dm, ob, scene);
    }
}

bool BKE_fracture_dynamic_lookup_mesh_state(FractureModifierData *fmd, int frame, int do_lookup, Scene* scene)
{
    bool changed = false;
    bool forward = false;
    bool backward = false;

    backward = ((fmd->last_frame > frame) && fmd->shared->current_mi_entry && fmd->shared->current_mi_entry->prev);
    forward = ((fmd->last_frame < frame) && (fmd->shared->current_mi_entry) && (fmd->shared->current_mi_entry->next != NULL) &&
               (fmd->shared->current_mi_entry->is_new == false));

    if (backward)
    {
        if (do_lookup)
        {
            while (fmd->shared->current_mi_entry && fmd->shared->current_mi_entry->prev &&
                   frame <= fmd->shared->current_mi_entry->prev->frame)
            {
                printf("Jumping backward because %d is smaller than %d\n", frame, fmd->shared->current_mi_entry->prev->frame);
                changed = true;
                get_prev_entries(fmd);
            }
        }
    }
    else if (forward)
    {
        if (do_lookup)
        {
            while ((fmd->shared->current_mi_entry) && (fmd->shared->current_mi_entry->next != NULL) &&
                   (fmd->shared->current_mi_entry->is_new == false) &&
                   frame > fmd->shared->current_mi_entry->frame)
            {
                printf("Jumping forward because %d is greater than %d\n", frame, fmd->shared->current_mi_entry->frame);
                changed = true;
                get_next_entries(fmd);
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
            scene->rigidbody_world->flag |= RBW_FLAG_REFRESH_MODIFIERS;
            scene->rigidbody_world->flag |= RBW_FLAG_OBJECT_CHANGED;
        }

        return forward || backward;
    }
}

static void fracture_dynamic_islands_create(FractureModifierData *fmd, Object *ob,
                                         Mesh *orig_dm, Mesh *old_cached, Scene* scene)
{
    double start = 0.0;
    MDeformVert *ivert = NULL;

  //  copy_m4_m4(fmd->origmat, ob->obmat);

    /* refracture, convert the fracture shards to new meshislands here *
    * shards = fracture datastructure
    * meshisland = simulation datastructure */

    if (fmd->fix_normals)
    {
        start = PIL_check_seconds_timer();
    }

    ivert = BKE_fracture_shards_to_islands(fmd, ob, orig_dm, scene);

    if (fmd->fix_normals) {
        printf("Fixing normals done, %g\n", PIL_check_seconds_timer() - start);
    }

    BKE_fracture_fill_vgroup(fmd, fmd->shared->visible_mesh_cached, ivert, ob, old_cached);

    printf("Islands: %d\n", BLI_listbase_count(&fmd->shared->meshIslands));
    /* Grrr, due to stupid design of mine (listbase as value in struct instead of pointer)
     * we have to synchronize the lists here again */

    /* need to ensure(!) old pointers keep valid, else the whole meshisland concept is broken */
    fmd->shared->current_mi_entry->visible_dm = fmd->shared->visible_mesh_cached;
    fmd->shared->current_mi_entry->meshIslands = fmd->shared->meshIslands;
}

static ShardSequence* fracture_dynamic_shard_sequence_add(FractureModifierData* fmd, float frame, Mesh* dm, Object *ob, Scene* scene)
{
    ShardSequence *ssq = MEM_mallocN(sizeof(ShardSequence), "shard_sequence_add");

    /*copy last state, to be modified now */
    if (fmd->shared->frac_mesh == NULL) {
        Shard *s = NULL;
        bool temp = fmd->shards_to_islands;
        fmd->shared->frac_mesh = BKE_fracture_fracmesh_create();

        /*fill with initial shards*/
        if (fmd->shards_to_islands) {
            BKE_fracture_do_halving(fmd, ob, dm, dm, true, -1, scene);
        }
        else {
            /* create first shard covering the entire mesh */
            s = BKE_fracture_shard_create(dm->mvert,
                                          dm->mpoly,
                                          dm->mloop,
                                          dm->medge,
                                          dm->totvert,
                                          dm->totpoly,
                                          dm->totloop,
                                          dm->totedge,
                                          true);

            BKE_fracture_custom_data_mesh_to_shard(s, dm);
            s->flag = SHARD_INTACT;
            s->shard_id = 0;
            BLI_addtail(&fmd->shared->frac_mesh->shard_map, s);
            fmd->shared->frac_mesh->shard_count = 1;
        }

        //build fmd->dm here !
        fmd->shards_to_islands = false;
        BKE_fracture_assemble_mesh_from_shards(fmd, true, false);
        fmd->shards_to_islands = temp;

        ssq->frac_mesh = fmd->shared->frac_mesh;
    }
    else {
        ssq->frac_mesh = BKE_fracture_fracmesh_copy(fmd->shared->frac_mesh);
    }

    ssq->is_new = true;
    ssq->frame = frame;
    BLI_addtail(&fmd->shared->shard_sequence, ssq);

    return ssq;
}

static MeshIslandSequence* fracture_dynamic_meshisland_sequence_add(FractureModifierData* fmd,
                                                                 float frame, Object *ob, Mesh *dm, Scene* scene)
{
    MeshIslandSequence *msq = MEM_mallocN(sizeof(MeshIslandSequence), "meshisland_sequence_add");
    msq->frame = (int)frame;

    if (BLI_listbase_is_empty(&fmd->shared->meshIslands)) {
        msq->meshIslands.first = NULL;
        msq->meshIslands.last = NULL;
        fmd->shared->visible_mesh_cached = BKE_fracture_mesh_copy(fmd->shared->dm, ob);
        BKE_fracture_shards_to_islands(fmd, ob, dm, scene);
        msq->meshIslands = fmd->shared->meshIslands;
        msq->visible_dm = fmd->shared->visible_mesh_cached;
        fmd->refresh = false;
        msq->is_new = false;
    }
    else {
        msq->meshIslands.first = NULL;
        msq->meshIslands.last = NULL;
        msq->visible_dm = NULL;
        msq->is_new = true;
    }

    BLI_addtail(&fmd->shared->meshIsland_sequence, msq);

    return msq;
}

static void fracture_dynamic_new_entries_add(FractureModifierData* fmd, Mesh *dm, Object* ob, Scene *scene)
{
    int frame = (int)BKE_scene_frame_get(scene);
    int end = 250; //TODO might be problematic with longer sequences, take proper end value ?

    if (scene->rigidbody_world)
    {
        end = scene->rigidbody_world->shared->pointcache->endframe;
    }

    if (fmd->shared->current_shard_entry)
    {
        fmd->shared->current_shard_entry->is_new = false;
        fmd->shared->current_shard_entry->frame = frame;
    }
    fmd->shared->current_shard_entry = fracture_dynamic_shard_sequence_add(fmd, end, dm, ob, scene);
    fmd->shared->frac_mesh = fmd->shared->current_shard_entry->frac_mesh;

    if (fmd->shared->current_mi_entry) {
        fmd->shared->current_mi_entry->frame = frame;
    }

    fmd->shared->current_mi_entry = fracture_dynamic_meshisland_sequence_add(fmd, end, ob, dm, scene);
    fmd->shared->meshIslands = fmd->shared->current_mi_entry->meshIslands;
}

static void fracture_dynamic_sequences_free(FractureModifierData *fmd, Scene* scene)
{
    /* in dynamic mode we have to get rid of the entire Meshisland sequence */
    MeshIslandSequence *msq;
    ShardSequence *ssq;

    while (fmd->shared->meshIsland_sequence.first) {
        msq = fmd->shared->meshIsland_sequence.first;
        BLI_remlink(&fmd->shared->meshIsland_sequence, msq);
        BKE_fracture_meshislands_free(fmd, &msq->meshIslands, true, scene);
        MEM_freeN(msq);
        msq = NULL;
    }

    fmd->shared->meshIsland_sequence.first = NULL;
    fmd->shared->meshIsland_sequence.last = NULL;

    fmd->shared->meshIslands.first = NULL;
    fmd->shared->meshIslands.last = NULL;

    fmd->shared->current_mi_entry = NULL;

    while (fmd->shared->shard_sequence.first)
    {
        ssq = fmd->shared->shard_sequence.first;
        BLI_remlink(&fmd->shared->shard_sequence, ssq);
        BKE_fracture_fracmesh_free(ssq->frac_mesh, true);
        MEM_freeN(ssq->frac_mesh);
        MEM_freeN(ssq);
    }

    fmd->shared->shard_sequence.first = NULL;
    fmd->shared->shard_sequence.last = NULL;
    fmd->shared->current_shard_entry = NULL;
    fmd->shared->frac_mesh = NULL;
}

void BKE_fracture_dynamic_free(FractureModifierData *fmd, bool do_free_seq, bool do_free_rigidbody, Scene* scene)
{
    if ((!fmd->refresh && !fmd->refresh_constraints)) {
        /* free entire modifier */
        BKE_fracture_modifier_free(fmd, do_free_seq, do_free_rigidbody, scene);
    }
    else if (fmd->refresh_constraints && !fmd->is_dynamic_external) {
        /* refresh constraints only */
        BKE_fracture_constraints_free(fmd, scene);
    }
}

static void get_next_entries(FractureModifierData *fmd)
{
    /*meshislands and shards SHOULD be synchronized !!!!*/
    if (fmd->shared->current_mi_entry && fmd->shared->current_mi_entry->next)
    { // && fmd->current_mi_entry->next->is_new == false) {
        fmd->shared->current_mi_entry = fmd->shared->current_mi_entry->next;
        fmd->shared->meshIslands = fmd->shared->current_mi_entry->meshIslands;
        fmd->shared->visible_mesh_cached = fmd->shared->current_mi_entry->visible_dm;
    }

    if (fmd->shared->current_shard_entry && fmd->shared->current_shard_entry->next)
    {
        fmd->shared->current_shard_entry = fmd->shared->current_shard_entry->next;
        fmd->shared->frac_mesh = fmd->shared->current_shard_entry->frac_mesh;
    }
}

static void get_prev_entries(FractureModifierData *fmd)
{
    /*meshislands and shards SHOULD be synchronized !!!!*/
    if (fmd->shared->current_mi_entry && fmd->shared->current_mi_entry->prev)
    {
        fmd->shared->current_mi_entry = fmd->shared->current_mi_entry->prev;
        fmd->shared->meshIslands = fmd->shared->current_mi_entry->meshIslands;
        fmd->shared->visible_mesh_cached = fmd->shared->current_mi_entry->visible_dm;
    }

    if (fmd->shared->current_shard_entry && fmd->shared->current_shard_entry->prev)
    {
        fmd->shared->current_shard_entry = fmd->shared->current_shard_entry->prev;
        fmd->shared->frac_mesh = fmd->shared->current_shard_entry->frac_mesh;
    }
}
