#include "BKE_fracture.h"
#include "BKE_pointcache.h"
#include "BKE_rigidbody.h"

#include "DEG_depsgraph_query.h"

#include "DNA_modifier_types.h"
#include "DNA_object_types.h"
#include "DNA_mesh_types.h"
#include "DNA_meshdata_types.h"
#include "DNA_fracture_types.h"
#include "DNA_rigidbody_types.h"
#include "DNA_scene_types.h"


Mesh *BKE_fracture_prefractured_apply(FractureModifierData *fmd, Object *ob, Mesh *derivedData, Depsgraph* depsgraph)
{
    bool do_refresh = (fmd->auto_execute) || (fmd->dm_group && fmd->use_constraint_group && fmd->refresh_constraints);
    Scene *scene = fmd->scene;

    Mesh *final_dm = derivedData;
  //  Mesh *group_dm = BKE_fracture_group_dm(fmd, derivedData, ob, do_refresh || fmd->refresh);

    if (do_refresh) {
        fmd->refresh = true;
    }

    if (fmd->refresh)
    {
        BKE_fracture_initialize(fmd, ob, derivedData, depsgraph);
    }

    /* TODO_5, get rid of fmd->dm and perhaps of fmd->visible_mesh (BMESH!) too, the latter should be runtime data for creating islands ONLY */
    /* we should ideally only have one cached derivedmesh */
    if (fmd->dm && fmd->frac_mesh && (fmd->dm->totpoly > 0)) {
        final_dm = BKE_fracture_prefractured_do(fmd, ob, fmd->dm, derivedData, NULL, 0, scene, depsgraph);
    }
    else {
        final_dm = BKE_fracture_prefractured_do(fmd, ob, derivedData, derivedData, NULL, 0, scene, depsgraph);
    }

    return final_dm;
}
