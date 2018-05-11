
struct ViewLayer;
struct Scene;

void *EEVEE_lightcache_job_data_alloc(struct Main *bmain, struct ViewLayer *viewlayer, struct Scene *scene);
void EEVEE_lightcache_job_data_free(void *custom_data);
void EEVEE_lightcache_bake_job(void *custom_data, short *stop, short *do_update, float *progress);