#include "sos_filter.h"
#include <stdlib.h>

void sos_filter_init(sos_filter_t *filter, uint8_t num_sections)
{
    filter->num_sections = num_sections;
    filter->sections = (biquad_filter_t *)malloc(num_sections * sizeof(biquad_filter_t));
    if (filter->sections == NULL) {
        filter->num_sections = 0;
        return;
    }

    for (uint8_t i = 0; i < num_sections; i++) {
        biquad_filter_set_empty(&filter->sections[i]);
    }
}

void sos_filter_destroy(sos_filter_t *filter)
{
    free(filter->sections);
    filter->sections = NULL;
    filter->num_sections = 0;
}

void sos_filter_set_section(sos_filter_t *filter, uint8_t section_index,
                             const float num_z[3], const float den_z[3])
{
    if (section_index >= filter->num_sections) {
        return;
    }
    biquad_filter_init(&filter->sections[section_index], num_z, den_z);
}

float sos_filter_update(sos_filter_t *filter, float input)
{
    float x = input;
    for (uint8_t i = 0; i < filter->num_sections; i++) {
        x = biquad_filter_update(&filter->sections[i], x);
    }
    return x;
}

float sos_filter_get_output(const sos_filter_t *filter)
{
    if (filter->num_sections == 0) {
        return 0.0f;
    }
    return biquad_filter_get_output(&filter->sections[filter->num_sections - 1]);
}

float sos_filter_get_input(const sos_filter_t *filter)
{
    if (filter->num_sections == 0) {
        return 0.0f;
    }
    return biquad_filter_get_input(&filter->sections[0]);
}

void sos_filter_reset(sos_filter_t *filter, float equilibrium)
{
    float x = equilibrium;
    for (uint8_t i = 0; i < filter->num_sections; i++) {
        biquad_filter_reset(&filter->sections[i], x);
        x = biquad_filter_get_output(&filter->sections[i]);
    }
}
