#ifndef SPECTRAL_MIX_H
#define SPECTRAL_MIX_H

#define SPD_N 36

void rgb_to_spectral(float r, float g, float b, float *out_spd);
void spectral_mix_geometric(const float *A, const float *B,
                            float wA, float wB, float *out);
void spectral_to_rgb(const float *spd,
                     float *out_r, float *out_g, float *out_b,
                     int apply_srgb_gamma);

void pigment_mix_rgb(float r1,float g1,float b1,
                     float r2,float g2,float b2,
                     float t,
                     float *or,float *og,float *ob,
                     int apply_srgb_gamma);

#endif