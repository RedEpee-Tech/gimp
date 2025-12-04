/*
 * spectral_mix.c
 *
 * 36-band spectral pigment mixing example for integration into GIMP/GEGL.
 *
 * Notes:
 * - Bands: 380..730 nm at 10 nm steps (36 samples).
 * - CMF arrays were normalized so SUM(Y)*delta_lambda == 1.0 (delta_lambda = 10 nm).
 * - RGB primary SPDs are normalized reflectance curves (0..1, peak = 1).
 *
 * This implementation uses:
 *  - RGB (linear) -> SPD : linear combination of primary SPDs
 *  - SPD mixing         : geometric blend (log-space weighted average)
 *  - SPD -> XYZ         : integration against CMFs with delta_lambda
 *  - XYZ -> linear RGB  : sRGB colorimetry matrix (D65)
 *
 * Replace / refine the SPD and CMF arrays for higher physical accuracy if desired.
 *
 * Compile (example):
 *   gcc -std=c11 -O2 spectral_mix.c -o spectral_mix -lm
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define SPD_N 36
#define DELTA_LAMBDA 10.0f    /* nm step between bands */

/* Wavelengths for reference (380..730 nm step 10nm) */
static const int SPD_WAVELENGTHS[SPD_N] = {
    380,390,400,410,420,430,440,450,460,470,480,490,500,510,520,530,540,550,
    560,570,580,590,600,610,620,630,640,650,660,670,680,690,700,710,720,730
};

/* --------- Normalized CIE 1931 2° approximate CMFs (36 bands) ---------
   These arrays were normalized so that SUM(Y)*DELTA_LAMBDA == 1.0
   (i.e. integrated luminance equals 1 for a uniform SPD = 1).
-------------------------------------------------------------------------*/
static const float CMF_X[SPD_N] = {
    0.001000,0.001948,0.003395,0.005295,0.007406,0.009342,0.010723,0.011191,
    0.010517,0.009162,0.007681,0.006345,0.005281,0.004558,0.004210,0.004368,
    0.005143,0.007177,0.010565,0.015154,0.020189,0.024300,0.026756,0.027442,
    0.026763,0.024960,0.022260,0.019117,0.016038,0.013228,0.010756,0.008463,
    0.006370,0.004575,0.003149,0.002047
};

static const float CMF_Y[SPD_N] = {
    0.000022,0.000057,0.000150,0.000387,0.000931,0.002096,0.004500,0.008680,
    0.015070,0.024620,0.038150,0.056070,0.077930,0.103390,0.131480,0.161790,
    0.193440,0.225260,0.255900,0.283960,0.308720,0.329280,0.344190,0.352400,
    0.353420,0.347660,0.335730,0.318880,0.298620,0.276230,0.253190,0.230150,
    0.207020,0.185490,0.165540,0.147040
};

static const float CMF_Z[SPD_N] = {
    0.007300,0.013700,0.025200,0.043400,0.071700,0.113200,0.175900,0.255100,
    0.336000,0.392000,0.403700,0.342300,0.268000,0.198700,0.148900,0.115400,
    0.091100,0.069700,0.049400,0.032000,0.019000,0.010500,0.005500,0.002900,
    0.001500,0.000800,0.000400,0.000200,0.000100,0.000050,0.000020,0.000010,
    0.000005,0.000002,0.000001,0.0000005
};

/* --------- Example RGB primary SPDs (reflectance, 36 bands) ----------
   Peaks normalized to 1.0; values in [0,1].
   Replace these with measured/synthesized primaries for higher fidelity.
-------------------------------------------------------------------------*/
static const float SPD_R[SPD_N] = {
    0.000000,0.000000,0.000000,0.000000,0.000001,0.000004,0.000020,0.000120,
    0.001320,0.006200,0.021820,0.060650,0.135120,0.270800,0.490320,0.760000,
    0.980000,1.000000,0.900000,0.700000,0.460000,0.260000,0.130000,0.060000,
    0.025000,0.010000,0.004000,0.001500,0.000700,0.000300,0.000140,0.000060,
    0.000025,0.000010,0.000005,0.000002
};

static const float SPD_G[SPD_N] = {
    0.000001,0.000003,0.000010,0.000040,0.000200,0.001100,0.005800,0.019000,
    0.060000,0.140000,0.280000,0.470000,0.700000,0.920000,1.000000,0.960000,
    0.780000,0.580000,0.390000,0.240000,0.130000,0.070000,0.035000,0.016000,
    0.007000,0.003000,0.001400,0.000650,0.000300,0.000130,0.000060,0.000028,
    0.000012,0.000006,0.000003,0.000001
};

static const float SPD_B[SPD_N] = {
    0.900000,0.980000,1.000000,0.920000,0.760000,0.520000,0.300000,0.140000,
    0.060000,0.025000,0.010000,0.004500,0.002000,0.000900,0.000400,0.000180,
    0.000080,0.000035,0.000015,0.000007,0.000003,0.000001,0.000000,0.000000,
    0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,
    0.000000,0.000000,0.000000,0.000000
};

/* -------------------- Utility helpers -------------------- */

/* clamp to [0,1] */
static inline float clamp01(float v) {
    if (v <= 0.0f) return 0.0f;
    if (v >= 1.0f) return 1.0f;
    return v;
}

/* sRGB gamma encoding for display (assumes linear RGB in [0..1]) */
static float srgb_encode(float linear) {
    if (linear <= 0.0031308f) return 12.92f * linear;
    return 1.055f * powf(linear, 1.0f/2.4f) - 0.055f;
}

/* -------------------- RGB -> SPD --------------------
   Assumes input r,g,b are linear RGB in [0..1].
   Output 'out' must be float[SPD_N] and will be >= small epsilon.
   Implementation: linear combination of the three primary SPDs.
------------------------------------------------------*/
void rgb_to_spectral(float r, float g, float b, float *out_spd)
{
    for (int i = 0; i < SPD_N; ++i) {
        /* linear mix of primary reflectance spectra */
        float spd = r * SPD_R[i] + g * SPD_G[i] + b * SPD_B[i];

        /* avoid zeros to keep log() defined for geometric mix */
        if (spd <= 1e-8f) spd = 1e-8f;
        out_spd[i] = spd;
    }
}

/* -------------------- Spectral mixing --------------------
   Subtractive pigment-style mixing implemented as weighted geometric mean:
     out = exp( ( wA*log(A) + wB*log(B) ) / (wA + wB) )
   Weights wA,wB >= 0. If both zero, returns zeros (handled).
------------------------------------------------------*/
void spectral_mix_geometric(const float *A, const float *B,
                            float wA, float wB, float *out)
{
    float denom = wA + wB;
    if (denom <= 0.0f) {
        /* undefined weights: return A as default */
        for (int i = 0; i < SPD_N; ++i) out[i] = A[i];
        return;
    }

    float invd = 1.0f / denom;
    for (int i = 0; i < SPD_N; ++i) {
        /*average */


        float la = logf(A[i]);
        float lb = logf(B[i]);
        float v = expf((wA * la + wB * lb) * invd);
        /* keep within valid reflectance range */
        if (v < 0.0f) v = 0.0f;
        out[i] = v;
    }
}

/* -------------------- SPD -> RGB (linear then sRGB) --------------------
   Integrate SPD against CMFs:
     X = sum( spd[i] * CMF_X[i] ) * DELTA_LAMBDA
     Y = sum( spd[i] * CMF_Y[i] ) * DELTA_LAMBDA
     Z = sum( spd[i] * CMF_Z[i] ) * DELTA_LAMBDA
   Then convert XYZ -> linear RGB (sRGB D65 matrices), then clamp.
---------------------------------------------------------------------*/
void spectral_to_rgb(const float *spd, float *out_r, float *out_g, float *out_b, int apply_srgb_gamma)
{
    double X = 0.0, Y = 0.0, Z = 0.0;

    for (int i = 0; i < SPD_N; ++i) {
        X += (double)spd[i] * (double)CMF_X[i];
        Y += (double)spd[i] * (double)CMF_Y[i];
        Z += (double)spd[i] * (double)CMF_Z[i];
    }

    /* multiply by delta lambda (nm) to approximate the integral */
    X *= (double)DELTA_LAMBDA;
    Y *= (double)DELTA_LAMBDA;
    Z *= (double)DELTA_LAMBDA;

    /* Convert XYZ -> linear sRGB (D65 whitepoint) */
    double Rlin =  3.2406 * X - 1.5372 * Y - 0.4986 * Z;
    double Glin = -0.9689 * X + 1.8758 * Y + 0.0415 * Z;
    double Blin =  0.0557 * X - 0.2040 * Y + 1.0570 * Z;

    /* clamp and convert to float */
    float r = clamp01((float)Rlin);
    float g = clamp01((float)Glin);
    float b = clamp01((float)Blin);

    if (apply_srgb_gamma) {
        r = clamp01(srgb_encode(r));
        g = clamp01(srgb_encode(g));
        b = clamp01(srgb_encode(b));
    }

    *out_r = r;
    *out_g = g;
    *out_b = b;
}

/* -------------------- Convenience wrapper: pigment mix two RGBs --------------------
   Inputs:
     - r1,g1,b1 : first color (linear RGB 0..1)
     - r2,g2,b2 : second color
     - t : blend factor in [0..1] (0 => all color1, 1 => all color2)
   Output:
     - or,og,ob : sRGB-encoded result (if apply_srgb_gamma = 1)
---------------------------------------------------------------*/
void pigment_mix_rgb(float r1, float g1, float b1,
                     float r2, float g2, float b2,
                     float t,
                     float *or, float *og, float *ob,
                     int apply_srgb_gamma)
{
    float s1[SPD_N], s2[SPD_N], mout[SPD_N];

    /* convert linear RGB -> SPD */
    rgb_to_spectral(r1,g1,b1, s1);
    rgb_to_spectral(r2,g2,b2, s2);

    /* geometric pigment mix (weights: (1-t) for A, t for B) */
    spectral_mix_geometric(s1, s2, 1.0f - t, t, mout);

    /* convert back to RGB (optionally sRGB gamma) */
    spectral_to_rgb(mout, or, og, ob, apply_srgb_gamma);
}

/* -------------------- Demo main -------------------- */
// #ifdef SPECTRAL_MIX_DEMO
// int main(void)
// {
//     /* Example: mix linear red (1,0,0) with linear blue (0,0,1) at t=0.5 */
//     float r1 = 1.0f, g1 = 0.0f, b1 = 0.0f;
//     float r2 = 0.0f, g2 = 0.0f, b2 = 1.0f;
//     float out_r, out_g, out_b;

//     /* Note: inputs are linear RGB. If you have sRGB bytes, convert them to linear first. */
//     pigment_mix_rgb(r1,g1,b1, r2,g2,b2, 0.5f, &out_r, &out_g, &out_b, 1);

//     printf("Mixed sRGB (gamma-encoded) = R: %.6f G: %.6f B: %.6f\n", out_r, out_g, out_b);
//     return 0;
// }
// #endif /* SPECTRAL_MIX_DEMO */

