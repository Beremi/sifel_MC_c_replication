#include "matmodel.h"
#include "vector.h"
#include "matrix.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace
{

    struct spectral_data
    {
        vector eig;
        vector proj[3];
        matrix hess[3];

        spectral_data() : eig(3)
        {
            for (long i = 0; i < 3; i++)
            {
                reallocv(MATMODEL_NCOMP_STRAIN, proj[i]);
                reallocm(MATMODEL_NCOMP_STRAIN, MATMODEL_NCOMP_STRAIN, hess[i]);
            }
        }
    };

    struct matmodel_cache
    {
        vector epsp;
        int return_type;

        vector proj[3];
        matrix hess[3];
        vector sig_princ;
        vector stress;

        matmodel_cache() : epsp(MATMODEL_NCOMP_EQOTHER), return_type(MATMODEL_RETURN_ELASTIC), sig_princ(3), stress(MATMODEL_NCOMP_STRESS)
        {
            for (long i = 0; i < 3; i++)
            {
                reallocv(MATMODEL_NCOMP_STRAIN, proj[i]);
                reallocm(MATMODEL_NCOMP_STRAIN, MATMODEL_NCOMP_STRAIN, hess[i]);
            }
        }
    };

    void init_params(matmodel_params *par,
                     double young,
                     double poisson,
                     double cohesion,
                     double phi)
    {
        if (par == NULL)
            return;

        par->young = young;
        par->poisson = poisson;
        par->c = cohesion;
        par->phi = phi;
        par->shear = young / (2.0 * (1.0 + poisson));
        par->bulk = young / (3.0 * (1.0 - 2.0 * poisson));
        par->lame = par->bulk - 2.0 * par->shear / 3.0;
        par->sin_phi = sin(phi);
        par->cos_phi = cos(phi);
        par->c_bar = 2.0 * cohesion * par->cos_phi;
    }

    void null_cache(matmodel_cache &cache)
    {
        nullv(cache.epsp);
        cache.return_type = MATMODEL_RETURN_ELASTIC;
        nullv(cache.sig_princ);
        nullv(cache.stress);
        for (long i = 0; i < 3; i++)
        {
            nullv(cache.proj[i]);
            nullm(cache.hess[i]);
        }
    }

    void copy_mode(vector &dst_proj,
                   matrix &dst_hess,
                   double &dst_eig,
                   const vector &src_proj,
                   const matrix &src_hess,
                   double src_eig)
    {
        copyv(src_proj, dst_proj);
        copym(src_hess, dst_hess);
        dst_eig = src_eig;
    }

} // namespace

matmodel::matmodel()
    : cached_strain(MATMODEL_NCOMP_STRAIN),
      cached_eqother(MATMODEL_NCOMP_EQOTHER),
      cached_other(MATMODEL_NCOMP_OTHER)
{
    init_params(&par, 1.0, 0.25, 0.0, 0.52359877559829887308);
    nullv(cached_strain);
    nullv(cached_eqother);
    nullv(cached_other);
    cached_response_valid = 0;
}

long matmodel::read(FILE *in)
{
    matmodel_params tmp;
    double young, poisson, cohesion, phi;

    if (in == NULL)
        return 1;

    if (fscanf(in, "%lf %lf %lf %lf", &young, &poisson, &cohesion, &phi) != 4)
        return 1;

    init_params(&tmp, young, poisson, cohesion, phi);
    if (tmp.young <= 0.0)
        return 1;
    if ((tmp.poisson <= -1.0) || (tmp.poisson >= 0.5))
        return 1;
    if (tmp.c < 0.0)
        return 1;
    if ((tmp.phi <= 1.0e-14) || (tmp.phi >= 0.5 * M_PI - 1.0e-14))
        return 1;

    par = tmp;
    cached_response_valid = 0;
    return 0;
}

void matmodel::print(FILE *out)
{
    if (out == NULL)
        return;

    fprintf(out, "MC perfect plastic associative 2D model\n");
    fprintf(out, "  E      = %.15g\n", par.young);
    fprintf(out, "  nu     = %.15g\n", par.poisson);
    fprintf(out, "  c      = %.15g\n", par.c);
    fprintf(out, "  phi    = %.15g rad\n", par.phi);
    fprintf(out, "  G      = %.15g\n", par.shear);
    fprintf(out, "  K      = %.15g\n", par.bulk);
    fprintf(out, "  lambda = %.15g\n", par.lame);
    fprintf(out, "  other  = %d components\n", MATMODEL_NCOMP_OTHER);
    fprintf(out, "  eqother= %d components\n", MATMODEL_NCOMP_EQOTHER);
    fprintf(out, "  buffer = epsp(4), return_type(1), proj(12), hess(48), sigma_princ(3)\n");
}

void matmodel::fill_response(const vector &strain,
                             const vector &eqstatev,
                             vector &stress,
                             vector &statev)
{
    const long nstrain = (strain.n < MATMODEL_NCOMP_STRAIN) ? strain.n : static_cast<long>(MATMODEL_NCOMP_STRAIN);
    const long neqother = (eqstatev.n < MATMODEL_NCOMP_EQOTHER) ? eqstatev.n : static_cast<long>(MATMODEL_NCOMP_EQOTHER);
    vector e(MATMODEL_NCOMP_STRAIN), eqp(MATMODEL_NCOMP_EQOTHER), trial(MATMODEL_NCOMP_STRAIN);
    vector eps_el(MATMODEL_NCOMP_STRAIN), eig(3);
    matrix metric(MATMODEL_NCOMP_STRAIN, MATMODEL_NCOMP_STRAIN);
    spectral_data trial_spec, ordered_spec;
    matmodel_cache cache;
    double trace_E, f_tr, gamma_sl, gamma_sr, gamma_la, gamma_ra;
    double denom_s, denom_l, denom_r;
    double lambda_s, lambda_l, lambda_r;
    long idx;

    reallocv(MATMODEL_NCOMP_STRESS, stress);
    reallocv(MATMODEL_NCOMP_OTHER, statev);
    nullv(e);
    nullv(eqp);

    if (nstrain > 0)
        rcopyv(strain, 0, e, 0, nstrain);
    if (neqother > 0)
        rcopyv(eqstatev, 0, eqp, 0, neqother);

    null_cache(cache);
    copyv(eqp, cache.epsp);
    addmultv(e, 1.0, cache.epsp, -1.0, trial);

    nullv(trial_spec.eig);
    nullm(metric);
    for (long i = 0; i < 3; i++)
    {
        nullv(trial_spec.proj[i]);
        nullm(trial_spec.hess[i]);
        nullv(ordered_spec.proj[i]);
        nullm(ordered_spec.hess[i]);
    }
    metric(0, 0) = 1.0;
    metric(1, 1) = 1.0;
    metric(2, 2) = 0.5;

    {
        const double ex = trial[0];
        const double ey = trial[1];
        const double g = trial[2];
        const double d = ex - ey;
        const double i1 = ex + ey;
        const double i2 = sqrt(d * d + g * g);
        const double tol = 1.0e-14;

        trial_spec.eig[0] = 0.5 * (i1 + i2);
        trial_spec.eig[1] = 0.5 * (i1 - i2);
        trial_spec.eig[2] = trial[3];

        if (i2 <= tol)
        {
            trial_spec.proj[0][0] = 1.0;
            trial_spec.proj[0][1] = 1.0;
        }
        else
        {
            trial_spec.proj[0][0] = (ex - trial_spec.eig[1]) / i2;
            trial_spec.proj[0][1] = (ey - trial_spec.eig[1]) / i2;
            trial_spec.proj[0][2] = g / (2.0 * i2);

            trial_spec.proj[1][0] = 1.0 - trial_spec.proj[0][0];
            trial_spec.proj[1][1] = 1.0 - trial_spec.proj[0][1];
            trial_spec.proj[1][2] = -trial_spec.proj[0][2];

            for (long a = 0; a < MATMODEL_NCOMP_STRAIN; a++)
            {
                for (long b = 0; b < MATMODEL_NCOMP_STRAIN; b++)
                {
                    trial_spec.hess[0](a, b) = (metric(a, b) - trial_spec.proj[0][a] * trial_spec.proj[0][b] - trial_spec.proj[1][a] * trial_spec.proj[1][b]) / i2;
                    trial_spec.hess[1](a, b) = -trial_spec.hess[0](a, b);
                }
            }

            for (long a = 0; a < MATMODEL_NCOMP_STRAIN; a++)
            {
                trial_spec.hess[0](3, a) = 0.0;
                trial_spec.hess[0](a, 3) = 0.0;
                trial_spec.hess[1](3, a) = 0.0;
                trial_spec.hess[1](a, 3) = 0.0;
            }
        }

        trial_spec.proj[2][3] = 1.0;
    }

    copy_mode(ordered_spec.proj[0], ordered_spec.hess[0], eig[0], trial_spec.proj[0], trial_spec.hess[0], trial_spec.eig[0]);
    copy_mode(ordered_spec.proj[1], ordered_spec.hess[1], eig[1], trial_spec.proj[1], trial_spec.hess[1], trial_spec.eig[1]);
    copy_mode(ordered_spec.proj[2], ordered_spec.hess[2], eig[2], trial_spec.proj[2], trial_spec.hess[2], trial_spec.eig[2]);

    if ((trial_spec.eig[0] >= trial_spec.eig[2]) && (trial_spec.eig[2] > trial_spec.eig[1]))
    {
        copy_mode(ordered_spec.proj[1], ordered_spec.hess[1], eig[1], trial_spec.proj[2], trial_spec.hess[2], trial_spec.eig[2]);
        copy_mode(ordered_spec.proj[2], ordered_spec.hess[2], eig[2], trial_spec.proj[1], trial_spec.hess[1], trial_spec.eig[1]);
    }
    if (trial_spec.eig[2] > trial_spec.eig[0])
    {
        copy_mode(ordered_spec.proj[0], ordered_spec.hess[0], eig[0], trial_spec.proj[2], trial_spec.hess[2], trial_spec.eig[2]);
        copy_mode(ordered_spec.proj[1], ordered_spec.hess[1], eig[1], trial_spec.proj[0], trial_spec.hess[0], trial_spec.eig[0]);
        copy_mode(ordered_spec.proj[2], ordered_spec.hess[2], eig[2], trial_spec.proj[1], trial_spec.hess[1], trial_spec.eig[1]);
    }

    for (long i = 0; i < 3; i++)
    {
        copyv(ordered_spec.proj[i], cache.proj[i]);
        copym(ordered_spec.hess[i], cache.hess[i]);
    }

    trace_E = eig[0] + eig[1] + eig[2];
    f_tr = 2.0 * par.shear * ((1.0 + par.sin_phi) * eig[0] - (1.0 - par.sin_phi) * eig[2]) + 2.0 * par.lame * par.sin_phi * trace_E - par.c_bar;
    gamma_sl = (eig[0] - eig[1]) / (1.0 + par.sin_phi);
    gamma_sr = (eig[1] - eig[2]) / (1.0 - par.sin_phi);
    gamma_la = (eig[0] + eig[1] - 2.0 * eig[2]) / (3.0 - par.sin_phi);
    gamma_ra = (2.0 * eig[0] - eig[1] - eig[2]) / (3.0 + par.sin_phi);

    denom_s = 4.0 * par.lame * par.sin_phi * par.sin_phi + 2.0 * par.shear * (1.0 + par.sin_phi) * (1.0 + par.sin_phi) + 2.0 * par.shear * (1.0 - par.sin_phi) * (1.0 - par.sin_phi);
    denom_l = 4.0 * par.lame * par.sin_phi * par.sin_phi + par.shear * (1.0 + par.sin_phi) * (1.0 + par.sin_phi) + 2.0 * par.shear * (1.0 - par.sin_phi) * (1.0 - par.sin_phi);
    denom_r = 4.0 * par.lame * par.sin_phi * par.sin_phi + 2.0 * par.shear * (1.0 + par.sin_phi) * (1.0 + par.sin_phi) + par.shear * (1.0 - par.sin_phi) * (1.0 - par.sin_phi);

    lambda_s = f_tr / denom_s;
    lambda_l = (par.shear * ((1.0 + par.sin_phi) * (eig[0] + eig[1]) - 2.0 * (1.0 - par.sin_phi) * eig[2]) + 2.0 * par.lame * par.sin_phi * trace_E - par.c_bar) / denom_l;
    lambda_r = (par.shear * (2.0 * (1.0 + par.sin_phi) * eig[0] - (1.0 - par.sin_phi) * (eig[1] + eig[2])) + 2.0 * par.lame * par.sin_phi * trace_E - par.c_bar) / denom_r;

    if (f_tr <= 0.0)
    {
        cache.return_type = MATMODEL_RETURN_ELASTIC;
        cache.sig_princ[0] = par.lame * trace_E + 2.0 * par.shear * eig[0];
        cache.sig_princ[1] = par.lame * trace_E + 2.0 * par.shear * eig[1];
        cache.sig_princ[2] = par.lame * trace_E + 2.0 * par.shear * eig[2];
    }
    else if (lambda_s <= ((gamma_sl < gamma_sr) ? gamma_sl : gamma_sr))
    {
        cache.return_type = MATMODEL_RETURN_SMOOTH;
        cache.sig_princ[0] = par.lame * trace_E + 2.0 * par.shear * eig[0] - lambda_s * (2.0 * par.lame * par.sin_phi + 2.0 * par.shear * (1.0 + par.sin_phi));
        cache.sig_princ[1] = par.lame * trace_E + 2.0 * par.shear * eig[1] - lambda_s * (2.0 * par.lame * par.sin_phi);
        cache.sig_princ[2] = par.lame * trace_E + 2.0 * par.shear * eig[2] - lambda_s * (2.0 * par.lame * par.sin_phi - 2.0 * par.shear * (1.0 - par.sin_phi));
    }
    else if ((gamma_sl < gamma_sr) && (lambda_l >= gamma_sl) && (lambda_l <= gamma_la))
    {
        cache.return_type = MATMODEL_RETURN_LEFT_EDGE;
        cache.sig_princ[0] = par.lame * trace_E + par.shear * (eig[0] + eig[1]) - lambda_l * (2.0 * par.lame * par.sin_phi + par.shear * (1.0 + par.sin_phi));
        cache.sig_princ[1] = cache.sig_princ[0];
        cache.sig_princ[2] = par.lame * trace_E + 2.0 * par.shear * eig[2] - lambda_l * (2.0 * par.lame * par.sin_phi - 2.0 * par.shear * (1.0 - par.sin_phi));
    }
    else if ((gamma_sl > gamma_sr) && (lambda_r >= gamma_sr) && (lambda_r <= gamma_ra))
    {
        cache.return_type = MATMODEL_RETURN_RIGHT_EDGE;
        cache.sig_princ[0] = par.lame * trace_E + 2.0 * par.shear * eig[0] - lambda_r * (2.0 * par.lame * par.sin_phi + 2.0 * par.shear * (1.0 + par.sin_phi));
        cache.sig_princ[2] = par.lame * trace_E + par.shear * (eig[1] + eig[2]) - lambda_r * (2.0 * par.lame * par.sin_phi - par.shear * (1.0 - par.sin_phi));
        cache.sig_princ[1] = cache.sig_princ[2];
    }
    else
    {
        cache.return_type = MATMODEL_RETURN_APEX;
        cache.sig_princ[0] = par.c_bar / (2.0 * par.sin_phi);
        cache.sig_princ[1] = cache.sig_princ[0];
        cache.sig_princ[2] = cache.sig_princ[0];
    }

    nullv(cache.stress);
    addmultv(cache.stress, cache.proj[0], cache.sig_princ[0]);
    addmultv(cache.stress, cache.proj[1], cache.sig_princ[1]);
    addmultv(cache.stress, cache.proj[2], cache.sig_princ[2]);

    eps_el[0] = (cache.stress[0] - par.poisson * (cache.stress[1] + cache.stress[3])) / par.young;
    eps_el[1] = (cache.stress[1] - par.poisson * (cache.stress[0] + cache.stress[3])) / par.young;
    eps_el[2] = cache.stress[2] / par.shear;
    eps_el[3] = (cache.stress[3] - par.poisson * (cache.stress[0] + cache.stress[1])) / par.young;
    addmultv(e, 1.0, eps_el, -1.0, cache.epsp);

    nullv(statev);
    copyv(cache.epsp, &statev[MATMODEL_IO_EP_XX]);
    statev[MATMODEL_IO_RETURN_TYPE] = static_cast<double>(cache.return_type);
    copyv(cache.proj[0], &statev[MATMODEL_IO_PROJ_1]);
    copyv(cache.proj[1], &statev[MATMODEL_IO_PROJ_2]);
    copyv(cache.proj[2], &statev[MATMODEL_IO_PROJ_3]);

    idx = MATMODEL_IO_HESS_1;
    for (long k = 0; k < 3; k++)
    {
        for (long j = 0; j < MATMODEL_NCOMP_STRAIN; j++)
        {
            for (long i = 0; i < MATMODEL_NCOMP_STRAIN; i++)
                statev[idx++] = cache.hess[k](i, j);
        }
    }
    copyv(cache.sig_princ, &statev[MATMODEL_IO_SIGMA_1]);
    copyv(cache.stress, stress);
    copyv(e, cached_strain);
    copyv(eqp, cached_eqother);
    copyv(statev, cached_other);
    cached_response_valid = 1;
}

void matmodel::nlstresses(const vector &strain,
                          const vector &eqstatev,
                          vector &stress,
                          vector &statev)
{
    fill_response(strain, eqstatev, stress, statev);
}

void matmodel::stiffmat(const vector &strain,
                        const vector &eqstatev,
                        const vector &stress,
                        matrix &d)
{
    matmodel_cache cache;
    vector statev;
    vector hstress;
    vector iota(4), aux(4), proj0(4), proj1(4), proj2(4), proj12(4), proj23(4);
    matrix zero(4, 4), outer(4, 4), hess0(4, 4), hess1(4, 4), hess2(4, 4), hess12(4, 4), hess23(4, 4);
    const double tol = 1.0e-14;
    double denom_s, denom_l, denom_r;
    int match = 1;
    long idx;
    (void)stress;

    if (!cached_response_valid)
        match = 0;

    for (long i = 0; match && (i < MATMODEL_NCOMP_STRAIN); i++)
    {
        const double value = (i < strain.n) ? strain[i] : 0.0;
        if (fabs(value - cached_strain[i]) > tol)
            match = 0;
    }

    for (long i = 0; match && (i < MATMODEL_NCOMP_EQOTHER); i++)
    {
        const double value = (i < eqstatev.n) ? eqstatev[i] : 0.0;
        if (fabs(value - cached_eqother[i]) > tol)
            match = 0;
    }

    if (match)
    {
        reallocv(MATMODEL_NCOMP_OTHER, statev);
        copyv(cached_other, statev);
    }
    else
    {
        fill_response(strain, eqstatev, hstress, statev);
    }

    null_cache(cache);
    if (statev.n > MATMODEL_IO_EP_ZZ)
        rcopyv(statev, MATMODEL_IO_EP_XX, cache.epsp, 0, MATMODEL_NCOMP_EQOTHER);
    if (statev.n > MATMODEL_IO_RETURN_TYPE)
        cache.return_type = static_cast<int>(statev[MATMODEL_IO_RETURN_TYPE] + 0.5);
    if (statev.n > MATMODEL_IO_PROJ_3 + MATMODEL_NCOMP_STRAIN - 1)
    {
        rcopyv(statev, MATMODEL_IO_PROJ_1, cache.proj[0], 0, MATMODEL_NCOMP_STRAIN);
        rcopyv(statev, MATMODEL_IO_PROJ_2, cache.proj[1], 0, MATMODEL_NCOMP_STRAIN);
        rcopyv(statev, MATMODEL_IO_PROJ_3, cache.proj[2], 0, MATMODEL_NCOMP_STRAIN);
    }

    idx = MATMODEL_IO_HESS_1;
    for (long k = 0; k < 3; k++)
    {
        for (long j = 0; j < MATMODEL_NCOMP_STRAIN; j++)
        {
            for (long i = 0; i < MATMODEL_NCOMP_STRAIN; i++)
            {
                if (idx < statev.n)
                    cache.hess[k](i, j) = statev[idx];
                idx++;
            }
        }
    }

    if (statev.n > MATMODEL_IO_SIGMA_3)
        rcopyv(statev, MATMODEL_IO_SIGMA_1, cache.sig_princ, 0, 3);

    reallocm(MATMODEL_NCOMP_STRAIN, MATMODEL_NCOMP_STRAIN, d);
    nullm(d);
    nullm(zero);
    fillv(0.0, iota);
    iota[0] = 1.0;
    iota[1] = 1.0;
    iota[3] = 1.0;

    copyv(cache.proj[0], proj0);
    copyv(cache.proj[1], proj1);
    copyv(cache.proj[2], proj2);
    copym(cache.hess[0], hess0);
    copym(cache.hess[1], hess1);
    copym(cache.hess[2], hess2);

    denom_s = 4.0 * par.lame * par.sin_phi * par.sin_phi + 2.0 * par.shear * (1.0 + par.sin_phi) * (1.0 + par.sin_phi) + 2.0 * par.shear * (1.0 - par.sin_phi) * (1.0 - par.sin_phi);
    denom_l = 4.0 * par.lame * par.sin_phi * par.sin_phi + par.shear * (1.0 + par.sin_phi) * (1.0 + par.sin_phi) + 2.0 * par.shear * (1.0 - par.sin_phi) * (1.0 - par.sin_phi);
    denom_r = 4.0 * par.lame * par.sin_phi * par.sin_phi + 2.0 * par.shear * (1.0 + par.sin_phi) * (1.0 + par.sin_phi) + par.shear * (1.0 - par.sin_phi) * (1.0 - par.sin_phi);

    switch (cache.return_type)
    {
    case MATMODEL_RETURN_ELASTIC:
        d(0, 0) = par.lame + 2.0 * par.shear;
        d(0, 1) = par.lame;
        d(0, 3) = par.lame;
        d(1, 0) = par.lame;
        d(1, 1) = par.lame + 2.0 * par.shear;
        d(1, 3) = par.lame;
        d(2, 2) = par.shear;
        d(3, 0) = par.lame;
        d(3, 1) = par.lame;
        d(3, 3) = par.lame + 2.0 * par.shear;
        break;

    case MATMODEL_RETURN_SMOOTH:
        addmultm(zero, 0.0, hess0, cache.sig_princ[0], d);
        addmultm(zero, 0.0, hess1, cache.sig_princ[1], d);
        addmultm(zero, 0.0, hess2, cache.sig_princ[2], d);
        vxv(iota, iota, outer);
        addmultm(zero, 0.0, outer, par.lame, d);
        vxv(proj0, proj0, outer);
        addmultm(zero, 0.0, outer, 2.0 * par.shear, d);
        vxv(proj1, proj1, outer);
        addmultm(zero, 0.0, outer, 2.0 * par.shear, d);
        vxv(proj2, proj2, outer);
        addmultm(zero, 0.0, outer, 2.0 * par.shear, d);
        addmultv(proj0, 2.0 * par.shear * (1.0 + par.sin_phi),
                 proj2, -2.0 * par.shear * (1.0 - par.sin_phi), aux);
        addmultv(aux, iota, 2.0 * par.lame * par.sin_phi);
        vxv(aux, aux, outer);
        addmultm(zero, 0.0, outer, -1.0 / denom_s, d);
        break;

    case MATMODEL_RETURN_LEFT_EDGE:
        addmultv(proj0, 1.0, proj1, 1.0, proj12);
        addm(hess0, hess1, hess12);
        addmultm(zero, 0.0, hess12, cache.sig_princ[0], d);
        addmultm(zero, 0.0, hess2, cache.sig_princ[2], d);
        vxv(iota, iota, outer);
        addmultm(zero, 0.0, outer, par.lame, d);
        vxv(proj12, proj12, outer);
        addmultm(zero, 0.0, outer, par.shear, d);
        vxv(proj2, proj2, outer);
        addmultm(zero, 0.0, outer, 2.0 * par.shear, d);
        addmultv(proj12, par.shear * (1.0 + par.sin_phi),
                 proj2, -2.0 * par.shear * (1.0 - par.sin_phi), aux);
        addmultv(aux, iota, 2.0 * par.lame * par.sin_phi);
        vxv(aux, aux, outer);
        addmultm(zero, 0.0, outer, -1.0 / denom_l, d);
        break;

    case MATMODEL_RETURN_RIGHT_EDGE:
        addmultv(proj1, 1.0, proj2, 1.0, proj23);
        addm(hess1, hess2, hess23);
        addmultm(zero, 0.0, hess0, cache.sig_princ[0], d);
        addmultm(zero, 0.0, hess23, cache.sig_princ[2], d);
        vxv(iota, iota, outer);
        addmultm(zero, 0.0, outer, par.lame, d);
        vxv(proj0, proj0, outer);
        addmultm(zero, 0.0, outer, 2.0 * par.shear, d);
        vxv(proj23, proj23, outer);
        addmultm(zero, 0.0, outer, par.shear, d);
        addmultv(proj0, 2.0 * par.shear * (1.0 + par.sin_phi),
                 proj23, -par.shear * (1.0 - par.sin_phi), aux);
        addmultv(aux, iota, 2.0 * par.lame * par.sin_phi);
        vxv(aux, aux, outer);
        addmultm(zero, 0.0, outer, -1.0 / denom_r, d);
        break;

    case MATMODEL_RETURN_APEX:
    default:
        break;
    }
}

void matmodel::updateval(const vector &statev, vector &eqstatev)
{
    const long nstatev = (statev.n < MATMODEL_NCOMP_EQOTHER) ? statev.n : static_cast<long>(MATMODEL_NCOMP_EQOTHER);
    reallocv(MATMODEL_NCOMP_EQOTHER, eqstatev);
    nullv(eqstatev);
    if (nstatev > 0)
        rcopyv(statev, 0, eqstatev, 0, nstatev);
}
