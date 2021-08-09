#define PI 3.14159265358979
#define DtoR (PI / 180)	/* conversion of degrees to radians */

/* Lumme-Bowell (1981) phase function */
/* phase_angle in radians */
double
lumme_bowell(double phase_angle)
{
    double
       Phi1,		/* single-scattering phase function */
       Phi_m,		/* multiple-scattering phase function */
       phase_factor,	/* Moon phase factor */
       Q;		/* multiple-scattering fraction */

    Q = 0.108;		/* multiple-scattering fraction */

    /* TODO: find more precise angle for smooth curve? */
    if (phase_angle / DtoR < 25.07)
	Phi1 = 1.0 - sin(phase_angle)
	      / (0.124 + 1.407 * sin(phase_angle) - 0.758 * square(sin(phase_angle)));
    else
	Phi1 = exp(-3.343 * pow(tan(phase_angle / 2.0), 0.632));

    /* multiple scattering */
    Phi_m = 1.0 / PI * (sin(phase_angle) + (PI - phase_angle) * cos(phase_angle));

    phase_factor = (1.0 - Q) * Phi1 + Q * Phi_m;

    return phase_factor;
}

/* air mass: Kasten and Young (1989); apparent altitude in degrees */
double
calc_X_kasten_young(double altitude)
{
    return 1.0 / (sin(altitude * DtoR) + 0.50572 * pow(altitude + 6.07995, -1.6364));
}
