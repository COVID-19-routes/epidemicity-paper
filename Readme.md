# *SEPIAR* epidemiological model

The MATLAB® code in this repository runs the *SEPIAR* epidemiological model.

*SEPIAR* is an acronym for *S*usceptible-*E*xposed-*P*ost
latent-*I*nfectious-*A*symptomatic-*R*ecovered, the most important
epidemiological compartments to describe the transmission of SARS-CoV-2, the
viral agent of COVID-19. The model is spatially explicit, and can be used to
describe the spread of infection over a geographic landscape with multiple human
communities connected via mobility. Implementations of the *SEPIAR* model have
been validated through applications to the early phase of the COVID-19 pandemic
in Italy (Gatto et al. 2020, <https://doi.org/10.1073/pnas.2004978117>; Bertuzzo
et al. 2020, <https://doi.org/10.1038/s41467-020-18050-2>).

The code in this repository evaluates the reproduction number and epidemicity
index for the *SEPIAR* model, two threshold quantities that determine,
respectively, whether pathogen endemicity and short-term epidemicity are
possible. Short-term epidemicity is the appearance of transient epidemic waves
in the presence of an asymptotically stable disease-free equilibrium that
prevents the pathogen from permanently establishing in the community (Mari et
al. 2017, <https://doi.org/10.1111/2041-210X.12805>; 2018,
<https://doi.org/10.1016/j.jtbi.2018.03.024>; 2019,
<https://doi.org/10.1098/rsos.181517>).

Different variations of these quantities can be evaluated with the code, namely
in the absence/presence of containment measures (basic/control reproduction
number and epidemicity index), or for an ongoing outbreak (effective
reproduction number and epidemicity index). The code also performs numerical
simulations of the *SEPIAR* model.

## Running instructions

Launch `SEPIAR_main.m` to run the model.

The code reads the data (`data_provinces.mat`) for the application of the
*SEPIAR* model to the Italian case study at the scale of second-level
administrative units (provinces and metropolitan cities) and set the parameter
values in the absence of containment measures (after Bertuzzo et al. 2020,
*cit*).

Afterwards, it evaluates the basic reproduction number and epidemicity index, as
well as the control reproduction number and epidemicity index for a given choice
of the control parameters. To that end, the external function `SEPIAR_eigen_t.m`
is called within the main script.

Finally, the model compares the outcomes of the two different set-ups
(absence/presence of controls) by numerically integrating the *SEPIAR* model as
described in the external function `SEPIAR_ode.m`. Parameter values, control
implementations, and simulation details can all be easily defined by the user.

The user can also call the `SEPIAR_eigen_t.m` function to evaluate the effective
reproduction number and epidemicity index, namely by using as inputs the
abundances of susceptible (`St`), exposed (`Et`), post-latent (`Pt`),
symptomatic infectious (`It`), asymptomatic (`At`), and recovered individuals
(`Rt`) at a given point in time during the course of the epidemic.

The code has been tested on MATLAB versions R2020b and runs in around 1.2
seconds on MATLAB Online.

### License

The code is licensed under the MIT license.
