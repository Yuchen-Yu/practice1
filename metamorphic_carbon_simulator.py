import time

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import streamlit as st

from matplotlib.path import Path
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D


# ============================================================
# 1. PAGE SETTINGS
# ============================================================

st.set_page_config(
    page_title="Metamorphic carbon simulator",
    page_icon="🪨",
    layout="wide",
)

st.title(
    "Metamorphic carbon reaction–transport simulator"
)

st.caption(
    "A simplified P–T–XCO₂ constrained model of decarbonation, "
    "carbon transport and carbon fixation."
)


# ============================================================
# 2. PUBLICATION-STYLE FIGURE SETTINGS
# ============================================================

plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": [
            "Arial",
            "DejaVu Sans",
        ],
        "mathtext.fontset": "dejavusans",

        "font.size": 8,
        "axes.titlesize": 8,
        "axes.labelsize": 8,

        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "legend.fontsize": 7,

        "axes.linewidth": 0.6,
        "xtick.major.width": 0.5,
        "ytick.major.width": 0.5,

        "lines.linewidth": 1.2,

        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",

        "axes.grid": False,
    }
)


# ============================================================
# 3. COLOR PALETTE
# ============================================================

# Red series
RED_DARK = "#780018"
RED = "#AA0033"
RED_MID = "#CC0044"
RED_LIGHT = "#FA8095"

# Green series
GREEN_DARK = "#25362B"
GREEN = "#276E3D"
GREEN_LIGHT = "#56B76A"

# Grey
GREY = "#8E8E99"

# Brown
BROWN_DARK = "#532F1A"
BROWN_LIGHT = "#BE9A62"

# Blue
BLUE_DARK = "#242F4D"
BLUE = "#6382AA"

# Purple
PURPLE_DARK = "#413A4C"

# Olive
OLIVE = "#5A5745"


# ============================================================
# 4. VISUAL ENCODING
# ============================================================

COLOR_CALCITE = BROWN_LIGHT
COLOR_QUARTZ = GREY
COLOR_WOLLASTONITE = RED_MID

COLOR_CO2 = BLUE

COLOR_FORSTERITE = OLIVE
COLOR_MAGNESITE = GREEN_LIGHT
COLOR_ENSTATITE = GREEN_DARK

COLOR_EXPORT = RED
COLOR_TOTAL = PURPLE_DARK

COLOR_LIGHT_BORDER = "#B9B9B9"


# ------------------------------------------------------------
# Restore the original coloured process-zone backgrounds
# ------------------------------------------------------------

COLOR_ZONE_DECARB = "#F8E9ED"
COLOR_ZONE_TRANSPORT = "#EDF1F6"
COLOR_ZONE_FIXATION = "#EAF3EC"


# ============================================================
# 5. PHYSICAL CONSTANTS
# ============================================================

R_GAS = 8.314462618


# ============================================================
# 6. REFERENCE CONDITIONS
# ============================================================

# Reference fluid compositions associated with the
# experimental/reference curves.

XCO2_REF_DECARB = 0.96
XCO2_REF_FIX = 1.00


# Effective entropy changes used only to estimate how
# the reference P–T boundaries move with XCO2.

DS_DECARB_EFFECTIVE = 160.0
DS_FIX_EFFECTIVE = 170.0


# ============================================================
# 7. REACTION-DRIVING-FORCE SCALES
# ============================================================

# Once a point enters the thermodynamically favourable region,
# the reaction strength increases smoothly with distance
# from the equilibrium boundary.

DEC_DRIVE_SCALE_KBAR = 0.35
FIX_DRIVE_SCALE_KBAR = 1.50


# ============================================================
# 8. DECARBONATION REFERENCE BOUNDARY
# ============================================================

# Approximate digitization of the classical:
#
# Calcite + Quartz
# =
# Wollastonite + CO2
#
# experimental boundary.

_DEC_T_REF = np.array(
    [
        480.0,
        550.0,
        600.0,
        650.0,
        700.0,
        750.0,
        790.0,
        820.0,
    ]
)

_DEC_P_REF = np.array(
    [
        0.10,
        0.30,
        0.55,
        0.95,
        1.55,
        2.35,
        3.05,
        3.70,
    ]
)


_DEC_COEF = np.polyfit(
    _DEC_T_REF,
    _DEC_P_REF,
    3,
)

_DEC_DCOEF = np.polyder(
    _DEC_COEF
)


def decarbonation_boundary_reference(
    temperature_c
):
    """
    Reference WQC boundary near XCO2 = 0.96.

    Reaction:

        Calcite + Quartz
        =
        Wollastonite + CO2

    Pressure returned in kbar.
    """

    T = np.asarray(
        temperature_c,
        dtype=float,
    )

    pressure = np.polyval(
        _DEC_COEF,
        T,
    )

    return np.maximum(
        pressure,
        0.0,
    )


def decarbonation_boundary_dP_dT(
    temperature_c
):
    """
    Local slope dP/dT of the reference
    decarbonation boundary.

    Unit:
        kbar K-1
    """

    T = np.asarray(
        temperature_c,
        dtype=float,
    )

    slope = np.polyval(
        _DEC_DCOEF,
        T,
    )

    return np.maximum(
        slope,
        1.0e-5,
    )


# ============================================================
# 9. CARBONATION REFERENCE BOUNDARY
# ============================================================

def fixation_boundary_reference(
    temperature_c
):
    """
    Koziol & Newton (1998) experimental boundary:

        Magnesite + Enstatite
        =
        Forsterite + CO2

    Carbon fixation is the reverse:

        Forsterite + CO2
        =
        Magnesite + Enstatite

    Published equation:

        P =
        -3.215
        -4.784e-3 T
        +2.542e-5 T^2

    P:
        kbar

    T:
        °C
    """

    T = np.asarray(
        temperature_c,
        dtype=float,
    )

    return (
        -3.215
        - 4.784e-3 * T
        + 2.542e-5 * T**2
    )


def fixation_boundary_dP_dT(
    temperature_c
):
    """
    Local dP/dT slope of the
    Koziol–Newton boundary.
    """

    T = np.asarray(
        temperature_c,
        dtype=float,
    )

    slope = (
        -4.784e-3
        + 2.0
        * 2.542e-5
        * T
    )

    return np.maximum(
        slope,
        1.0e-5,
    )


# ============================================================
# 10. XCO2 EFFECT ON REACTION BOUNDARIES
# ============================================================

def _xco2_pressure_shift(
    temperature_c,
    xco2,
    xco2_reference,
    dP_dT_kbar_per_k,
    ds_effective,
):
    """
    Teaching-level ideal-fluid correction.

    At the reference boundary:

        DeltaG = 0

    Changing XCO2 contributes:

        RT ln(XCO2 / Xref)

    The local reaction volume is approximated
    using the Clapeyron relation:

        dP/dT = DeltaS / DeltaV

    The resulting pressure correction shifts both
    the plotted boundary and the numerical reaction
    gate together.

    This prevents the visual P-T field and the
    numerical model from disagreeing.
    """

    T_k = (
        np.asarray(
            temperature_c,
            dtype=float,
        )
        + 273.15
    )

    x = np.clip(
        float(xco2),
        1.0e-6,
        1.0,
    )

    x_ref = float(
        xco2_reference
    )

    slope_pa_per_k = (
        np.asarray(
            dP_dT_kbar_per_k,
            dtype=float,
        )
        * 1.0e8
    )

    delta_v = (
        ds_effective
        /
        np.maximum(
            slope_pa_per_k,
            1.0e3,
        )
    )

    delta_p_kbar = (
        -R_GAS
        * T_k
        * np.log(
            x / x_ref
        )
        /
        (
            delta_v
            * 1.0e8
        )
    )

    return delta_p_kbar


# ============================================================
# 11. XCO2-ADJUSTED DECARBONATION BOUNDARY
# ============================================================

def decarbonation_boundary(
    temperature_c,
    xco2,
):
    """
    P-T-XCO2 boundary actually used by the model.

    Below this boundary:
        decarbonation is favourable

    Above this boundary:
        decarbonation is disabled.
    """

    base = (
        decarbonation_boundary_reference(
            temperature_c
        )
    )

    shift = (
        _xco2_pressure_shift(
            temperature_c=
            temperature_c,

            xco2=
            xco2,

            xco2_reference=
            XCO2_REF_DECARB,

            dP_dT_kbar_per_k=
            decarbonation_boundary_dP_dT(
                temperature_c
            ),

            ds_effective=
            DS_DECARB_EFFECTIVE,
        )
    )

    return np.maximum(
        base + shift,
        0.0,
    )


# ============================================================
# 12. XCO2-ADJUSTED CARBONATION BOUNDARY
# ============================================================

def fixation_boundary(
    temperature_c,
    xco2,
):
    """
    P-T-XCO2 boundary actually used by the model.

    Above this boundary:
        carbonation is favourable

    Below this boundary:
        carbonation is disabled.
    """

    base = (
        fixation_boundary_reference(
            temperature_c
        )
    )

    shift = (
        _xco2_pressure_shift(
            temperature_c=
            temperature_c,

            xco2=
            xco2,

            xco2_reference=
            XCO2_REF_FIX,

            dP_dT_kbar_per_k=
            fixation_boundary_dP_dT(
                temperature_c
            ),

            ds_effective=
            DS_FIX_EFFECTIVE,
        )
    )

    return np.maximum(
        base + shift,
        0.0,
    )


# ============================================================
# 13. THERMODYNAMIC STATE
# ============================================================

def reaction_thermodynamics(
    temperature_c,
    decarb_pressure_kbar,
    fixation_pressure_kbar,
    xco2,
):
    """
    IMPORTANT:

    The numerical reaction gate uses EXACTLY the
    same boundaries displayed in panel a.

    Decarbonation:
        favourable only if

        P_decarb < P_eq_decarbonation

    Carbon fixation:
        favourable only if

        P_fixation > P_eq_carbonation

    Therefore a point located in the grey field
    cannot continue the corresponding forward reaction.
    """

    p_dec_eq = float(
        decarbonation_boundary(
            temperature_c,
            xco2,
        )
    )

    p_fix_eq = float(
        fixation_boundary(
            temperature_c,
            xco2,
        )
    )


    # --------------------------------------------------------
    # Signed distance from reaction boundaries
    # --------------------------------------------------------

    dec_margin = (
        p_dec_eq
        - decarb_pressure_kbar
    )

    fix_margin = (
        fixation_pressure_kbar
        - p_fix_eq
    )


    # --------------------------------------------------------
    # Exact thermodynamic gates
    # --------------------------------------------------------

    dec_favourable = (
        dec_margin > 0.0
    )

    fix_favourable = (
        fix_margin > 0.0
    )


    # --------------------------------------------------------
    # Reaction strength inside favourable region
    #
    # Outside:
    #     exactly 0
    # --------------------------------------------------------

    if dec_favourable:

        dec_drive = (
            1.0
            - np.exp(
                -dec_margin
                /
                DEC_DRIVE_SCALE_KBAR
            )
        )

    else:

        dec_drive = 0.0


    if fix_favourable:

        fix_drive = (
            1.0
            - np.exp(
                -fix_margin
                /
                FIX_DRIVE_SCALE_KBAR
            )
        )

    else:

        fix_drive = 0.0


    # ========================================================
    # Approximate signed Delta G for display
    # ========================================================

    dec_slope_pa_per_k = (
        float(
            decarbonation_boundary_dP_dT(
                temperature_c
            )
        )
        * 1.0e8
    )

    dec_dv = (
        DS_DECARB_EFFECTIVE
        /
        max(
            dec_slope_pa_per_k,
            1.0e3,
        )
    )

    dg_dec_kj = (
        dec_dv
        * (
            decarb_pressure_kbar
            - p_dec_eq
        )
        * 1.0e8
        / 1000.0
    )


    fix_slope_pa_per_k = (
        float(
            fixation_boundary_dP_dT(
                temperature_c
            )
        )
        * 1.0e8
    )

    fix_dv = (
        DS_FIX_EFFECTIVE
        /
        max(
            fix_slope_pa_per_k,
            1.0e3,
        )
    )

    # For carbonation:
    # pressure above equilibrium gives negative DeltaG.

    dg_fix_kj = (
        -fix_dv
        * (
            fixation_pressure_kbar
            - p_fix_eq
        )
        * 1.0e8
        / 1000.0
    )


    return {

        "p_dec_eq":
            p_dec_eq,

        "p_fix_eq":
            p_fix_eq,

        "dec_margin_kbar":
            dec_margin,

        "fix_margin_kbar":
            fix_margin,

        "decarbonation_favourable":
            dec_favourable,

        "carbonation_favourable":
            fix_favourable,

        "decarbonation_drive":
            dec_drive,

        "carbonation_drive":
            fix_drive,

        "dg_decarbonation":
            dg_dec_kj,

        "dg_carbonation":
            dg_fix_kj,
    }


# ============================================================
# 14. ARRHENIUS KINETICS
# ============================================================

def arrhenius_rate(
    k_reference,
    temperature_c,
    reference_temperature_c,
    activation_energy_kj,
):
    """
    Arrhenius scaling:

        k(T) =
        k_ref exp[
            -Ea/R
            (1/T - 1/Tref)
        ]
    """

    T = (
        temperature_c
        + 273.15
    )

    T_ref = (
        reference_temperature_c
        + 273.15
    )

    Ea = (
        activation_energy_kj
        * 1000.0
    )

    exponent = (
        -Ea
        / R_GAS
        * (
            1.0 / T
            - 1.0 / T_ref
        )
    )

    exponent = np.clip(
        exponent,
        -50,
        50,
    )

    return (
        k_reference
        * np.exp(
            exponent
        )
    )


# ============================================================
# 15. EFFECTIVE REACTION RATES
# ============================================================

def calculate_reaction_rates(
    temperature_c,
    decarb_pressure_kbar,
    fixation_pressure_kbar,
    xco2,
    k_dec_ref,
    k_fix_ref,
):
    """
    Model logic:

        P-T-XCO2
             ↓
        equilibrium boundary
             ↓
        favourable / unfavourable
             ↓
        thermodynamic driving force
             ↓
        Arrhenius kinetics
             ↓
        effective reaction rate

    Critical rule:

        grey field → rate = 0
    """

    thermo = (
        reaction_thermodynamics(
            temperature_c=
            temperature_c,

            decarb_pressure_kbar=
            decarb_pressure_kbar,

            fixation_pressure_kbar=
            fixation_pressure_kbar,

            xco2=
            xco2,
        )
    )


    k_dec_kinetic = (
        arrhenius_rate(
            k_reference=
            k_dec_ref,

            temperature_c=
            temperature_c,

            reference_temperature_c=
            800.0,

            activation_energy_kj=
            100.0,
        )
    )


    k_fix_kinetic = (
        arrhenius_rate(
            k_reference=
            k_fix_ref,

            temperature_c=
            temperature_c,

            reference_temperature_c=
            800.0,

            activation_energy_kj=
            80.0,
        )
    )


    # --------------------------------------------------------
    # THIS IS THE CRITICAL FIX
    # --------------------------------------------------------

    k_dec_effective = (
        k_dec_kinetic
        * thermo[
            "decarbonation_drive"
        ]
    )

    k_fix_effective = (
        k_fix_kinetic
        * thermo[
            "carbonation_drive"
        ]
    )


    return {

        **thermo,

        "kinetic_decarbonation":
            k_dec_kinetic,

        "kinetic_carbonation":
            k_fix_kinetic,

        "effective_decarbonation":
            k_dec_effective,

        "effective_carbonation":
            k_fix_effective,
    }


# ============================================================
# 16. INITIAL 1-D MODEL STATE
# ============================================================

def build_initial_state(
    nx=180,
    length_km=30.0,
    initial_carbon=100.0,
    fixation_capacity_fraction=0.70,
):

    x = np.linspace(
        0.0,
        length_km,
        nx,
    )

    dx = (
        x[1]
        - x[0]
    )


    source_mask = (
        x <= 9.0
    )

    transport_mask = (
        (x > 9.0)
        &
        (x < 20.0)
    )

    fixation_mask = (
        x >= 20.0
    )


    S = np.zeros(
        nx
    )

    C = np.zeros(
        nx
    )

    F = np.zeros(
        nx
    )


    # --------------------------------------------------------
    # Initial carbonate carbon
    # --------------------------------------------------------

    source_length = (
        source_mask.sum()
        * dx
    )

    S[
        source_mask
    ] = (
        initial_carbon
        /
        source_length
    )


    # --------------------------------------------------------
    # Fixation capacity
    # --------------------------------------------------------

    total_fix_capacity = (
        initial_carbon
        * fixation_capacity_fraction
    )

    fixation_capacity = (
        np.zeros(
            nx
        )
    )

    fixation_length = (
        fixation_mask.sum()
        * dx
    )

    if fixation_length > 0:

        fixation_capacity[
            fixation_mask
        ] = (
            total_fix_capacity
            /
            fixation_length
        )


    masks = {

        "source":
            source_mask,

        "transport":
            transport_mask,

        "fixation":
            fixation_mask,
    }


    return (
        x,
        dx,
        S,
        C,
        F,
        fixation_capacity,
        masks,
    )


# ============================================================
# 17. NUMERICAL SIMULATION
# ============================================================

@st.cache_data(
    show_spinner=False
)
def simulate(
    temperature_c=750.0,
    decarb_pressure_kbar=1.0,
    fixation_pressure_kbar=15.0,
    xco2=0.20,
    velocity_km_per_kyr=0.18,
    diffusivity_km2_per_kyr=0.05,
    k_dec_ref=0.045,
    k_fix_ref=0.060,
    fixation_capacity_fraction=0.70,
    total_time_kyr=100.0,
    n_output=181,
    nx=180,
    length_km=30.0,
):

    rates = (
        calculate_reaction_rates(
            temperature_c=
            temperature_c,

            decarb_pressure_kbar=
            decarb_pressure_kbar,

            fixation_pressure_kbar=
            fixation_pressure_kbar,

            xco2=
            xco2,

            k_dec_ref=
            k_dec_ref,

            k_fix_ref=
            k_fix_ref,
        )
    )


    k_dec = (
        rates[
            "effective_decarbonation"
        ]
    )

    k_fix = (
        rates[
            "effective_carbonation"
        ]
    )


    (
        x,
        dx,
        S,
        C,
        F,
        fix_capacity,
        masks,
    ) = build_initial_state(
        nx=
        nx,

        length_km=
        length_km,

        initial_carbon=
        100.0,

        fixation_capacity_fraction=
        fixation_capacity_fraction,
    )


    source = (
        masks[
            "source"
        ]
    )

    fixation = (
        masks[
            "fixation"
        ]
    )


    # ========================================================
    # NUMERICAL TIME STEP
    # ========================================================

    dt_limits = [
        0.20
    ]


    if velocity_km_per_kyr > 0:

        dt_limits.append(
            0.40
            * dx
            /
            velocity_km_per_kyr
        )


    if diffusivity_km2_per_kyr > 0:

        dt_limits.append(
            0.20
            * dx**2
            /
            diffusivity_km2_per_kyr
        )


    if k_dec > 0:

        dt_limits.append(
            0.10
            /
            k_dec
        )


    if k_fix > 0:

        dt_limits.append(
            0.10
            /
            k_fix
        )


    dt = max(
        1.0e-5,
        min(
            dt_limits
        ),
    )


    n_steps = int(
        np.ceil(
            total_time_kyr
            /
            dt
        )
    )


    dt = (
        total_time_kyr
        /
        n_steps
    )


    output_times = (
        np.linspace(
            0.0,
            total_time_kyr,
            n_output,
        )
    )


    history = {

        "time": [],

        "solid_original": [],

        "fluid": [],

        "fixed": [],

        "exported": [],

        "total_accounted": [],

        "C_profile": [],

        "S_profile": [],

        "F_profile": [],
    }


    exported = 0.0


    initial_total = (
        (
            S.sum()
            + C.sum()
            + F.sum()
        )
        * dx
    )


    # ========================================================
    # RECORD MODEL STATE
    # ========================================================

    def record(
        t
    ):

        solid = (
            S.sum()
            * dx
        )

        fluid = (
            C.sum()
            * dx
        )

        fixed = (
            F.sum()
            * dx
        )

        total = (
            solid
            + fluid
            + fixed
            + exported
        )


        history[
            "time"
        ].append(
            t
        )

        history[
            "solid_original"
        ].append(
            solid
        )

        history[
            "fluid"
        ].append(
            fluid
        )

        history[
            "fixed"
        ].append(
            fixed
        )

        history[
            "exported"
        ].append(
            exported
        )

        history[
            "total_accounted"
        ].append(
            total
        )

        history[
            "C_profile"
        ].append(
            C.copy()
        )

        history[
            "S_profile"
        ].append(
            S.copy()
        )

        history[
            "F_profile"
        ].append(
            F.copy()
        )


    record(
        0.0
    )


    next_output = 1


    # ========================================================
    # TIME INTEGRATION
    # ========================================================

    for step in range(
        1,
        n_steps + 1,
    ):

        # ====================================================
        # A. DECARBONATION
        # ====================================================

        release = (
            np.zeros_like(
                S
            )
        )


        if k_dec > 0:

            release[
                source
            ] = (
                S[
                    source
                ]
                * (
                    1.0
                    - np.exp(
                        -k_dec
                        * dt
                    )
                )
            )


        S -= release
        C += release


        # ====================================================
        # B. CARBON FIXATION
        # ====================================================

        remaining_capacity = (
            np.maximum(
                fix_capacity
                - F,
                0.0,
            )
        )


        capture = (
            np.zeros_like(
                C
            )
        )


        if k_fix > 0:

            kinetic_capture = (
                C[
                    fixation
                ]
                * (
                    1.0
                    - np.exp(
                        -k_fix
                        * dt
                    )
                )
            )


            capture[
                fixation
            ] = (
                np.minimum(
                    kinetic_capture,

                    remaining_capacity[
                        fixation
                    ],
                )
            )


        C -= capture
        F += capture


        # ====================================================
        # C. ADVECTION
        # ====================================================

        if velocity_km_per_kyr > 0:

            courant = (
                velocity_km_per_kyr
                * dt
                / dx
            )


            move = (
                courant
                * C
            )


            C -= move


            C[
                1:
            ] += (
                move[
                    :-1
                ]
            )


            exported += (
                move[
                    -1
                ]
                * dx
            )


        # ====================================================
        # D. DIFFUSION
        # ====================================================

        if diffusivity_km2_per_kyr > 0:

            r = (
                diffusivity_km2_per_kyr
                * dt
                / dx**2
            )


            transfer = (
                r
                * (
                    C[
                        :-1
                    ]
                    -
                    C[
                        1:
                    ]
                )
            )


            C[
                :-1
            ] -= (
                transfer
            )


            C[
                1:
            ] += (
                transfer
            )


        # ----------------------------------------------------
        # Numerical cleanup
        # ----------------------------------------------------

        C[
            C < 0
        ] = 0.0

        S[
            S < 0
        ] = 0.0

        F[
            F < 0
        ] = 0.0


        t = (
            step
            * dt
        )


        while (
            next_output
            < n_output

            and

            t + 1.0e-12
            >=
            output_times[
                next_output
            ]
        ):

            record(
                output_times[
                    next_output
                ]
            )

            next_output += 1


    while (
        len(
            history[
                "time"
            ]
        )
        <
        n_output
    ):

        record(
            output_times[
                len(
                    history[
                        "time"
                    ]
                )
            ]
        )


    for key in (
        "C_profile",
        "S_profile",
        "F_profile",
    ):

        history[
            key
        ] = (
            np.asarray(
                history[
                    key
                ]
            )
        )


    for key in (
        "time",
        "solid_original",
        "fluid",
        "fixed",
        "exported",
        "total_accounted",
    ):

        history[
            key
        ] = (
            np.asarray(
                history[
                    key
                ]
            )
        )


    metadata = {

        "x":
            x,

        "dx":
            dx,

        "masks":
            masks,

        "initial_total":
            initial_total,

        "fix_capacity_total":
            fix_capacity.sum()
            * dx,

        "rates":
            rates,

        "temperature_c":
            temperature_c,

        "decarb_pressure_kbar":
            decarb_pressure_kbar,

        "fixation_pressure_kbar":
            fixation_pressure_kbar,

        "xco2":
            xco2,
    }


    return (
        history,
        metadata,
    )


# ============================================================
# 18. MINERAL SYMBOLS
# ============================================================

CALCITE_MARKER = (
    Path(
        np.array(
            [
                [-0.80, -0.30],
                [-0.25, -0.85],
                [0.80, -0.30],
                [0.25, 0.85],
                [-0.80, -0.30],
            ]
        )
    )
)


WOLLASTONITE_MARKER = (
    Path(
        np.array(
            [
                [-0.17, -1.00],
                [0.17, -1.00],
                [0.32, -0.15],
                [0.19, 1.00],
                [-0.19, 1.00],
                [-0.32, -0.15],
                [-0.17, -1.00],
            ]
        )
    )
)


MAGNESITE_MARKER = (
    Path(
        np.array(
            [
                [-0.72, 0.00],
                [0.00, -0.82],
                [0.72, 0.00],
                [0.00, 0.82],
                [-0.72, 0.00],
            ]
        )
    )
)


# ============================================================
# 19. FIGURE UTILITIES
# ============================================================

def clean_axes(
    ax
):

    ax.spines[
        "top"
    ].set_visible(
        False
    )

    ax.spines[
        "right"
    ].set_visible(
        False
    )

    ax.spines[
        "left"
    ].set_linewidth(
        0.6
    )

    ax.spines[
        "bottom"
    ].set_linewidth(
        0.6
    )

    ax.tick_params(
        direction=
        "out",

        length=
        3,

        width=
        0.5,
    )


def crystal_positions(
    xmin,
    xmax,
    n=20,
    seed=1,
):

    rng = (
        np.random.default_rng(
            seed
        )
    )


    x = (
        rng.uniform(
            xmin,
            xmax,
            n,
        )
    )


    y = (
        rng.uniform(
            0.34,
            0.64,
            n,
        )
    )


    return (
        x,
        y,
    )


# ============================================================
# 20. DECARBONATION MINERALS
# ============================================================

def draw_decarbonation_minerals(
    ax,
    progress,
):

    n = 20


    x, y = (
        crystal_positions(
            0.85,
            8.10,

            n=
            n,

            seed=
            100,
        )
    )


    rng = (
        np.random.default_rng(
            101
        )
    )


    order = (
        rng.permutation(
            n
        )
    )


    n_reacted = int(
        np.round(
            progress
            * n
        )
    )


    reacted = (
        order[
            :n_reacted
        ]
    )


    unreacted = (
        order[
            n_reacted:
        ]
    )


    # --------------------------------------------------------
    # Calcite + Quartz
    # --------------------------------------------------------

    if len(
        unreacted
    ) > 0:

        ax.scatter(
            x[
                unreacted
            ] - 0.10,

            y[
                unreacted
            ],

            s=
            72,

            marker=
            CALCITE_MARKER,

            facecolor=
            COLOR_CALCITE,

            edgecolor=
            "white",

            linewidth=
            0.45,

            zorder=
            5,
        )


        ax.scatter(
            x[
                unreacted
            ] + 0.13,

            y[
                unreacted
            ] + 0.018,

            s=
            24,

            marker=
            "D",

            facecolor=
            COLOR_QUARTZ,

            edgecolor=
            "white",

            linewidth=
            0.35,

            zorder=
            6,
        )


    # --------------------------------------------------------
    # Wollastonite
    # --------------------------------------------------------

    if len(
        reacted
    ) > 0:

        ax.scatter(
            x[
                reacted
            ],

            y[
                reacted
            ],

            s=
            90,

            marker=
            WOLLASTONITE_MARKER,

            facecolor=
            COLOR_WOLLASTONITE,

            edgecolor=
            "white",

            linewidth=
            0.45,

            zorder=
            5,
        )


# ============================================================
# 21. FIXATION MINERALS
# ============================================================

def draw_fixation_minerals(
    ax,
    progress,
):

    n = 20


    x, y = (
        crystal_positions(
            20.85,
            28.75,

            n=
            n,

            seed=
            200,
        )
    )


    rng = (
        np.random.default_rng(
            201
        )
    )


    order = (
        rng.permutation(
            n
        )
    )


    n_reacted = int(
        np.round(
            progress
            * n
        )
    )


    reacted = (
        order[
            :n_reacted
        ]
    )


    unreacted = (
        order[
            n_reacted:
        ]
    )


    # --------------------------------------------------------
    # Forsterite
    # --------------------------------------------------------

    if len(
        unreacted
    ) > 0:

        ax.scatter(
            x[
                unreacted
            ],

            y[
                unreacted
            ],

            s=
            72,

            marker=
            "h",

            facecolor=
            COLOR_FORSTERITE,

            edgecolor=
            "white",

            linewidth=
            0.45,

            zorder=
            5,
        )


    # --------------------------------------------------------
    # Magnesite + Enstatite
    # --------------------------------------------------------

    if len(
        reacted
    ) > 0:

        ax.scatter(
            x[
                reacted
            ] - 0.09,

            y[
                reacted
            ],

            s=
            66,

            marker=
            MAGNESITE_MARKER,

            facecolor=
            COLOR_MAGNESITE,

            edgecolor=
            "white",

            linewidth=
            0.45,

            zorder=
            5,
        )


        ax.scatter(
            x[
                reacted
            ] + 0.13,

            y[
                reacted
            ] + 0.018,

            s=
            27,

            marker=
            "^",

            facecolor=
            COLOR_ENSTATITE,

            edgecolor=
            "white",

            linewidth=
            0.35,

            zorder=
            6,
        )


# ============================================================
# 22. REACTION PROGRESS BAR
# ============================================================

def draw_progress_bar(
    ax,
    x0,
    x1,
    progress,
    color,
    label,
):

    progress = float(
        np.clip(
            progress,
            0.0,
            1.0,
        )
    )


    y = 0.145
    height = 0.018


    ax.add_patch(
        Rectangle(
            (
                x0,
                y,
            ),

            x1 - x0,
            height,

            facecolor=
            "white",

            edgecolor=
            COLOR_LIGHT_BORDER,

            linewidth=
            0.5,

            zorder=
            2,
        )
    )


    ax.add_patch(
        Rectangle(
            (
                x0,
                y,
            ),

            (
                x1 - x0
            )
            * progress,

            height,

            facecolor=
            color,

            edgecolor=
            "none",

            zorder=
            3,
        )
    )


    ax.text(
        (
            x0 + x1
        )
        / 2,

        y - 0.025,

        (
            f"{label} "
            f"{progress * 100:.0f}%"
        ),

        ha=
        "center",

        va=
        "top",

        fontsize=
        7,

        color=
        color,
    )


# ============================================================
# 23. CO2 PARTICLES
# ============================================================

def sample_co2_particles(
    history,
    meta,
    idx,
    max_particles=130,
):

    C = (
        history[
            "C_profile"
        ][
            idx
        ]
    )


    x = (
        meta[
            "x"
        ]
    )


    dx = (
        meta[
            "dx"
        ]
    )


    fluid_inventory = (
        history[
            "fluid"
        ][
            idx
        ]
    )


    initial_total = (
        meta[
            "initial_total"
        ]
    )


    fraction = (
        fluid_inventory
        /
        max(
            initial_total,
            1.0e-12,
        )
    )


    n_particles = int(
        np.clip(
            np.round(
                max_particles
                * fraction
                * 3.2
            ),
            0,
            max_particles,
        )
    )


    if (
        n_particles == 0

        or

        np.sum(
            C
        ) <= 0
    ):

        return (
            np.array([]),
            np.array([]),
        )


    probability = (
        C
        /
        np.sum(
            C
        )
    )


    rng = (
        np.random.default_rng(
            5000
            + idx
        )
    )


    selected = (
        rng.choice(
            len(
                x
            ),

            size=
            n_particles,

            replace=
            True,

            p=
            probability,
        )
    )


    particle_x = (
        x[
            selected
        ]

        + rng.uniform(
            -0.40
            * dx,

            0.40
            * dx,

            n_particles,
        )
    )


    particle_y = (
        rng.uniform(
            0.34,
            0.64,
            n_particles,
        )
    )


    return (
        particle_x,
        particle_y,
    )


# ============================================================
# 24. MAIN COMPOSITE FIGURE
#
#       a                  b
#   P-T fields      carbon profiles
#
#              c
#        process model
#
# ============================================================

def make_process_figure(
    history,
    meta,
    idx,
):

    # ========================================================
    # MODEL STATE
    # ========================================================

    x = (
        meta[
            "x"
        ]
    )


    S = (
        history[
            "S_profile"
        ][
            idx
        ]
    )


    C = (
        history[
            "C_profile"
        ][
            idx
        ]
    )


    F = (
        history[
            "F_profile"
        ][
            idx
        ]
    )


    t = (
        history[
            "time"
        ][
            idx
        ]
    )


    temperature_c = (
        meta[
            "temperature_c"
        ]
    )


    decarb_pressure = (
        meta[
            "decarb_pressure_kbar"
        ]
    )


    fixation_pressure = (
        meta[
            "fixation_pressure_kbar"
        ]
    )


    xco2 = (
        meta[
            "xco2"
        ]
    )


    rates = (
        meta[
            "rates"
        ]
    )


    # ========================================================
    # REACTION PROGRESS
    # ========================================================

    solid_initial = (
        history[
            "solid_original"
        ][
            0
        ]
    )


    solid_current = (
        history[
            "solid_original"
        ][
            idx
        ]
    )


    fixed_current = (
        history[
            "fixed"
        ][
            idx
        ]
    )


    fixation_capacity = (
        meta[
            "fix_capacity_total"
        ]
    )


    decarb_progress = (
        1.0
        -
        solid_current
        /
        max(
            solid_initial,
            1.0e-12,
        )
    )


    decarb_progress = float(
        np.clip(
            decarb_progress,
            0.0,
            1.0,
        )
    )


    if fixation_capacity > 0:

        fixation_progress = (
            fixed_current
            /
            fixation_capacity
        )

    else:

        fixation_progress = 0.0


    fixation_progress = float(
        np.clip(
            fixation_progress,
            0.0,
            1.0,
        )
    )


    # ========================================================
    # FIGURE
    # ========================================================

    fig = (
        plt.figure(
            figsize=(
                7.1,
                6.15,
            ),

            facecolor=
            "white",
        )
    )


    gs = (
        fig.add_gridspec(
            nrows=
            2,

            ncols=
            2,

            height_ratios=[
                1.00,
                1.12,
            ],

            width_ratios=[
                1.00,
                1.00,
            ],

            hspace=
            0.42,

            wspace=
            0.30,
        )
    )


    ax_pt = (
        fig.add_subplot(
            gs[
                0,
                0
            ]
        )
    )


    ax_curve = (
        fig.add_subplot(
            gs[
                0,
                1
            ]
        )
    )


    ax_process = (
        fig.add_subplot(
            gs[
                1,
                :
            ]
        )
    )


    # ========================================================
    # PANEL A — P-T-XCO2 REACTION FIELDS
    # ========================================================

    T_MIN = 480.0
    T_MAX = 1120.0

    P_MIN = 0.0
    P_MAX = 26.0


    T_all = (
        np.linspace(
            T_MIN,
            T_MAX,
            700,
        )
    )


    # --------------------------------------------------------
    # IMPORTANT:
    #
    # The same two functions control BOTH:
    #
    # 1. the plotted boundaries
    # 2. the actual numerical reaction gates
    # --------------------------------------------------------

    P_dec_all = (
        decarbonation_boundary(
            T_all,
            xco2,
        )
    )


    P_fix_all = (
        fixation_boundary(
            T_all,
            xco2,
        )
    )


    P_dec_fill = (
        np.clip(
            P_dec_all,
            P_MIN,
            P_MAX,
        )
    )


    P_fix_fill = (
        np.clip(
            P_fix_all,
            P_MIN,
            P_MAX,
        )
    )


    # --------------------------------------------------------
    # GREY BASE
    # --------------------------------------------------------

    ax_pt.fill_between(
        T_all,

        P_MIN,
        P_MAX,

        color=
        GREY,

        alpha=
        0.50,

        linewidth=
        0,

        zorder=
        0,
    )


    # --------------------------------------------------------
    # RED — DECARBONATION FAVOURABLE
    # --------------------------------------------------------

    ax_pt.fill_between(
        T_all,

        P_MIN,

        P_dec_fill,

        color=
        RED_LIGHT,

        alpha=
        0.50,

        linewidth=
        0,

        zorder=
        1,
    )


    # --------------------------------------------------------
    # GREEN — CARBONATION FAVOURABLE
    # --------------------------------------------------------

    ax_pt.fill_between(
        T_all,

        P_fix_fill,

        P_MAX,

        color=
        GREEN_LIGHT,

        alpha=
        0.50,

        linewidth=
        0,

        zorder=
        1,
    )


    # ========================================================
    # BOUNDARY CURVES
    # ========================================================

    # --------------------------------------------------------
    # Full decarbonation curve
    # --------------------------------------------------------

    ax_pt.plot(
        T_all,
        P_dec_all,

        color=
        RED,

        linewidth=
        1.10,

        linestyle=
        "--",

        alpha=
        0.90,

        zorder=
        3,
    )


    # --------------------------------------------------------
    # Full carbonation curve
    # --------------------------------------------------------

    ax_pt.plot(
        T_all,
        P_fix_all,

        color=
        GREEN,

        linewidth=
        1.10,

        linestyle=
        "--",

        alpha=
        0.90,

        zorder=
        3,
    )


    # --------------------------------------------------------
    # Reference experimental temperature ranges
    # --------------------------------------------------------

    mask_dec_valid = (
        (T_all >= 480.0)
        &
        (T_all <= 820.0)
    )


    mask_fix_valid = (
        (T_all >= 700.0)
        &
        (T_all <= 1100.0)
    )


    ax_pt.plot(
        T_all[
            mask_dec_valid
        ],

        P_dec_all[
            mask_dec_valid
        ],

        color=
        RED,

        linewidth=
        1.8,

        linestyle=
        "-",

        label=
        "Decarbonation boundary",

        zorder=
        4,
    )


    ax_pt.plot(
        T_all[
            mask_fix_valid
        ],

        P_fix_all[
            mask_fix_valid
        ],

        color=
        GREEN,

        linewidth=
        1.8,

        linestyle=
        "-",

        label=
        "Carbonation boundary",

        zorder=
        4,
    )


    # ========================================================
    # CURRENT MODEL P-T POINTS
    # ========================================================

    ax_pt.scatter(
        temperature_c,
        decarb_pressure,

        s=
        46,

        marker=
        "o",

        facecolor=
        RED,

        edgecolor=
        "white",

        linewidth=
        0.70,

        zorder=
        10,

        label=
        "Decarbonation P–T",
    )


    ax_pt.scatter(
        temperature_c,
        fixation_pressure,

        s=
        48,

        marker=
        "s",

        facecolor=
        GREEN,

        edgecolor=
        "white",

        linewidth=
        0.70,

        zorder=
        10,

        label=
        "Fixation P–T",
    )


    ax_pt.axvline(
        temperature_c,

        color=
        "#666666",

        linestyle=
        "--",

        linewidth=
        0.65,

        alpha=
        0.70,

        zorder=
        2,
    )


    # ========================================================
    # FIELD LABELS
    # ========================================================

    ax_pt.text(
        690,
        0.75,

        "Decarbonation\nfavourable",

        ha=
        "center",

        va=
        "center",

        fontsize=
        6.2,

        fontweight=
        "bold",

        color=
        RED_DARK,

        zorder=
        5,
    )


    ax_pt.text(
        755,
        7.5,

        "Forward reactions\n"
        "thermodynamically\n"
        "unfavourable",

        ha=
        "center",

        va=
        "center",

        fontsize=
        5.8,

        color=
        "#404040",

        zorder=
        5,
    )


    ax_pt.text(
        795,
        23.0,

        "Carbonation favourable",

        ha=
        "center",

        va=
        "center",

        fontsize=
        6.2,

        fontweight=
        "bold",

        color=
        GREEN_DARK,

        zorder=
        5,
    )


    ax_pt.text(
        800,
        21.4,

        "Magnesite + enstatite",

        ha=
        "center",

        fontsize=
        5.8,

        color=
        GREEN_DARK,
    )


    ax_pt.text(
        1010,
        8.0,

        r"Forsterite + CO$_2$",

        ha=
        "center",

        fontsize=
        5.8,

        color=
        GREEN_DARK,
    )


    ax_pt.text(
        720,
        2.5,

        r"Wollastonite + CO$_2$",

        ha=
        "center",

        fontsize=
        5.8,

        color=
        RED_DARK,
    )


    ax_pt.text(
        1110,
        25.3,

        (
            rf"$X_{{CO_2}}$ = "
            f"{xco2:.2f}"
        ),

        ha=
        "right",

        va=
        "top",

        fontsize=
        6.1,

        color=
        BLUE_DARK,
    )


    ax_pt.text(
        1110,
        0.55,

        "Dashed = T extrapolation",

        ha=
        "right",

        va=
        "bottom",

        fontsize=
        5.4,

        color=
        "#555555",
    )


    # ========================================================
    # AXES
    # ========================================================

    ax_pt.set_xlim(
        T_MIN,
        T_MAX,
    )


    ax_pt.set_ylim(
        P_MIN,
        P_MAX,
    )


    ax_pt.set_xlabel(
        "Temperature (°C)"
    )


    ax_pt.set_ylabel(
        "Pressure (kbar)"
    )


    ax_pt.set_title(
        "P–T reaction boundaries",
        pad=
        6,
    )


    ax_pt.text(
        -0.14,
        1.05,

        "a",

        transform=
        ax_pt.transAxes,

        fontsize=
        8,

        fontweight=
        "bold",

        va=
        "bottom",
    )


    ax_pt.legend(
        frameon=
        False,

        fontsize=
        5.6,

        loc=
        "upper left",

        handlelength=
        1.8,

        labelspacing=
        0.40,
    )


    clean_axes(
        ax_pt
    )


    # ========================================================
    # PANEL B — SPATIAL CARBON CONCENTRATIONS
    # ========================================================

    ax_curve.plot(
        x,
        S,

        color=
        BROWN_DARK,

        linewidth=
        1.4,

        label=
        "Original solid carbon",
    )


    ax_curve.plot(
        x,
        C,

        color=
        BLUE,

        linewidth=
        1.4,

        label=
        "Mobile fluid carbon",
    )


    ax_curve.plot(
        x,
        F,

        color=
        GREEN,

        linewidth=
        1.4,

        label=
        "Newly fixed carbon",
    )


    ax_curve.axvline(
        9,

        color=
        GREY,

        linestyle=
        "--",

        linewidth=
        0.6,
    )


    ax_curve.axvline(
        20,

        color=
        GREY,

        linestyle=
        "--",

        linewidth=
        0.6,
    )


    ymax = max(
        np.max(
            S
        ),

        np.max(
            C
        ),

        np.max(
            F
        ),

        1.0,
    )


    ax_curve.set_xlim(
        0,
        30,
    )


    ax_curve.set_ylim(
        0,
        ymax
        * 1.20,
    )


    ax_curve.text(
        4.5,

        ymax
        * 1.10,

        "Decarbonation",

        ha=
        "center",

        fontsize=
        6.2,

        color=
        RED_DARK,
    )


    ax_curve.text(
        14.5,

        ymax
        * 1.10,

        "Transport",

        ha=
        "center",

        fontsize=
        6.2,

        color=
        BLUE_DARK,
    )


    ax_curve.text(
        25.0,

        ymax
        * 1.10,

        "Fixation",

        ha=
        "center",

        fontsize=
        6.2,

        color=
        GREEN_DARK,
    )


    ax_curve.set_xlabel(
        "Distance along metamorphic system (km)"
    )


    ax_curve.set_ylabel(
        "Carbon concentration\n"
        "(units km$^{-1}$)"
    )


    ax_curve.set_title(
        (
            "Spatial carbon distribution "
            f"at {t:.1f} kyr"
        ),

        pad=
        6,
    )


    ax_curve.text(
        -0.14,
        1.05,

        "b",

        transform=
        ax_curve.transAxes,

        fontsize=
        8,

        fontweight=
        "bold",

        va=
        "bottom",
    )


    ax_curve.legend(
        frameon=
        False,

        fontsize=
        6.1,

        loc=
        "upper right",
    )


    clean_axes(
        ax_curve
    )


    # ========================================================
    # PANEL C — PROCESS VISUALIZATION
    # ========================================================

    ax_process.set_xlim(
        0,
        30,
    )


    ax_process.set_ylim(
        0,
        1,
    )


    ax_process.axis(
        "off"
    )


    zone_y = 0.26
    zone_h = 0.44


    # --------------------------------------------------------
    # RESTORED COLOURED BACKGROUNDS
    # --------------------------------------------------------

    ax_process.add_patch(
        Rectangle(
            (
                0.25,
                zone_y,
            ),

            8.50,
            zone_h,

            facecolor=
            COLOR_ZONE_DECARB,

            edgecolor=
            RED_DARK,

            linewidth=
            0.8,
        )
    )


    ax_process.add_patch(
        Rectangle(
            (
                9.10,
                zone_y,
            ),

            10.55,
            zone_h,

            facecolor=
            COLOR_ZONE_TRANSPORT,

            edgecolor=
            BLUE_DARK,

            linewidth=
            0.8,
        )
    )


    ax_process.add_patch(
        Rectangle(
            (
                20.00,
                zone_y,
            ),

            9.60,
            zone_h,

            facecolor=
            COLOR_ZONE_FIXATION,

            edgecolor=
            GREEN_DARK,

            linewidth=
            0.8,
        )
    )


    # ========================================================
    # TITLES
    # ========================================================

    ax_process.text(
        4.5,
        0.94,

        "Decarbonation",

        ha=
        "center",

        fontsize=
        8,

        fontweight=
        "bold",

        color=
        RED_DARK,
    )


    ax_process.text(
        14.4,
        0.94,

        "Carbon transport",

        ha=
        "center",

        fontsize=
        8,

        fontweight=
        "bold",

        color=
        BLUE_DARK,
    )


    ax_process.text(
        24.8,
        0.94,

        "Carbon fixation",

        ha=
        "center",

        fontsize=
        8,

        fontweight=
        "bold",

        color=
        GREEN_DARK,
    )


    # ========================================================
    # CHEMICAL REACTIONS
    # ========================================================

    ax_process.text(
        4.5,
        0.83,

        r"$\mathrm{CaCO_3+SiO_2"
        r"\rightarrow CaSiO_3+CO_2}$",

        ha=
        "center",

        fontsize=
        7,
    )


    ax_process.text(
        14.4,
        0.83,

        "Advection + diffusion",

        ha=
        "center",

        fontsize=
        7,
    )


    ax_process.text(
        24.8,
        0.83,

        r"$\mathrm{Mg_2SiO_4+CO_2"
        r"\rightarrow MgCO_3+MgSiO_3}$",

        ha=
        "center",

        fontsize=
        7,
    )


    # ========================================================
    # THERMODYNAMIC STATUS
    # ========================================================

    ax_process.text(
        4.5,
        0.745,

        (
            "Thermodynamically favourable"
            if rates[
                "decarbonation_favourable"
            ]
            else
            "Thermodynamically unfavourable"
        ),

        ha=
        "center",

        fontsize=
        5.8,

        color=(
            RED_DARK
            if rates[
                "decarbonation_favourable"
            ]
            else
            "#555555"
        ),
    )


    ax_process.text(
        24.8,
        0.745,

        (
            "Thermodynamically favourable"
            if rates[
                "carbonation_favourable"
            ]
            else
            "Thermodynamically unfavourable"
        ),

        ha=
        "center",

        fontsize=
        5.8,

        color=(
            GREEN_DARK
            if rates[
                "carbonation_favourable"
            ]
            else
            "#555555"
        ),
    )


    # ========================================================
    # MINERAL TRANSFORMATION
    # ========================================================

    draw_decarbonation_minerals(
        ax_process,
        decarb_progress,
    )


    draw_fixation_minerals(
        ax_process,
        fixation_progress,
    )


    # ========================================================
    # FLUID PATHWAY
    # ========================================================

    ax_process.plot(
        [
            7.8,
            29.0,
        ],

        [
            0.49,
            0.49,
        ],

        color=
        BLUE_DARK,

        linewidth=
        0.65,

        zorder=
        2,
    )


    for arrow_x in [
        9.6,
        12.4,
        15.2,
        18.0,
        21.5,
        24.5,
        27.0,
    ]:

        ax_process.annotate(
            "",

            xy=(
                arrow_x
                + 0.8,

                0.49,
            ),

            xytext=(
                arrow_x,
                0.49,
            ),

            arrowprops=
            dict(
                arrowstyle=
                "-|>",

                color=
                BLUE_DARK,

                lw=
                0.65,

                mutation_scale=
                8,
            ),
        )


    # ========================================================
    # CO2 PARTICLES
    # ========================================================

    (
        particle_x,
        particle_y,
    ) = (
        sample_co2_particles(
            history,
            meta,
            idx,
        )
    )


    if len(
        particle_x
    ) > 0:

        ax_process.scatter(
            particle_x,
            particle_y,

            s=
            14,

            marker=
            "o",

            facecolor=
            COLOR_CO2,

            edgecolor=
            "white",

            linewidth=
            0.30,

            alpha=
            0.92,

            zorder=
            7,
        )


    # ========================================================
    # EXPORT
    # ========================================================

    ax_process.annotate(
        "",

        xy=(
            30.0,
            0.49,
        ),

        xytext=(
            29.0,
            0.49,
        ),

        arrowprops=
        dict(
            arrowstyle=
            "-|>",

            color=
            COLOR_EXPORT,

            lw=
            1.0,

            mutation_scale=
            10,
        ),
    )


    ax_process.text(
        29.65,
        0.39,

        "Export",

        ha=
        "right",

        fontsize=
        7,

        color=
        RED_DARK,
    )


    # ========================================================
    # PROGRESS BARS
    # ========================================================

    draw_progress_bar(
        ax_process,

        0.90,
        8.10,

        decarb_progress,

        RED,

        "Reaction progress",
    )


    draw_progress_bar(
        ax_process,

        20.70,
        28.90,

        fixation_progress,

        GREEN,

        "Fixation progress",
    )


    ax_process.text(
        14.4,
        0.18,

        r"CO$_2$-bearing metamorphic fluid",

        ha=
        "center",

        fontsize=
        7,

        color=
        BLUE_DARK,
    )


    ax_process.text(
        29.65,
        0.985,

        f"{t:.1f} kyr",

        ha=
        "right",

        va=
        "top",

        fontsize=
        8,

        fontweight=
        "bold",
    )


    ax_process.text(
        0.0,
        0.995,

        "c",

        fontsize=
        8,

        fontweight=
        "bold",

        va=
        "top",
    )


    # ========================================================
    # LEGEND
    # ========================================================

    legend_handles = [

        Line2D(
            [0],
            [0],

            marker=
            CALCITE_MARKER,

            linestyle=
            "",

            markerfacecolor=
            COLOR_CALCITE,

            markeredgecolor=
            "white",

            markersize=
            6,

            label=
            "Calcite",
        ),


        Line2D(
            [0],
            [0],

            marker=
            "D",

            linestyle=
            "",

            markerfacecolor=
            COLOR_QUARTZ,

            markeredgecolor=
            "white",

            markersize=
            4,

            label=
            "Quartz",
        ),


        Line2D(
            [0],
            [0],

            marker=
            WOLLASTONITE_MARKER,

            linestyle=
            "",

            markerfacecolor=
            COLOR_WOLLASTONITE,

            markeredgecolor=
            "white",

            markersize=
            7,

            label=
            "Wollastonite",
        ),


        Line2D(
            [0],
            [0],

            marker=
            "o",

            linestyle=
            "",

            markerfacecolor=
            COLOR_CO2,

            markeredgecolor=
            "white",

            markersize=
            5,

            label=
            r"CO$_2$",
        ),


        Line2D(
            [0],
            [0],

            marker=
            "h",

            linestyle=
            "",

            markerfacecolor=
            COLOR_FORSTERITE,

            markeredgecolor=
            "white",

            markersize=
            6,

            label=
            "Forsterite",
        ),


        Line2D(
            [0],
            [0],

            marker=
            MAGNESITE_MARKER,

            linestyle=
            "",

            markerfacecolor=
            COLOR_MAGNESITE,

            markeredgecolor=
            "white",

            markersize=
            6,

            label=
            "Magnesite",
        ),


        Line2D(
            [0],
            [0],

            marker=
            "^",

            linestyle=
            "",

            markerfacecolor=
            COLOR_ENSTATITE,

            markeredgecolor=
            "white",

            markersize=
            5,

            label=
            "Enstatite",
        ),
    ]


    ax_process.legend(
        handles=
        legend_handles,

        frameon=
        False,

        ncol=
        7,

        loc=
        "lower center",

        bbox_to_anchor=(
            0.5,
            0.005,
        ),

        fontsize=
        6,

        columnspacing=
        0.65,

        handletextpad=
        0.25,
    )


    fig.subplots_adjust(
        left=
        0.075,

        right=
        0.985,

        top=
        0.96,

        bottom=
        0.055,
    )


    return fig


# ============================================================
# 25. CARBON RESERVOIRS THROUGH TIME
# ============================================================

def make_reservoir_figure(
    history,
    current_idx,
):

    fig, ax = (
        plt.subplots(
            figsize=(
                5.3,
                3.4,
            )
        )
    )


    t = (
        history[
            "time"
        ]
    )


    ax.plot(
        t,

        history[
            "solid_original"
        ],

        color=
        BROWN_DARK,

        linewidth=
        1.4,

        label=
        "Original solid carbon",
    )


    ax.plot(
        t,

        history[
            "fluid"
        ],

        color=
        BLUE,

        linewidth=
        1.4,

        label=
        "Mobile fluid carbon",
    )


    ax.plot(
        t,

        history[
            "fixed"
        ],

        color=
        GREEN,

        linewidth=
        1.4,

        label=
        "Newly fixed carbon",
    )


    ax.plot(
        t,

        history[
            "exported"
        ],

        color=
        RED,

        linewidth=
        1.4,

        label=
        "Exported carbon",
    )


    ax.axvline(
        t[
            current_idx
        ],

        color=
        GREY,

        linestyle=
        "--",

        linewidth=
        0.7,
    )


    ax.set_xlabel(
        "Time (kyr)"
    )


    ax.set_ylabel(
        "Carbon inventory (units)"
    )


    ax.set_title(
        "Carbon reservoirs through time",
        pad=
        6,
    )


    ax.text(
        -0.10,
        1.04,

        "d",

        transform=
        ax.transAxes,

        fontsize=
        8,

        fontweight=
        "bold",
    )


    ax.legend(
        frameon=
        False,

        loc=
        "upper right",
    )


    clean_axes(
        ax
    )


    fig.tight_layout()


    return fig


# ============================================================
# 26. MASS BALANCE FIGURE
# ============================================================

def make_balance_figure(
    history
):

    fig, ax = (
        plt.subplots(
            figsize=(
                5.3,
                3.4,
            )
        )
    )


    t = (
        history[
            "time"
        ]
    )


    total = (
        history[
            "total_accounted"
        ]
    )


    initial = (
        total[
            0
        ]
    )


    ax.plot(
        t,
        total,

        color=
        COLOR_TOTAL,

        linewidth=
        1.5,

        label=
        "Accounted carbon",
    )


    ax.axhline(
        initial,

        color=
        GREY,

        linestyle=
        "--",

        linewidth=
        0.8,

        label=
        "Initial carbon",
    )


    ax.set_xlabel(
        "Time (kyr)"
    )


    ax.set_ylabel(
        "Total carbon (units)"
    )


    ax.set_title(
        "Mass-balance check",
        pad=
        6,
    )


    ax.text(
        -0.10,
        1.04,

        "e",

        transform=
        ax.transAxes,

        fontsize=
        8,

        fontweight=
        "bold",
    )


    ax.legend(
        frameon=
        False
    )


    clean_axes(
        ax
    )


    fig.tight_layout()


    return fig


# ============================================================
# 27. SIDEBAR CONTROLS
# ============================================================

st.sidebar.header(
    "Model parameters"
)


temperature_c = (
    st.sidebar.slider(
        "Metamorphic temperature (°C)",
        500,
        1050,
        750,
        10,
    )
)


decarb_pressure = (
    st.sidebar.slider(
        "Decarbonation-zone pressure (kbar)",
        0.1,
        15.0,
        1.0,
        0.1,
    )
)


fixation_pressure = (
    st.sidebar.slider(
        "Fixation-zone pressure (kbar)",
        0.1,
        25.0,
        15.0,
        0.5,
    )
)


xco2 = (
    st.sidebar.slider(
        "Fluid XCO₂",
        0.05,
        1.00,
        0.20,
        0.01,
        format=
        "%.2f",
    )
)


velocity = (
    st.sidebar.slider(
        "Fluid velocity (km kyr⁻¹)",
        0.00,
        0.50,
        0.18,
        0.01,
    )
)


diffusivity = (
    st.sidebar.slider(
        "Effective diffusivity (km² kyr⁻¹)",
        0.00,
        0.20,
        0.05,
        0.01,
    )
)


k_decarb_ref = (
    st.sidebar.slider(
        "Reference decarbonation rate (kyr⁻¹)",
        0.005,
        0.100,
        0.045,
        0.005,
        format=
        "%.3f",
    )
)


k_fix_ref = (
    st.sidebar.slider(
        "Reference fixation rate (kyr⁻¹)",
        0.005,
        0.120,
        0.060,
        0.005,
        format=
        "%.3f",
    )
)


fix_capacity_pct = (
    st.sidebar.slider(
        "Maximum fixation capacity (% of initial carbon)",
        0,
        120,
        70,
        5,
    )
)


total_time = (
    st.sidebar.slider(
        "Simulation duration (kyr)",
        20,
        200,
        100,
        10,
    )
)


# ============================================================
# 28. RUN MODEL
# ============================================================

history, meta = (
    simulate(

        temperature_c=
        float(
            temperature_c
        ),

        decarb_pressure_kbar=
        float(
            decarb_pressure
        ),

        fixation_pressure_kbar=
        float(
            fixation_pressure
        ),

        xco2=
        float(
            xco2
        ),

        velocity_km_per_kyr=
        float(
            velocity
        ),

        diffusivity_km2_per_kyr=
        float(
            diffusivity
        ),

        k_dec_ref=
        float(
            k_decarb_ref
        ),

        k_fix_ref=
        float(
            k_fix_ref
        ),

        fixation_capacity_fraction=(
            float(
                fix_capacity_pct
            )
            / 100.0
        ),

        total_time_kyr=
        float(
            total_time
        ),
    )
)


rates = (
    meta[
        "rates"
    ]
)


# ============================================================
# 29. CORE PROCESSES
# ============================================================

st.subheader(
    "1. Core processes"
)


c1, c2, c3 = (
    st.columns(
        3
    )
)


with c1:

    st.markdown(
        "**Decarbonation**"
    )

    st.latex(
        r"""
        \mathrm{
        CaCO_3+SiO_2
        \rightarrow
        CaSiO_3+CO_2
        }
        """
    )


with c2:

    st.markdown(
        "**Carbon transport**"
    )

    st.latex(
        r"""
        \frac{\partial C}{\partial t}
        +
        v\frac{\partial C}{\partial x}
        =
        D\frac{\partial^2 C}{\partial x^2}
        +
        R_{\rm dec}
        -
        R_{\rm fix}
        """
    )


with c3:

    st.markdown(
        "**Carbon fixation**"
    )

    st.latex(
        r"""
        \mathrm{
        Mg_2SiO_4+CO_2
        \rightarrow
        MgCO_3+MgSiO_3
        }
        """
    )


# ============================================================
# 30. THERMODYNAMIC AND KINETIC STATE
# ============================================================

st.subheader(
    "2. Thermodynamic and kinetic state"
)


t1, t2, t3, t4 = (
    st.columns(
        4
    )
)


if rates[
    "decarbonation_favourable"
]:

    dec_state = (
        "Favourable"
    )

else:

    dec_state = (
        "Unfavourable"
    )


if rates[
    "carbonation_favourable"
]:

    fix_state = (
        "Favourable"
    )

else:

    fix_state = (
        "Unfavourable"
    )


t1.metric(
    "Decarbonation ΔG",

    (
        f"{rates['dg_decarbonation']:.2f} "
        "kJ mol⁻¹"
    ),

    dec_state,
)


t2.metric(
    "Carbonation ΔG",

    (
        f"{rates['dg_carbonation']:.2f} "
        "kJ mol⁻¹"
    ),

    fix_state,
)


t3.metric(
    "Effective decarbonation rate",

    (
        f"{rates['effective_decarbonation']:.4f} "
        "kyr⁻¹"
    ),
)


t4.metric(
    "Effective fixation rate",

    (
        f"{rates['effective_carbonation']:.4f} "
        "kyr⁻¹"
    ),
)


# ============================================================
# 31. DYNAMIC SYSTEM
# ============================================================

st.subheader(
    "3. Dynamic metamorphic carbon system"
)


time_idx = (
    st.slider(
        "Simulation time",

        min_value=
        0,

        max_value=(
            len(
                history[
                    "time"
                ]
            )
            - 1
        ),

        value=
        min(
            90,

            len(
                history[
                    "time"
                ]
            )
            - 1,
        ),

        step=
        1,
    )
)


# ============================================================
# 32. CURRENT CARBON STATE
# ============================================================

solid_now = (
    history[
        "solid_original"
    ][
        time_idx
    ]
)


fluid_now = (
    history[
        "fluid"
    ][
        time_idx
    ]
)


fixed_now = (
    history[
        "fixed"
    ][
        time_idx
    ]
)


exported_now = (
    history[
        "exported"
    ][
        time_idx
    ]
)


m1, m2, m3, m4 = (
    st.columns(
        4
    )
)


m1.metric(
    "Original solid carbon",
    f"{solid_now:.2f}",
)


m2.metric(
    "Mobile fluid carbon",
    f"{fluid_now:.2f}",
)


m3.metric(
    "Newly fixed carbon",
    f"{fixed_now:.2f}",
)


m4.metric(
    "Exported carbon",
    f"{exported_now:.2f}",
)


# ============================================================
# 33. MAIN COMPOSITE FIGURE
# ============================================================

fig_process = (
    make_process_figure(
        history,
        meta,
        time_idx,
    )
)


st.pyplot(
    fig_process,
    use_container_width=
    True,
)


plt.close(
    fig_process
)


# ============================================================
# 34. ANIMATION
# ============================================================

control_1, control_2 = (
    st.columns(
        [
            1,
            4,
        ]
    )
)


with control_1:

    animate = (
        st.button(
            "▶ Run simulation",
            type=
            "primary",
        )
    )


with control_2:

    animation_speed = (
        st.select_slider(
            "Animation speed",

            options=[
                "Slow",
                "Normal",
                "Fast",
            ],

            value=
            "Normal",
        )
    )


speed_map = {

    "Slow":
        0.12,

    "Normal":
        0.06,

    "Fast":
        0.025,
}


animation_placeholder = (
    st.empty()
)


if animate:

    frame_ids = (
        np.linspace(
            0,

            len(
                history[
                    "time"
                ]
            )
            - 1,

            min(
                80,

                len(
                    history[
                        "time"
                    ]
                ),
            ),

            dtype=
            int,
        )
    )


    for i in frame_ids:

        fig_animation = (
            make_process_figure(
                history,
                meta,
                int(
                    i
                ),
            )
        )


        animation_placeholder.pyplot(
            fig_animation,
            use_container_width=
            True,
        )


        plt.close(
            fig_animation
        )


        time.sleep(
            speed_map[
                animation_speed
            ]
        )


# ============================================================
# 35. QUANTITATIVE RESULTS
# ============================================================

st.subheader(
    "4. Carbon reservoirs and mass balance"
)


left, right = (
    st.columns(
        2
    )
)


with left:

    fig_reservoir = (
        make_reservoir_figure(
            history,
            time_idx,
        )
    )


    st.pyplot(
        fig_reservoir,
        use_container_width=
        True,
    )


    plt.close(
        fig_reservoir
    )


with right:

    fig_balance = (
        make_balance_figure(
            history
        )
    )


    st.pyplot(
        fig_balance,
        use_container_width=
        True,
    )


    plt.close(
        fig_balance
    )


# ============================================================
# 36. MASS BALANCE CHECK
# ============================================================

initial_total = (
    meta[
        "initial_total"
    ]
)


final_total = (
    history[
        "total_accounted"
    ][
        -1
    ]
)


relative_error = (
    abs(
        final_total
        - initial_total
    )
    /
    max(
        initial_total,
        1.0e-12,
    )
)


if relative_error < 1.0e-8:

    st.success(
        f"Mass balance conserved | "
        f"Initial carbon = {initial_total:.6f} | "
        f"Final accounted carbon = {final_total:.6f} | "
        f"Relative error = {relative_error:.2e}"
    )


else:

    st.warning(
        f"Numerical mass-balance error = "
        f"{relative_error:.2e}"
    )


# ============================================================
# 37. FINAL CARBON PARTITIONING
# ============================================================

st.subheader(
    "5. Final carbon partitioning"
)


summary = (
    pd.DataFrame(
        {
            "Reservoir": [

                "Original solid carbon",

                "Mobile fluid carbon",

                "Newly fixed carbon",

                "Exported carbon",
            ],

            "Final inventory": [

                history[
                    "solid_original"
                ][
                    -1
                ],

                history[
                    "fluid"
                ][
                    -1
                ],

                history[
                    "fixed"
                ][
                    -1
                ],

                history[
                    "exported"
                ][
                    -1
                ],
            ],
        }
    )
)


summary[
    "Fraction of initial carbon (%)"
] = (
    summary[
        "Final inventory"
    ]
    /
    initial_total
    * 100.0
)


st.dataframe(
    summary,
    use_container_width=
    True,

    hide_index=
    True,
)