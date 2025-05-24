import numpy as np
import streamlit as st
import plotly.graph_objects as go

# --- Parameter ---
kmod = 0.6
fc0k, fc90k, fvk, frs = 21, 2.5, 2, 1
bg, hg, bs, hs, ls = 200, 220, 140, 150, 4000
tv, bza, bzal, lv = 40, 50, 55, 90
mh, mv = 0, 0
lseff = 1.3 * hs
lveff = min(lv, 3 * tv)
tx = 0.5 * tv
kc90 = 1.5

fc0d = fc0k * kmod / 1.3
fc90d = fc90k * kmod / 1.3
fvd = fvk * kmod / 1.3
frsd = frs * kmod / 1.3

# --- Bereich ---
alpha_array = np.arange(1, 91, 1)
nd_array = np.arange(0, 20.1, 0.1)

nachweisarten = {
    "Torsion_Zx": "Schub aus Torsion in ZX-Richtung",
    "Torsion_Zy": "Schub aus Torsion in ZY-Richtung",
    "Druck_Stirn": "Versagen der STirnfläche des Zapfens",
    "Druck_Grund": "Versagen der Grundfläche",
    "Schub_Vorholz": "Versagen des Vorholzes"
}
farben = {
    "Torsion_Zx": "red",
    "Torsion_Zy": "blue",
    "Druck_Stirn": "green",
    "Druck_Grund": "orange",
    "Schub_Vorholz": "purple"
}

def berechne_nachweis(nd_array, alpha_deg):
    a = np.radians(alpha_deg)
    s, c = np.sin(a), np.cos(a)
    t = hs / (2 * s)
    b = 1 - mh * mv

    num = ((s - mh * c) / b) * (t - mv * tx) + c * (ls * s + tx) - s * (ls * c + t)
    den = ((c + mh * s) / b) * (t - mv * tx) - s * (ls * s + tx) - c * (ls * c + t)
    m = num / den

    Hcal = nd_array * c + nd_array * m * s
    Vcal = nd_array * ((s - mh * c) / b) - nd_array * m * ((c + mh * s) / b)
    Nza = -Hcal * c - mh * Hcal * s
    Vza = mh * Hcal * c - Hcal * s
    Mza = mh * Hcal * t - Hcal * tx

    x = 0.5 * lseff
    Nx = Nza - mv * Vcal * c - Vcal * s
    Vx = Vza + Vcal * c - mv * Vcal * s
    Mx = Mza + Vza * x + Vcal * x * c - mv * Vcal * x * s

    dNcal = 0.5 * (Nza - Nx / 3)
    dVcal = 0.5 * (Vza - Vx / 3)
    dMcal = 0.5 * (Mza - Mx / 3 + Vza * 0.5 * lseff) - dVcal * (lseff / 6)

    Wt = (0.18 / 0.8) * hs ** 2 * lseff
    c3 = 0.915
    tzx = np.abs(dNcal / (lseff * hs) * 1000) + np.abs(dMcal / Wt * 1000)
    tzy = np.abs(2 * dVcal / (lseff * hs) * 1000) + np.abs(c3 * dMcal / Wt * 1000)
    n1 = 100 * tzx / fvd
    n2 = 100 * tzy / frsd

    Avorh = tv * bza
    Aeff = (tv + 30 * s) * bza
    kca = 1 + (kc90 - 1) * s
    A = -0.5 - 0.5 * np.cos(2 * a)
    B = -0.5 + 0.5 * np.cos(2 * a)
    C = 0.5 * np.sin(2 * a)
    ka = 1 / np.sqrt(A**2 + (B * fc0k / fc90k)**2 + (C * fc0k / fvk)**2)
    fca = ka * kca * fc0k * Aeff / Avorh
    Hcal2 = Avorh * fca / 1000
    Fsfcal = Hcal2 / (c + m * s - mv * (((s - mh * c) / b) - m * ((c + mh * s) / b)))
    Fsfcald = Fsfcal * kmod / 1.3
    n3 = 100 * nd_array / Fsfcald

    Ageff = (bg - bzal) * (hs / s + 60)
    Vcal2 = Ageff * kc90 * fc90k / 1000
    Fgfcal = Vcal2 / (((s - mh * c) / b) - m * ((c + mh * s) / b))
    Fgfcald = Fgfcal * kmod / 1.3
    n4 = 100 * nd_array / Fgfcald

    Aveff = lveff * bza
    tavd = 1000 * nd_array * c / Aveff
    n5 = 100 * tavd / fvd

    return {
        "Torsion_Zx": n1,
        "Torsion_Zy": n2,
        "Druck_Stirn": n3,
        "Druck_Grund": n4,
        "Schub_Vorholz": n5
    }

# --- Streamlit UI ---
st.set_page_config(layout="wide")
st.title("Nachweisvergleich in Abhängigkeit vom α-Winkel")

selected_keys = st.multiselect(
    "Wähle eine oder mehrere Nachweisarten:",
    options=list(nachweisarten.keys()),
    default=["Torsion_Zx"]
)

alpha = st.slider("Winkel α [°]", min_value=1, max_value=90, value=45)
werte = berechne_nachweis(nd_array, alpha)

fig = go.Figure()
for key in selected_keys:
    fig.add_trace(go.Scatter(
        x=nd_array,
        y=werte[key],
        mode="lines",
        name=nachweisarten[key],
        line=dict(color=farben[key])
    ))

fig.update_layout(
    title=f"Ausnutzung bei α = {alpha}°",
    xaxis_title="Kraft Nd [kN]",
    yaxis_title="Ausnutzung [%]",
    yaxis_range=[0, 150]
)

st.plotly_chart(fig, use_container_width=True)

st.info("Die Kurven zeigen die Ausnutzung der gewählten Nachweise in Abhängigkeit von Nd und dem α-Winkel.")
