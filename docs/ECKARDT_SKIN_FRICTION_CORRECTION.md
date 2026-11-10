# Eckardt A skin-friction formula correction

The baseline impeller skin-friction calculation now follows the published
Jansen relation reproduced by Khoshkalam, Mojaddam, and Pullen, *Energies*
12(14), 2711 (2019), equations 8–12, DOI: 10.3390/en12142711.

The published specific loss is

    dh_sf = 2 Cf (L_H / d_H) Wbar^2

with

    Wbar = (2 W2 + W_h1 + W_s1) / 4.

The previous implementation used `4 Cf` in the loss equation and an RMS
velocity proxy. Its reported Darcy factor remains `4 Cf`; the specific-loss
equation uses `2 Cf`. The prior RMS/shroud proxy remains available only for
the side-by-side diagnostic scenario.

This corrects the source equation implementation. It does not establish that
the Jansen meanline correlation is sufficient for Eckardt A. Re-run the map
and energy-residual diagnostics before drawing conclusions about the remaining
flow-dependent mismatch. A distributed passage integration would be a separate
model extension and would need its own geometry assumptions and validation.
