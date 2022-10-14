### Summary of Data Products:

Simulated data files are calculated for the full 2016 COSI balloon campaign: <br />

start time: 1463443400.0 s <br />
stop time: 1467475400.0 s <br />
total time: 4032000.0 s = 46.67 days <br />

mass model = COSIBalloon.9Detector.geo.setup <br />

simulated energy range: 100 keV - 10 MeV <br />
mimrec energy range selection (for continuum sources): 0 - 5 MeV <br />

Transmission probability calculated for 33 km. <br />
 
Individual files: <br />

- Scaled_Ling_BG_1x.npz: the BG model (.npz) generated from C. Karwin's scaled 1x Ling BG simulation. 
- 511keV_imaging_response.npz: the 6deg imaging response required for RL imaging of Galactic positron annihilation.
- 1809keV_imaging_response.npz: the 6deg imaging response required for RL imaging of Galactic 26Al. 
- **Al26_10xFlux with no background**: Al26_10xFlux_Only.inc1.id1.extracted.tra.gz
- **Al26_10xFlux with Ling BG**: DC1_Al26_10xFlux_and_Ling.inc1.id1.extracted.tra.gz
- **511_10xFlux with no background**: GC_511_10xFlux_only.inc1.id1.extracted.tra.gz
- **511_10xFlux with Ling BG**: GC511_10xFlux_and_Ling.inc1.id1.extracted.tra.gz
- DC1_combined_10x.tra.gz: 4 point sources with 10x flux (Crab, Cyg X1, Cen A, Vela), 511 kev & Al26 lines, 1x Ling background
- Point_sources_10x.tra.gz: 4 point sources with 10x flux (Crab, Cyg X1, Cen A, Vela)
- Crab_only_10x.tra.gz: Crab with 10x flux
- CygX1_only_10x.tra.gz: Cyg X1 with 10x flux
- CenA_only_10x.tra.gz: Cen A with 10x flux
- Vela_only_10x.tra.gz: Vela with 10x flux
- Crab_BG_10x.tra.gz: Crab with 10x flux and 1x Ling background
- CenA_BG_10x.tra.gz: Cen A with 10x flux and 1x Ling background
- Point_sources_10x_BG.tra.gz: 4 point sources with 10x flux (Crab, Cyg X1, Cen A, Vela) and 1x Ling background
- Scaled_Ling_BG_1x.tra.gz: C. Karwin's scaled 1x Ling background simulation
- Scaled_Ling_BG_3x.npz: background model generated from C. Karwin's scaled 3x Ling BG simulation
- Continuum_Response.npz: 6deg response used for spectral analysis


