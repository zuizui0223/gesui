# Authentic 2019 HydroATLAS local and upstream land-use: minimal original source audit before new ecological confounder fits

2026-11-12. **POSTHOC after observing** the EEA H7 n454 spatial static-WWTP/EQR and original annual reported N/P/EQR results. This is a new original source technical availability check, not an ecological result.

Source HydroATLAS/BasinATLAS v1.0, publisher DOI https://doi.org/10.6084/m9.figshare.9890531 , HydroSHEDS https://www.hydrosheds.org/hydroatlas . The official original full `BasinATLAS_Data_v10_shp.zip` is approximately 3.98GB containing 12 Level shapefile databases. We need only the *original level7 dBASE attribute table*, not the entire 4GB polygon data.

**Target environmentally meaningful covariates, chosen before opening original source values**:
- `crp_pc_sse`/original document-backed cropland % in the local Level7 polygon; prefer original explicit cropland extent to generic land cover mosaic.
- `urb_pc_sse` original local urban extent and `glc_pc_s22` original artificial/built landcover percentage as an additional consistency crosscheck.
- `crp_pc_use` and `glc_pc_u22` original upstream cumulative watershed proportions if officially present.
- `dis_m3_pyr` original modelled mean river discharge at basin outlet, explicitly not a time-varying 2018/22 measured station river flow.
- Original `HYBAS_ID` must be unique one-to-one with the source-fixed 1,219 river-monitored Level7 basins and 232 EU macro study basins. Do not join by rounded latitude or alias when IDs missing.
- Original static 2019 source covariates fixed in both 2017–19 and 2020–22 rows. Do NOT model them as actual land use changes during the study.

## Technical first stage: do NOT download the whole 4GB file unless required

Use author/publisher Figshare API `https://api.figshare.com/v2/articles/9890531`, identify exactly named `BasinATLAS_Data_v10_shp.zip`, file ID, published source MD5 and size, native publisher download URL, response SHA. Query HTTP HEAD and two tiny, explicit Range requests (`bytes=0-1023`, `bytes=-65536`) to determine whether original publisher supports partial source reading. Close responses after bounded content. **If Range returns full HTTP200 or the original source is inaccessible, SOURCE STOP**, report route and never claim agricultural/urban controls fitted.

If HTTP206 works and content-length correct, download only the original Level07 DBF entry via HTTP Range seeking (e.g., fsspec HTTPFileSystem and zipfile), verify the original ZIP central-directory/member metadata and original record IDs before any covariate magnitude analysis. Record original uncompressed DBF SHA256, precise native file name, and list original cropland/artificial/river discharge columns. No conversion from unsupported source or invented defaults. Source-level missingness and numeric coding must be audited and original scale/units checked before entering any model.

No merging into original main/chocho, no eBird. Positive EQR/N from preceding PR #66 must remain unchanged; land-cover adjustment may weaken or reverse it and that must be reported.


**2026-11-12 source-schema correction, before any new land-use effect coefficients:** Official Level07 DBF ZIP-member download was successful (original `BasinATLAS_v10_shp/BasinATLAS_v10_lev07.dbf`, original uncompressed 80,137,382 bytes / compressed 15,674,749). The first original raw field check found `crp_pc_sse` and `crp_pc_use`, not the mistaken guessed `crp_pc_syr`/`crp_pc_uyr`. Independent published BasinATLAS catalog confirms `crp_pc_sse` cropland extent, `crp_pc_use` upstream cropland extent, and `urb_pc_sse` urban extent. This is a SOURCE naming correction, before inspecting any landcover/river biotic association. Original river ID, years, 1219 source sample and fitted prior outcome signs remain unchanged. Original pre-fix source failed closed, not silently substituted. No land cover effect claim until the corrected original-source QA passes.
