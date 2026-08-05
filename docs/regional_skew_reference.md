# USGS regional (generalized) skew coefficients — state reference

Research deliverable backing `hydrolib.regional_skew`. Every row below is either
a value confirmed from a cited USGS (or state DOT) report, or an explicit
record of why no single number could responsibly be assigned — nothing here
is estimated or interpolated.

Bulletin 17C generalized skew is not a nationwide lookup by state: it comes
from separate, periodically-revised state or multi-state regression studies,
with the Bulletin 17B nationwide map (&minus;0.302, SE 0.55) as the fallback
where no study exists. Only the **Constant** and **Defers to nationwide**
rows are implemented in `hydrolib.regional_skew.STATE_SKEW`; everything else
is documented here for reference but intentionally left out of the module's
lookup table rather than guessed.

This research pass could not fetch `pubs.usgs.gov` PDFs directly (network
policy in the research environment), so values come from indexed report
abstracts/figures surfaced via search, not full-text extraction. Treat every
row as a starting point for verification against the source, not a final
engineering value. The live, authoritative index is the USGS
[Flood Frequency Reports](https://www.usgs.gov/streamstats/science/flood-frequency-reports)
page.

An interactive, searchable/sortable version of this table (with expandable
citations and notes) is published at:
https://claude.ai/code/artifact/94f4f722-3b4c-4b5b-ac76-44993fc0ba09

## Summary

| Status | Count |
|---|---|
| Constant (in module) | 12 |
| Defers to nationwide (in module) | 4 |
| Value found, SE unconfirmed | 2 |
| Partial coverage only | 7 |
| Spatially varying | 17 |
| Unconfirmed | 10 |
| **Total** | **52** (50 states + DC + PR) |

## Constant — implemented in the module

| State | Value | SE | MSE | Source | Notes |
|---|---:|---:|---:|---|---|
| Arizona (AZ) | -0.090 | 0.283 | 0.080 | [SIR 2014-5211](https://pubs.usgs.gov/sir/2014/5211/) | Bayesian GLS, 448 gages; no basin characteristic improved on the constant model. |
| Connecticut (CT) | 0.370 | 0.374 | 0.140 | [SIR 2017-5037](https://pubs.usgs.gov/publication/sir20175037) | Shared New England constant with RI, MA. |
| Georgia (GA) | 0.048 | 0.303 | 0.092 | [SIR 2023-5006](https://pubs.usgs.gov/publication/sir20235006) | Shared constant with SC, NC. |
| Maine (ME) | 0.029 | 0.297 | 0.088 | [WRI 99-4008](https://pubs.usgs.gov/publication/wri994008) | SE given directly in source. |
| Massachusetts (MA) | 0.370 | 0.374 | 0.140 | [SIR 2017-5037](https://pubs.usgs.gov/publication/sir20175037) | Shared New England constant with CT, RI. |
| Missouri (MO) | -0.300 | 0.374 | 0.140 | [SIR 2014-5165](https://pubs.usgs.gov/sir/2014/5165/) | Bayesian WLS/GLS, 108 long-term gages. |
| North Carolina (NC) | 0.048 | 0.303 | 0.092 | [SIR 2023-5006](https://pubs.usgs.gov/publication/sir20235006) | Shared constant with GA, SC. |
| Rhode Island (RI) | 0.370 | 0.374 | 0.140 | [SIR 2017-5037](https://pubs.usgs.gov/publication/sir20175037) | Shared New England constant with CT, MA. |
| South Carolina (SC) | 0.048 | 0.303 | 0.092 | [SIR 2023-5006](https://pubs.usgs.gov/publication/sir20235006) | Shared constant with GA, NC. |
| Vermont (VT) | 0.440 | 0.279 | 0.078 | [SIR 2014-5078](https://pubs.usgs.gov/sir/2014/5078/) | Newer SIR 2025-5088 may supersede — verify. |
| Virginia (VA) | 0.500 | 0.574 | 0.330 | [SIR 2025-5110](https://pubs.usgs.gov/publication/sir20255110) | Shared Mid-Atlantic/South-Atlantic-Gulf constant with WV. |
| West Virginia (WV) | 0.500 | 0.574 | 0.330 | [SIR 2025-5110](https://pubs.usgs.gov/publication/sir20255110) | Shared constant with VA. |

## Defers to nationwide — implemented in the module

| State | Value | SE | Source | Notes |
|---|---:|---:|---|---|
| Montana (MT) | -0.302 | 0.640 | [WRI 03-4308](https://pubs.usgs.gov/wri/wri03-4308) | Nationwide value kept, but SE recalibrated for MT gages (vs. published 0.55). |
| North Dakota (ND) | -0.302 | 0.550 | [SIR 2015-5096](https://pubs.usgs.gov/publication/sir20155096) | B17B map found adequate for ND. |
| South Dakota (SD) | -0.302 | 0.550 | [WRI 98-4055](https://pubs.usgs.gov/wri/wri98-4055/) | B17B map found adequate for SD. |
| Washington (WA) | -0.302 | 0.550 | [SIR 2016-5118](https://pubs.usgs.gov/sir/2016/5118/) | No dedicated state study exists. |

## Value found, SE unconfirmed — not in the module

| State | Value | Source | Notes |
|---|---:|---|---|
| Arkansas (AR) | -0.170 | [SIR 2016-5081](https://pubs.usgs.gov/sir/2016/5081/) | Developed from a regional dataset spanning AR/LA/parts of MO,OK. |
| Hawaii (HI) | -0.140 | [WRI 84-4027](https://pubs.usgs.gov/wri/1984/4027/report.pdf) | Average of 68 gages; retained in a later (through WY2008) analysis. |

## Partial coverage only — not in the module

| State | Value | SE | MSE | Source | Coverage |
|---|---:|---:|---:|---|---|
| Delaware (DE) | 0.380 | 0.380 | 0.144 | [MDOT SHA MD-20-SP910B4G-2](https://www.roads.maryland.gov/OPR_Research/MD-20-SP910B4G-2_FixedRegionRegressionEquations-EasternCoastalPlain_Report.pdf) | Eastern Coastal Plain only; not a USGS report. |
| Illinois (IL) | 0.086 | 0.361 | 0.130 | [SIR 2019-5105](https://pubs.usgs.gov/sir/2019/5105/sir20195105.pdf) | "Parts of" the Great Lakes/Ohio River basins — exact IL extent unconfirmed. |
| Indiana (IN) | 0.086 | 0.361 | 0.130 | [SIR 2019-5105](https://pubs.usgs.gov/sir/2019/5105/sir20195105.pdf) | Same candidate/caveat as Illinois. |
| Maryland (MD) | 0.380 | 0.380 | 0.144 | [MDOT SHA MD-20-SP910B4G-2](https://www.roads.maryland.gov/OPR_Research/MD-20-SP910B4G-2_FixedRegionRegressionEquations-EasternCoastalPlain_Report.pdf) | Eastern Coastal Plain only; western MD instead falls in the VA/WV SIR 2025-5110 area. |
| New York (NY) | 0.320 | 0.332 | 0.110 | [SIR 2021-5015](https://pubs.usgs.gov/publication/sir20215015) | Eastern NY only; rest of state uses an older statewide isoline map. |
| Ohio (OH) | 0.086 | 0.361 | 0.130 | [SIR 2019-5105](https://pubs.usgs.gov/sir/2019/5105/sir20195105.pdf) | Same candidate/caveat as Illinois. |
| Pennsylvania (PA) | 0.320 | 0.332 | 0.110 | [SIR 2021-5015](https://pubs.usgs.gov/publication/sir20215015) | Eastern PA only. |

## Spatially varying — confirmed not to be a single constant

| State | Source | Model |
|---|---|---|
| Alabama (AL) | [SIR 2020-5032](https://www.sciencebase.gov/catalog/item/5c0fe01de4b0c53ecb2d1bb3) | 4 flood regions; uses B17B nationwide skew as base. |
| Alaska (AK) | [SIR 2016-5024](https://pubs.usgs.gov/wri/wri034188/) | 2 zones: eastern-central 0.54 (AVP 0.45); Gulf coast 0.18 (AVP 0.12). |
| California (CA) | [SIR 2010-5260](https://pubs.usgs.gov/publication/sir20105260) | Skew as a nonlinear function of mean basin elevation (-0.62 to 0.61). |
| Florida (FL) | [SIR 2011-5034](https://pubs.usgs.gov/sir/2011/5034/) | 4 hydrologic flood regions. |
| Idaho (ID) | [OFR 81-909](https://pubs.usgs.gov/of/1981/0909/report.pdf) | 3 maps selected by flood type (snowmelt/rainstorm/mixed). |
| Iowa (IA) | [SIR 2013-5086](https://pubs.usgs.gov/publication/wri004233) | Kriged isoline map, 239 gages. |
| Kansas (KS) | [WRI 2000-4079](https://pubs.usgs.gov/wri/2000/4079/report.pdf) | Isoline map, 253 gages. |
| Michigan (MI) | [WRI 83-4194](https://pubs.usgs.gov/publication/wri834194) | 3 zones: Upper Peninsula 0.12; SW Lower Peninsula 0.081; rest -0.17. |
| Minnesota (MN) | [WRI 97-4089](https://pubs.usgs.gov/publication/wri974089) | Regression trend-surface map, MSE 0.182. |
| Mississippi (MS) | [SIR 2018-5148](https://pubs.usgs.gov/sir/2018/5148/sir20185148.pdf) | 4 flood regions; some weight B17B with station skew. |
| New Hampshire (NH) | [SIR 2008-5206](https://pubs.usgs.gov/sir/2008/5206/) | Statewide skew map (SE 0.298), not a constant. |
| New Mexico (NM) | [TX/OK/e.NM technique report](https://www.usgs.gov/publications/technique-estimate-generalized-skew-coefficients-annual-peak-streamflow-natural) | GAM spatial surface (MSE 0.216); eastern NM only, western NM uncovered. |
| Oklahoma (OK) | [WRI 84-4358](https://www.usgs.gov/publications/technique-estimate-generalized-skew-coefficients-annual-peak-streamflow-natural) | Older statewide grid; also touches the TX/OK/eNM surface and AR's dataset. |
| Oregon (OR) | [SIR 2005-5116](https://pubs.usgs.gov/sir/2005/5116/) | 3 flood regions (western Oregon). |
| Tennessee (TN) | [SIR 2024-5130](https://pubs.usgs.gov/publication/sir20245130) | 4 hydrologic areas. |
| Texas (TX) | [TX/OK/e.NM technique report](https://www.usgs.gov/publications/technique-estimate-generalized-skew-coefficients-annual-peak-streamflow-natural) | GAM spatial surface (MSE 0.216), 341 gages. |
| Wyoming (WY) | [WRIR 03-4107](https://pubs.usgs.gov/wri/wri034107) | 6 hydrologic regions, 364 gages. |

## Unconfirmed — report exists but no usable number found

| State | Source |
|---|---|
| Colorado (CO) | [SIR 2009-5136](https://pubs.usgs.gov/publication/sir20095136) |
| District of Columbia (DC) | No dedicated study found; likely within neighboring MD/VA study areas. |
| Kentucky (KY) | Grouped with VA/WV/TN in Feaster and others, 2023 |
| Louisiana (LA) | Louisiana DOTD Water Resources Technical Report No. 60 (1998) |
| Nebraska (NE) | [WRIR 99-4032](https://www.usgs.gov/streamstats/science/flood-frequency-reports) |
| Nevada (NV) | [WSP 2433](https://www.usgs.gov/streamstats/science/flood-frequency-reports) (multi-state Southwestern US, 1997) |
| New Jersey (NJ) | [SIR 2009-5167](https://pubs.usgs.gov/publication/sir20095167) |
| Puerto Rico (PR) | [SIR 2021-5062](https://doi.org/10.3133/sir20215062) |
| Utah (UT) | [SIR 2007-5158](https://www.usgs.gov/streamstats/science/flood-frequency-reports) |
| Wisconsin (WI) | [SIR 2022-5118](https://www.usgs.gov/streamstats/science/flood-frequency-reports) |

## Methodology notes

- "MSE" in this document is whatever variance term the source study reported
  as its regional-skew error (mean square error or average variance of
  prediction) — treated as equivalent, following standard Bulletin 17C
  usage. SE is either reported directly by the source or derived as
  &radic;MSE.
- A state appearing in more than one section (e.g. New York: partial +
  spatially varying) means different parts of the state are covered by
  different studies with different characteristics.
- This is a starting point for future data-gathering, not a finished
  nationwide dataset. Extending coverage further requires primary-source
  access to the cited PDFs (blocked from this research environment) or a
  contributor pasting in verified values from those reports.
