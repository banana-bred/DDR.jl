using DDR
using Plots

root = "data/BEND"   # folder that contains eigenph.all.geom*

# -- Load all records available recursively in root directory
records = DDR.load_eigenph_records(root)

# -- sort by Q
sort!(records, by = r -> (ismissing(r.meta.Q) ? Inf : r.meta.Q))

# -- preprocess; fix jumps in major resonances, and try to center
#    the resonances so theyre all more or less in the same branch
DDR.unwrap!(records; Emin=0.01, period=pi)
DDR.center_records_at!(records; Eref=0.01, period=pi, targtype=:minabs)

# -- define a selector for choosing which resonance(s) we want to look at.
#    This can be a function, a string, a vector of strings, or a regex (as below).
# selector = ["singlet A1", "triplet B2"]
# selector = "singlet A1"
# selector = r"(singlet A1|triplet B2)"
selector = r"(singlet A1)"

# -- get a list of indexed resonances matching our selector, with some
#    other variables used to filter resonances. Here, we only get the largest ones
idxres = DDR.find_resonances_across_records(records;
    selector=selector,
    Emin=0.05,
    smoothwin=5,
    min_height=40,
)

# -- plot δ as a function of E
plt_phase = DDR.plot_eigenph(records
    ; selector=selector
    , idxres=idxres
    , show=true
    , xlims=(0.05, :auto)
)
display(plt_phase)

# -- plot dδ/dE as a function of E
plt_deriv = DDR.plot_eigenph(records
    ; selector=selector
    , ytype=:deriv
    , smoothwin=5
    , idxres=idxres
    , show=true
    , xlims=(0.05, :auto)
    , ylims=(0, 200)
)
display(plt_deriv)

# -- generate a simple table of the resonances that we've found so far, and show it
longtable = resonance_table_long(records, idxres)
show_table(longtable)

################
### TRACKING ###
################

# -- Here we try and single out a resonance and track it across records (Q, geom, etc.)

# -- track a 1A1 resonance
selector = "singlet A1"
# get the 1A1 resonances
idxres_BEND_1A1 = DDR.find_resonances_across_records(records;
    selector=selector,
    Emin=0.05,
    smoothwin=5,
    min_height=40,
    min_prominence=0.5,
)

# -- get the 1A1 resonance in the window centered around 0.3 eV and follow it
tracked_BEND_1A1 =
  track_resonance(
      records
    , idxres_BEND_1A1
    ; selector = selector
    , center = 0.3
    , window = 0.4
    , follow_window = 0.1
  )

# -- track a 3B2 resonance
selector="triplet B2"
idxres_BEND_3B2 = DDR.find_resonances_across_records(records;
    selector=selector,
    Emin=0.05,
    smoothwin=5,
    min_height=40,
    min_prominence=0.5,
)

# -- get the 3B2 resonance in the given window, as above
tracked_BEND_3B2 =
  track_resonance(
      records
    , idxres_BEND_3B2
    ; selector = selector
    , center = 0.3
    , window = 0.4
    , follow_window = 0.1
  )

# -- this array of tuples defines our column pairs
specs = [
    (label = "BEND_1A1", tracked=tracked_BEND_1A1)
  , (label = "BEND_3B2", tracked=tracked_BEND_3B2)
]

# -- make a wide table and show it
widetable = resonance_table_wide(records, specs)
show_table(widetable)

# -- plot our tracked resonances
plt = plot_resonance(records, tracked_BEND_1A1, label="BEND 1A₁", connect=true, xflip=true, show=true)
plt = plot_resonance!(plt, records, tracked_BEND_3B2, label="BEND 3B₂", connect=true, xflip=true, show=true)
display(plt)

# -- write table to file
filename = DDR.write_table("important_resonances.tab", widetable; colsep="\t")
println("Wrote resonance table to $filename")
