# =============================================================================
# app.py — Streamlit UI for the Evolutionary DNA Matcher
#
# Handles multi-species comparison, heatmaps, and neighbor-joining summary.
# Run with: streamlit run app.py
# =============================================================================

import matplotlib.pyplot as plt
import streamlit as st

from dna_matcher.core import (
    search_species_options,
    fetch_common_marker_sequences,
    compare_species_matrix,
    build_neighbor_joining_tree,
)

st.set_page_config(page_title="DNA Matcher", page_icon="🧬", layout="wide")

st.title("🧬 Evolutionary DNA Matcher")
st.markdown(
    "Compare 3–5 species on the same DNA marker, view pairwise similarity, "
    "and inspect the best global and local matches."
)
st.markdown("---")


def get_species_flow(index: int):
    common = st.text_input(
        f"Species {index} common name",
        placeholder="e.g. goat, bison, chimaera",
        key=f"common_{index}",
    )
    org_type = st.selectbox(
        f"Organism type for species {index}",
        ["default", "animal", "plant", "fungus"],
        key=f"type_{index}",
        help="Helps select the right genetic marker. Use 'default' if unsure.",
    )

    sci_name = None
    if common:
        with st.spinner(f"Searching for '{common}'..."):
            options = search_species_options(common)

        if not options:
            st.error(f"No species found for '{common}'. Try a different name.")
            # If taxonomy debug info was recorded, show which queries were tried.
            dbg = st.session_state.get("tax_debug", {}).get(common)
            if dbg:
                st.info("Queries tried:")
                for q in dbg.get("queries_tried", []):
                    st.code(q)
                last = dbg.get("last_idlist") or []
                st.write(f"Last returned ID list: {last}")
        elif len(options) == 1:
            sci_name = options[0][1]
            st.info(f"Found: **{options[0][0]}** (*{sci_name}*)")
        else:
            labels = [f"{c}  —  {s}" for c, s in options]
            sci_names = [s for _, s in options]
            idx = st.selectbox(
                f"Multiple matches for '{common}' — pick one:",
                range(len(labels)),
                format_func=lambda i: labels[i],
                key=f"select_{index}",
            )
            sci_name = sci_names[idx]
            st.info(f"Selected: **{options[idx][0]}** (*{sci_name}*)")

    return sci_name, org_type


num_species = st.number_input(
    "How many species do you want to compare?",
    min_value=2,
    max_value=5,
    value=2,
    step=1,
    help="Choose the number of species to compare before entering their details.",
)

st.markdown("### Enter species names")

species_inputs = []
with st.expander(f"Enter details for {num_species} species", expanded=True):
    for idx in range(1, num_species + 1):
        species_inputs.append(get_species_flow(idx))

st.markdown("---")

with st.expander("⚙️ Advanced Settings", expanded=False):
    max_len = st.slider(
        "Sequence length to compare (base pairs)",
        min_value=1000,
        max_value=10000,
        value=1000,
        step=100,
        help="Higher = more accurate but slower. COI markers are ~650bp.",
    )

if st.button("🧬 Compare species", use_container_width=True):
    species_list = [sci for sci, _ in species_inputs if sci]
    types = [typ for sci, typ in species_inputs if sci]

    if len(species_list) < 2:
        st.warning("Please specify at least 2 species with valid scientific names.")
    else:
        progress = st.progress(0)
        status = st.empty()

        status.text("Fetching a common DNA marker for all species...")
        marker_result = fetch_common_marker_sequences(tuple(species_list), tuple(types))
        progress.progress(30)

        if not marker_result:
            st.error(
                "Could not find a shared marker for all selected species. "
                "Try removing one species or changing organism types."
            )
        else:
            marker = marker_result["marker"]
            sequences = marker_result["sequences"]
            status.text("Comparing species pairwise...")
            matrix_data = compare_species_matrix(sequences, marker, max_len=max_len)
            progress.progress(80)

            tree_text = build_neighbor_joining_tree(matrix_data)
            progress.progress(100)
            status.text("Done ✅")
            st.balloons()

            st.markdown("#### Genetic marker used")
            st.info(f"`{marker}` matched for all selected species")
            st.markdown("---")

            st.markdown("#### Pairwise similarity heatmap")
            species_names = matrix_data["species"]
            scores = [
                [matrix_data["matrix"][a][b] for b in species_names]
                for a in species_names
            ]

            fig, ax = plt.subplots(figsize=(8, 6))
            im = ax.imshow(scores, cmap="viridis", vmin=0, vmax=100)
            ax.set_xticks(range(len(species_names)))
            ax.set_yticks(range(len(species_names)))
            ax.set_xticklabels(species_names, rotation=45, ha="right")
            ax.set_yticklabels(species_names)
            for i, row in enumerate(scores):
                for j, value in enumerate(row):
                    text_color = "white" if value < 50 else "black"
                    ax.text(j, i, f"{value:.1f}", ha="center", va="center", color=text_color)
            fig.colorbar(im, ax=ax, label="Global similarity (%)")
            st.pyplot(fig)

            st.markdown("---")
            st.markdown("#### Best pairwise matches")
            best_global = matrix_data["best_global"]
            best_local = matrix_data["best_local"]
            g1, g2 = st.columns(2)
            g1.metric(
                "Best global alignment",
                f"{best_global['pair'][0]} vs {best_global['pair'][1]}",
                f"{best_global['score']:.2f}%",
            )
            g2.metric(
                "Best local alignment",
                f"{best_local['pair'][0]} vs {best_local['pair'][1]}",
                f"{best_local['score']:.2f}%",
            )

            st.markdown("---")
            st.markdown("#### Neighbor-joining tree approximation")
            st.code(tree_text, language="text")

            st.markdown("---")
            st.markdown("#### Pair detail explorer")
            pair_labels = [
                f"{a} vs {b}"
                for a in species_names
                for b in species_names
                if a != b
            ]
            selected_pair = st.selectbox("Select a pair to inspect", pair_labels)
            left, right = selected_pair.split(" vs ")
            detail = matrix_data["details"][(left, right)]

            with st.expander("Global alignment match segment"):
                st.code(detail["hirschberg_match"] or "(no match)", language="text")
            with st.expander("Local alignment match segment"):
                st.code(detail["sw_match"] or "(no match)", language="text")
            st.write(
                f"Compare length: {detail['compare_len']} bp, global: {detail['hirschberg_pct']:.2f}%, "
                f"local: {detail['sw_pct']:.2f}% (score {detail['sw_score']})"
            )
