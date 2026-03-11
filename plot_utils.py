"""
AlphaGenome Explorer — Plot building utilities.

Functions for constructing plot components from AlphaGenome prediction results.
"""

from alphagenome.visualization import plot_components

from config import OUTPUT_ATTR_MAP


def build_interval_components(result, selected_outputs, transcripts):
    """Build plot components for interval prediction results.

    Returns:
        list of plot components (TranscriptAnnotation + signal tracks).
    """
    components = [plot_components.TranscriptAnnotation(transcripts)]

    for key in selected_outputs:
        attr = OUTPUT_ATTR_MAP.get(key)
        tdata = getattr(result, attr, None) if attr else None
        if tdata is None or tdata.values.shape[-1] == 0:
            continue

        if key == "CONTACT_MAPS":
            components.append(plot_components.ContactMaps(
                tdata=tdata, ylabel_template=f"{key}: {{name}}",
            ))
        elif key == "SPLICE_JUNCTIONS":
            components.append(plot_components.Sashimi(
                junction_track=tdata, ylabel_template=f"{key}: {{name}} ({{strand}})",
            ))
        else:
            components.append(plot_components.Tracks(
                tdata=tdata,
                ylabel_template=f"{key}: {{name}} ({{strand}})\n{{name}}",
                filled=True,
            ))

    return components


def build_variant_components(result, selected_outputs, transcripts):
    """Build plot components for variant prediction results.

    Shows overlaid REF/ALT tracks + ALT-REF difference tracks.

    Returns:
        list of plot components.
    """
    components = [plot_components.TranscriptAnnotation(transcripts)]

    for key in selected_outputs:
        attr = OUTPUT_ATTR_MAP.get(key)
        ref_tdata = getattr(result.reference, attr, None) if attr else None
        alt_tdata = getattr(result.alternate, attr, None) if attr else None
        if ref_tdata is None or alt_tdata is None:
            continue
        if ref_tdata.values.shape[-1] == 0:
            continue

        if key == "CONTACT_MAPS":
            diff = alt_tdata - ref_tdata
            components.append(plot_components.ContactMapsDiff(
                tdata=diff, ylabel_template=f"{key}: {{name}}",
            ))
        elif key == "SPLICE_JUNCTIONS":
            components.append(plot_components.Sashimi(
                junction_track=ref_tdata,
                ylabel_template=f"REF {key}: {{name}} ({{strand}})",
            ))
            components.append(plot_components.Sashimi(
                junction_track=alt_tdata,
                ylabel_template=f"ALT {key}: {{name}} ({{strand}})",
            ))
        else:
            # Overlaid REF (gray) vs ALT (red)
            components.append(plot_components.OverlaidTracks(
                tdata={"REF": ref_tdata, "ALT": alt_tdata},
                colors={"REF": "dimgrey", "ALT": "red"},
                ylabel_template=f"{key}: {{name}} ({{strand}})\n{{name}}",
                shared_y_scale=True,
            ))
            # ALT - REF difference below
            diff_tdata = alt_tdata - ref_tdata
            components.append(plot_components.Tracks(
                tdata=diff_tdata,
                ylabel_template=f"{key} (ALT-REF): {{name}} ({{strand}})\n{{name}}",
                filled=True,
            ))

    return components


def count_interval_tracks(data):
    """Count number of track rows for dynamic plot height."""
    if data is None:
        return 1
    result, _, selected_outputs = data
    count = 1  # transcript annotation
    for key in selected_outputs:
        attr = OUTPUT_ATTR_MAP.get(key)
        tdata = getattr(result, attr, None) if attr else None
        if tdata is not None and tdata.values.shape[-1] > 0:
            count += tdata.values.shape[-1]
    return count


def count_variant_tracks(data):
    """Count track rows for dynamic variant plot height."""
    if data is None:
        return 1
    result, _, selected_outputs, _ = data
    count = 1  # transcript annotation
    for key in selected_outputs:
        attr = OUTPUT_ATTR_MAP.get(key)
        ref_tdata = getattr(result.reference, attr, None) if attr else None
        if ref_tdata is not None and ref_tdata.values.shape[-1] > 0:
            n = ref_tdata.values.shape[-1]
            if key == "SPLICE_JUNCTIONS":
                count += n * 2  # REF + ALT sashimi
            else:
                count += n * 2  # overlaid + diff
    return count
