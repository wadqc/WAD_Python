"""Functions for converting plugin results to XML format"""

#Dictionary for translating limits tags in config XML to result XML
translate_limit_tags = {
    "acc_low": "grens_acceptabel_onder",
    "acc_high": "grens_acceptabel_boven",
    "crit_low": "grens_kritisch_onder",
    "crit_high": "grens_kritisch_boven",
    "criterium": "criterium",
}


def getResultsXML(result):
    """Convert a SingleResult object to XML."""

    default_xml = (
        '\t<type>{}</type>\n'
        '\t<omschrijving>{}</omschrijving>\n'
        '\t<waarde>{}</waarde>\n'
        '\t<grootheid>{}</grootheid>\n'
        '\t<eenheid>{}</eenheid>\n'
        '\t<niveau>{}</niveau>\n'
        '\t<volgnummer>{}</volgnummer>\n'
    )

    object_xml = (
        '\t<type>{}</type>\n'
        '\t<omschrijving>{}</omschrijving>\n'
        '\t<object_naam_pad>{}</object_naam_pad>\n'
        '\t<niveau>{}</niveau>\n'
        '\t<volgnummer>{}</volgnummer>\n'
    )

    if result.category == "object":
        out_xml = object_xml.format(
            result.category,
            result.description,
            result.value,
            result.level,
            result.id,
        )
    else:
        out_xml = default_xml.format(
            result.category,
            result.description,
            result.value,
            result.quantity,
            result.units,
            result.level,
            result.id,
        )

    return out_xml


def getXMLFromResult(result, limits):
    """Return the XML string for a SingleResult object and possible bounds"""

    bound_xml = ""
    for element in limits:
        #Check if a result has corresponding limits and gather the XML tags.
        if element.attrib["description"] == result.description:
            for child in element:
                #Translate config XML tag to result XML tag
                label = translate_limit_tags[child.tag]
                bound_xml += "\t<{0}>{1}</{0}>\n".format(label, child.text)

    xml_str = "<results>\n"
    xml_str += getResultsXML(result)
    xml_str += bound_xml
    xml_str += "</results>"
    return xml_str