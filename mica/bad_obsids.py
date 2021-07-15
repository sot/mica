# Licensed under a 3-clause BSD sytle license - see LICENSE

# Venus
venus = [
    583,
    2411,
    2414,
    6395,
    7306,
    7307,
    7308,
    7309,
    7310,
    7311,
    7312,
    7313,
    7314,
    7315,
    7316,
    7406,
    9741,
    9742,
    9743,
    9744,
    9745,
    9746,
    9747,
    9748,
    9749,
    9752,
    9753,
    7316,
    15292,
    16499,
    16500,
    16501,
    16502,
    16503,
    16504,
    16505,
    16506,
    17439,
    18690,
    18691,
    18692,
    18693,
    18694,
    18695,
    18696,
    53250]

multi_obi = [
    943,
    897,
    60879,
    2042,
    60881,
    800,
    1900,
    2365,
    906,
    2010,
    3182,
    2547,
    380,
    3057,
    2077,
    60880,
    2783,
    1578,
    1561,
    4175,
    3764,
]

# Observations with nonload commands changing catalog or other oddities
weird = [50707]


def bad_obsids():
    """
    Returns a list of the obsids that are reasonable to exclude from most trending applications.
    The list includes observations of Venus and observations with multiple obis. The observations
    of Venus are problematic for star trending.  The observations with multiple obis may have more
    than one attitude for the same obsid.

    The lists of each (venus, multi_obi, weird) are maintained directly in the module and may be
    used directly if needed.

    :returns: list
    """
    return venus + multi_obi + weird
