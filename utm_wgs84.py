from math import floor, radians, sqrt, cos, sin, tan, atan, atanh, atan2, asinh, sinh, cosh


class LatLon():

    def __init__(self, lat, lon):
        if not (type(lat, float) and type(lon, float)):
            raise ValueError('Invalid point')
        self.lat = lat
        self.lon = lon

        self.datum.a = 6378137
        self.datum.b = 6356752.314245
        self.datum.f = 1 / 298.257223563

    def to_utm(self):
        if not -80 <= self.lat <= 84:
            raise ValueError('Outside UTM limits')

        false_easting = 500e3
        false_northing = 10000e3

        zone = floor((self.lon + 180) / 6) + 1  # longitudinal zone
        λ0 = radians((zone - 1) * 6 - 180 + 3)  # longitude of central meridian

        # handle Norway / Svalbard exceptions
        # grid zones are 8° tall; 0°N is offset 10 into latitude bands array
        mgrsLatBands = 'CDEFGHJKLMNPQRSTUVWXX'  # X is repeated for 80-84°N
        latBand = mgrsLatBands[floor(self.lat / 8 + 10)]
        # adjust zone & central meridian for Norway
        if (zone == 31 and latBand == 'V' and self.lon >= 3):
            zone += 1
            λ0 += radians(6)
        # adjust zone & central meridian for Svalbard
        if (zone == 32 and latBand == 'X' and self.lon < 9):
            zone -= 1
            λ0 -= radians(6)
        if (zone == 32 and latBand == 'X' and self.lon >= 9):
            zone += 1
            λ0 += radians(6)
        if (zone == 34 and latBand == 'X' and self.lon < 21):
            zone -= 1
            λ0 -= radians(6)
        if (zone == 34 and latBand == 'X' and self.lon >= 21):
            zone += 1
            λ0 += radians(6)
        if (zone == 36 and latBand == 'X' and self.lon < 33):
            zone -= 1
            λ0 -= radians(6)
        if (zone == 36 and latBand == 'X' and self.lon >= 33):
            zone += 1
            λ0 += radians(6)

        psi = radians(self.lat)         # latitude ± from equator
        λ = radians(self.lon) - λ0      # longitude ± from central meridian

        a = self.datum.a
        f = self.datum.f
        # WGS 84: a = 6378137, b = 6356752.314245, f = 1/298.257223563;

        k0 = 0.9996  # UTM scale on the central meridian

        #---- easting, northing: Karney 2011 Eq 7-14, 29, 35:

        e = sqrt(f * (2 - f))  # eccentricity
        n = f / (2 - f)  # 3rd flattening
        n2 = n * n
        n3 = n * n2
        n4 = n * n3
        n5 = n * n4
        n6 = n * n5

        cosλ = cos(λ)
        sinλ = sin(λ)
        tanλ = tan(λ)

        # τ ≡ tanφ, τʹ ≡ tanφʹ; prime (ʹ) indicates angles on the conformal sphere
        τ = tan(φ)
        σ = sinh(e * atanh(e * τ / sqrt(1 + τ * τ)))

        τʹ = τ * sqrt(1 + σ * σ) - σ * sqrt(1 + τ * τ)

        ξʹ = atan2(τʹ, cosλ)
        ηʹ = asinh(sinλ / sqrt(τʹ * τʹ + cosλ * cosλ))

        # 2πA is the circumference of a meridian
        A = a / (1 + n) * (1 + 1 / 4 * n2 + 1 / 64 * n4 + 1 / 256 * n6)

        α = [null,  # note α is one-based array (6th order Krüger expressions)
             1 / 2 * n - 2 / 3 * n2 + 5 / 16 * n3 + 41 / \
             180 * n4 - 127 / 288 * n5 + 7891 / 37800 * n6,
             13 / 48 * n2 - 3 / 5 * n3 + 557 / 1440 * n4 + \
             281 / 630 * n5 - 1983433 / 1935360 * n6,
             61 / 240 * n3 - 103 / 140 * n4 + 15061 / 26880 * n5 + 167603 / 181440 * n6,
             49561 / 161280 * n4 - 179 / 168 * n5 + 6601661 / 7257600 * n6,
             34729 / 80640 * n5 - 3418889 / 1995840 * n6,
             212378941 / 319334400 * n6]

        ξ = ξʹ
        η = ηʹ
        pʹ = 1
        qʹ = 0
        for j in range(1, 6):
            ξ += α[j] * sin(2 * j * ξʹ) * cosh(2 * j * ηʹ)
            η += α[j] * cos(2 * j * ξʹ) * sinh(2 * j * ηʹ)
            # ---- convergence: Karney 2011 Eq 23, 24
            pʹ += 2 * j * α[j] * cos(2 * j * ξʹ) * cosh(2 * j * ηʹ)
            qʹ += 2 * j * α[j] * sin(2 * j * ξʹ) * sinh(2 * j * ηʹ)

        x = k0 * A * η
        y = k0 * A * ξ

        γʹ = atan(τʹ / sqrt(1 + τʹ * τʹ) * tanλ)
        γʺ = atan2(qʹ, pʹ)

        γ = γʹ + γʺ

        # ---- scale: Karney 2011 Eq 25

        sinφ = sin(φ)
        kʹ = sqrt(1 - e * e * sinφ * sinφ) * \
            sqrt(1 + τ * τ) / sqrt(τʹ * τʹ + cosλ * cosλ)
        kʺ = A / a * sqrt(pʹ * pʹ + qʹ * qʹ)

        k = k0 * kʹ * kʺ

        # ------------

        # shift x/y to false origins
        x += falseEasting             # make x relative to false easting
        if y < 0:
            y += falseNorthing          # make y in southern hemisphere relative to false northing

        # round to reasonable precision
        x = round(x, 6)  # nm precision
        y = round(y, 6)  # nm precision
        convergence = round(degrees(γ), 9)
        scale = round(k, 12)

        if self.lat >= 0:
            h = 'N'
        else:
            h = 'S'

        return Utm(zone, h, x, y, this.datum, convergence, scale)


class UTM():

    def __init__(self, zone, easting, northing, datum, convergence=None, scale=None):
        self.zone = zone
        self.easting = easting
        self.northing = northing
        self.datum = datum
        self.convergence = convergence
        self.scale = scale
