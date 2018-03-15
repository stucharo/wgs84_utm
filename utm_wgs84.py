from math import floor, radians, degrees, sqrt, cos, sin, tan, atan, atanh, tanh, atan2, asinh, sinh, cosh
import numpy as np


class Ellipsoid():

    def __init__(self):
        self.a = 6378137
        self.b = 6356752.314245
        self.f = 1 / 298.257223563


class Datum():

    def __init__(self):
        self.ellipsoid = Ellipsoid()


class LatLon():

    def __init__(self, lat, lon, datum=None, scale=None, convergence=None):
        if not (type(lat) is np.float32 and type(lon) is np.float32):
            raise ValueError('Invalid point')
        self.lat = lat
        self.lon = lon
        self.scale = scale
        self.convergence = convergence

        if datum is None:
            self.datum = Datum()

    def __str__(self):
        return f'{self.lat}°E, {self.lon}°N'

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

        φ = radians(self.lat)         # latitude ± from equator
        λ = radians(self.lon) - λ0      # longitude ± from central meridian

        a = self.datum.ellipsoid.a
        f = self.datum.ellipsoid.f
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

        α = [None,  # note α is one-based array (6th order Krüger expressions)
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
        x += false_easting             # make x relative to false easting
        if y < 0:
            y += false_northing          # make y in southern hemisphere relative to false northing

        # round to reasonable precision
        x = round(x, 6)  # nm precision
        y = round(y, 6)  # nm precision
        convergence = round(degrees(γ), 9)
        scale = round(k, 12)

        if self.lon >= 0:
            h = 'N'
        else:
            h = 'S'

        return UTM(zone, h, x, y, self.datum, convergence, scale)


class UTM():

    def __init__(self, zone, hemisphere, easting, northing, datum, convergence=None, scale=None):
        self.zone = zone
        self.hemisphere = hemisphere
        self.easting = easting
        self.northing = northing
        self.datum = datum
        self.convergence = convergence
        self.scale = scale
        self.datum = Datum()

    def __str__(self):
        return f'{self.zone} {self.hemisphere.upper()} {round(self.easting, 2)} {round(self.northing, 2)}'

    def to_latlon(self):
        z = self.zone
        h = self.hemisphere.upper()
        x = self.easting
        y = self.northing

        false_easting = 500e3
        false_northing = 10000e3

        a = self.datum.ellipsoid.a
        f = self.datum.ellipsoid.f
        # WGS 84:  a = 6378137, b = 6356752.314245, f = 1/298.257223563

        k0 = 0.9996  # UTM scale on the central meridian

        x = x - false_easting               # make x ± relative to central meridian
        if h == 'S':                         # make y ± relative to equator
            y -= false_northing

        # ---- from Karney 2011 Eq 15-22, 36:

        e = sqrt(f * (2 - f))  # eccentricity
        n = f / (2 - f)        # 3rd flattening
        n2 = n * n
        n3 = n * n2
        n4 = n * n3
        n5 = n * n4
        n6 = n * n5

        # 2πA is the circumference of a meridian
        A = a / (1 + n) * (1 + 1 / 4 * n2 + 1 / 64 * n4 + 1 / 256 * n6)

        η = x / (k0 * A)
        ξ = y / (k0 * A)

        β = [None,  # note β is one-based array (6th order Krüger expressions)
             1 / 2 * n - 2 / 3 * n2 + 37 / 96 * n3 - 1 / 360 * \
             n4 - 81 / 512 * n5 + 96199 / 604800 * n6,
             1 / 48 * n2 + 1 / 15 * n3 - 437 / 1440 * n4 + \
             46 / 105 * n5 - 1118711 / 3870720 * n6,
             17 / 480 * n3 - 37 / 840 * n4 - 209 / 4480 * n5 + 5569 / 90720 * n6,
             4397 / 161280 * n4 - 11 / 504 * n5 - 830251 / 7257600 * n6,
             4583 / 161280 * n5 - 108847 / 3991680 * n6,
             20648693 / 638668800 * n6]

        ξʹ = ξ
        ηʹ = η
        p = 1
        q = 0
        for j in range(1, 6):
            ξʹ -= β[j] * sin(2 * j * ξ) * cosh(2 * j * η)
            ηʹ -= β[j] * cos(2 * j * ξ) * sinh(2 * j * η)
            # ---- convergence: Karney 2011 Eq 26, 27
            p -= 2 * j * β[j] * cos(2 * j * ξ) * cosh(2 * j * η)
            q += 2 * j * β[j] * sin(2 * j * ξ) * sinh(2 * j * η)

        sinhηʹ = sinh(ηʹ)
        sinξʹ = sin(ξʹ)
        cosξʹ = cos(ξʹ)

        τʹ = sinξʹ / sqrt(sinhηʹ * sinhηʹ + cosξʹ * cosξʹ)

        τi = τʹ
        δτi = 1

        # using IEEE 754 δτi -> 0 after 2-3 iterations
        # note relatively large convergence test as δτi toggles on ±1.12e-16 for eg 31 N 400000 5000000
        while (abs(δτi) > 1e-12):
            σi = sinh(e * atanh(e * τi / sqrt(1 + τi * τi)))
            τiʹ = τi * sqrt(1 + σi * σi) - σi * sqrt(1 + τi * τi)
            δτi = (τʹ - τiʹ) / sqrt(1 + τiʹ * τiʹ) * (1 + (1 - e * e)
                                                      * τi * τi) / ((1 - e * e) * sqrt(1 + τi * τi))
            τi += δτi

        τ = τi

        φ = atan(τ)

        λ = atan2(sinhηʹ, cosξʹ)

        γʹ = atan(tan(ξʹ) * tanh(ηʹ))
        γʺ = atan2(q, p)

        γ = γʹ + γʺ

        # ---- scale: Karney 2011 Eq 28

        sinφ = sin(φ)
        kʹ = sqrt(1 - e * e * sinφ * sinφ) * sqrt(1 + τ * τ) * \
            sqrt(sinhηʹ * sinhηʹ + cosξʹ * cosξʹ)
        kʺ = A / a / sqrt(p * p + q * q)

        k = k0 * kʹ * kʺ

        # ------------

        λ0 = radians((z - 1) * 6 - 180 + 3)  # longitude of central meridian
        λ += λ0  # move λ from zonal to global coordinates

        # round to reasonable precision
        lat = round(degrees(φ), 11)  # nm precision (1nm = 10^-11°)
        # (strictly lat rounding should be φ⋅cosφ!)
        lon = round(degrees(λ), 11)
        convergence = round(degrees(γ), 9)
        scale = round(k, 12)

        return LatLon(lat, lon, self.datum, convergence, scale)
