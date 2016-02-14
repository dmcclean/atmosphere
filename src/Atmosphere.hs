{-# OPTIONS_GHC -Wall #-}
{-# LANGUAGE NoImplicitPrelude #-}

module Atmosphere
       ( Atmos(..)
       , siAtmosphere
       , usAtmosphere
       , atmosphere
       , siAltitudeFromPressure
       ) where

import qualified Atmosphere.Dimensional as D
import qualified Atmosphere.Constants as A
import Numeric.Units.Dimensional.Prelude
import Numeric.Units.Dimensional.NonSI (degreeRankine, foot, slug, poundForce)

data Atmos a = Atmos { atmosTemperature :: a
                     , atmosPressure :: a
                     , atmosDensity :: a
                     , atmosSpeedOfSound :: a
                     , atmosViscosity :: a
                     , atmosKinematicViscosity :: a
                     }
  deriving (Show)

{- |
   atmosphere in SI units

   Input: altitude in meters

   Output: (pressure, density, speed of sound, viscosity, kinematic viscosity)

   > pressure            - N/m^2
   > density             - kg/m^3
   > speed of sound      - m/s
   > viscosity           - N-s/m^2
   > kinematic viscosity - m^2/s
-}
siAtmosphere :: (Floating a) => a -> Atmos a
siAtmosphere alt_m =
  Atmos { atmosTemperature = temp /~ kelvin
        , atmosPressure = pressure /~ (newton / square meter)
        , atmosDensity = density /~ (kilo gram / cubic meter)
        , atmosSpeedOfSound = asound /~ (meter / second)
        , atmosViscosity = viscosity /~ (newton * second / meter^pos2)
        , atmosKinematicViscosity = kinematicViscosity /~ (meter^pos2 / second)
        }
  where
    D.Atmos temp pressure density asound viscosity kinematicViscosity = D.atmosphere (alt_m *~ meter)

{- |
   atmosphere in imperial units

   Input: altitude in ft

   Output: (pressure, density, speed of sound, viscosity, kinematic viscosity)

   > pressure            - lb/ft^2
   > density             - slugs/ft^3
   > speed of sound      - ft/s
   > viscosity           - slugs/(ft-s)
   > kinematic viscosity - ft^2/s
-}
usAtmosphere :: (Floating a) => a -> Atmos a
usAtmosphere alt_ft =
  Atmos { atmosTemperature = temp /~ degreeRankine
        , atmosPressure = pressure /~ (poundForce / square foot)
        , atmosDensity = density /~ (slug / cubic foot)
        , atmosSpeedOfSound = asound /~ (foot / second)
        , atmosViscosity = viscosity /~ (slug / (foot * second))
        , atmosKinematicViscosity = kinematicViscosity /~ (foot^pos2 / second)
        }
  where
    D.Atmos temp pressure density asound viscosity kinematicViscosity = D.atmosphere (alt_ft *~ foot)

{- |
   Compute altitude at which the standard atmosphere has a certain pressure.

   Input: Pressure, N/m^2

   Output: Altitude in meters
-}
siAltitudeFromPressure :: (Floating a) => a -> a
siAltitudeFromPressure = (/~ meter) . D.altitudeFromPressure . (*~ (newton / square meter))

{- |
   Compute temperature, density, and pressure in standard atmosphere.

   Correct to 86 km.  Only approximate thereafter.

   Input: alt geometric altitude, km.

   Output: (sigma, delta, theta)

   > sigma - density/sea-level standard density
   > delta - pressure/sea-level standard pressure
   > theta - temperature/sea-level std. temperature
-}
atmosphere :: (Floating a) => a -> (a,a,a)
atmosphere alt = (sigma /~ one, delta /~ one, theta /~ one)
  where
    (sigma, delta, theta) = A.atmosphere (alt *~ kilo meter)
