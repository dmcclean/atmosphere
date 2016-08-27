{-# OPTIONS_GHC -Wall #-}

{-# LANGUAGE NegativeLiterals #-}

-- | This is a thin wrapper on top of the atmos package.
-- It provides the exaxt same functions but uses units
-- from the dimensional package.
module Atmosphere.Dimensional
       ( Atmos(..)
       , atmosphere
       , altitudeFromPressure
       , vaporPressureOfWaterAtSaturation
       , relativeHumidityFromTemperatureAndDewpoint
       ) where

import qualified Atmosphere as A

import qualified Prelude ()
import Numeric.Units.Dimensional.Prelude

data Atmos a = Atmos { atmosTemperature :: ThermodynamicTemperature a
                     , atmosPressure :: Pressure a
                     , atmosDensity :: Density a
                     , atmosSpeedOfSound :: Velocity a
                     , atmosViscosity :: DynamicViscosity a
                     , atmosKinematicViscosity :: KinematicViscosity a
                     }

atmosphere :: (Floating a, Ord a) => Length a -> Atmos a
atmosphere alt = Atmos
                 { atmosTemperature        = temp *~ kelvin
                 , atmosPressure           = pressure *~ pascal
                 , atmosDensity            = density *~ (kilo gram / meter ^ pos3)
                 , atmosSpeedOfSound       = asound *~ (meter / second)
                 , atmosViscosity          = viscosity *~ (newton * second / meter ^ pos2)
                 , atmosKinematicViscosity = kinematicViscosity *~ (meter ^ pos2 / second)
                 }
  where
    A.Atmos temp pressure density asound viscosity kinematicViscosity = A.siAtmosphere (alt /~ meter)

altitudeFromPressure :: (Floating a, Ord a) => Pressure a -> Length a
altitudeFromPressure dimPressure = siAltitude *~ meter
  where
    siAltitude = A.siAltitudeFromPressure siPressure
    siPressure = dimPressure /~ pascal

vaporPressureOfWaterAtSaturation :: (Floating a) => ThermodynamicTemperature a -> Pressure a
vaporPressureOfWaterAtSaturation t = (1 *~ pascal) * exp (a * t^pos2 + b * t + c + d * t^neg1)
  where
    a = 1.2378847e-5 *~ kelvin^neg2
    b = -1.9121316e-2 *~ kelvin^neg1
    c = 33.93711047 *~ one
    d = -6.3431645e3 *~ kelvin

enhancementFactor :: (Fractional a) => Pressure a -> ThermodynamicTemperature a -> Dimensionless a
enhancementFactor p t = alpha + beta * p + gamma * t^pos2
  where
    alpha = 1.00062 *~ one
    beta = 3.14e-8 *~ pascal^neg1
    gamma = 5.6e-7 *~ kelvin^neg2

relativeHumidityFromTemperatureAndDewpoint :: (Floating a) => Pressure a -> ThermodynamicTemperature a -> ThermodynamicTemperature a -> Dimensionless a
relativeHumidityFromTemperatureAndDewpoint p t td = (enhancementFactor p td * vaporPressureOfWaterAtSaturation td) / (enhancementFactor p t * vaporPressureOfWaterAtSaturation t)
