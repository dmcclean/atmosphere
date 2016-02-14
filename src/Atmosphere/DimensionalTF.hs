{-# OPTIONS_GHC -Wall #-}

-- | This is a thin wrapper on top of the atmos package.
-- It provides the exaxt same functions but uses units
-- from the dimensional-tf package.
module Atmosphere.DimensionalTF
       ( Atmos(..)
       , atmosphere
       , altitudeFromPressure
       ) where

import qualified Atmosphere as A

import qualified Prelude ()
import Numeric.Units.Dimensional.TF.Prelude

data Atmos a = Atmos { atmosTemperature :: ThermodynamicTemperature a
                     , atmosPressure :: Pressure a
                     , atmosDensity :: Density a
                     , atmosSpeedOfSound :: Velocity a
                     , atmosViscosity :: DynamicViscosity a
                     , atmosKinematicViscosity :: KinematicViscosity a
                     }
  deriving (Show)

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
