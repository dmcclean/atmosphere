{-# OPTIONS_GHC -Wall #-}

-- | This is a thin wrapper on top of the atmos package.
-- It provides the exaxt same functions but uses units
-- from the dimensional package.
module Atmosphere.Dimensional
       ( Atmos(..)
       , atmosphere
       , altitudeFromPressure
       ) where

import qualified Prelude ()
import Numeric.Units.Dimensional.Prelude
import qualified Atmosphere.Constants as A

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
                 { atmosTemperature        = A._TZERO * theta
                 , atmosPressure           = A._PZERO * delta
                 , atmosDensity            = density
                 , atmosSpeedOfSound       = A._AZERO * sqrt theta
                 , atmosViscosity          = viscosity
                 , atmosKinematicViscosity = viscosity/density
                 }
  where
    (sigma, delta, theta) = A.atmosphere alt
    t = (theta * A._TZERO) / (1 *~ kelvin)
    viscosity = A._BETAVISC*sqrt(t*t*t)/(t+(A._SUTH / (1 *~ kelvin)))
    density = A._RHOZERO * sigma

altitudeFromPressure :: (Floating a, Ord a) => Pressure a -> Length a
altitudeFromPressure = undefined

{-
altitudeFromPressure dimPressure = siAltitude *~ meter
  where
    siAltitude = A.siAltitudeFromPressure siPressure
    siPressure = dimPressure /~ pascal
-}
