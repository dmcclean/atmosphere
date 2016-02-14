{-# OPTIONS_GHC -Wall #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE NegativeLiterals #-}
{-# LANGUAGE NoImplicitPrelude #-}
{-# LANGUAGE PatternGuards #-}
{-# LANGUAGE TypeOperators #-}

module Atmosphere.Constants
       ( _TZERO
       , _PZERO
       , _RHOZERO
       , _AZERO
       , _BETAVISC
       , _SUTH
       , _REARTH
       , _GMR
       , Layer(..)
       , layerTable
       , atmosphere
       ) where

import Data.Maybe (fromMaybe)
import Numeric.Units.Dimensional.Prelude


-- | sea-level temperature
_TZERO :: (Floating a) => ThermodynamicTemperature a
_TZERO   = 288.15 *~ kelvin

-- | sea-level pressure
_PZERO :: (Floating a) => Pressure a
_PZERO   = 101325.0 *~ (newton / square meter)

-- | sea-level density
_RHOZERO :: (Floating a) => Density a
_RHOZERO = 1.225 *~ (kilo gram / cubic meter)

-- | speed of sound at sea-level
_AZERO :: (Floating a) => Velocity a
_AZERO   = 340.294 *~ (meter / second)

-- | viscosity constant
_BETAVISC :: (Floating a) => DynamicViscosity a
_BETAVISC = 1.458E-6 *~ (newton * second / meter^pos2)

-- | Sutherland's constant, kelvins
_SUTH :: (Floating a) => ThermodynamicTemperature a
_SUTH    = 110.4 *~ kelvin

-- | radius of the Earth, km
_REARTH :: (Floating a) => Length a
_REARTH = 6369.0 *~ kilo meter

-- | gas constant
_GMR :: (Floating a) => TemperatureGradient a
_GMR = 34.163195 *~ (kelvin / kilo meter)

type TemperatureGradient a = Quantity (DThermodynamicTemperature / DLength) a

data Layer a = Layer
             { baseHeight :: Length a
             , baseTemperature :: ThermodynamicTemperature a
             , baseSpecificPressure :: Dimensionless a
             , temperatureGradient' :: Maybe (TemperatureGradient a)
             , heightLimits :: (Maybe (Length a), Maybe (Length a))
             }
  deriving (Show)

layer :: (Floating a) => a -> a -> a -> Maybe a -> (Maybe a, Maybe a)-> Layer a
layer h t p g (bot,top) = Layer (h *~ kilo meter) (t *~ kelvin) (p *~ one) (g *~~ (kelvin / kilo meter)) (bot *~~ kilo meter, top *~~ kilo meter)

temperatureGradient :: (Floating a) => Layer a -> TemperatureGradient a
temperatureGradient = fromMaybe _0 . temperatureGradient'

layerTable :: (Floating a) => [Layer a]
layerTable = [ layer  0.0   (_TZERO/~kelvin)  1.0          (Just -6.5) (Nothing,     Just 11.0)
             , layer 11.0   216.65            2.2336110E-1 Nothing     (Just 11.0,   Just 20.0)
             , layer 20.0   216.65            5.4032950E-2 (Just 1.0)  (Just 20.0,   Just 32.0)
             , layer 32.0   228.65            8.5666784E-3 (Just 2.8)  (Just 32.0,   Just 47.0)
             , layer 47.0   270.65            1.0945601E-3 Nothing     (Just 47.0,   Just 51.0)
             , layer 51.0   270.65            6.6063531E-4 (Just -2.8) (Just 51.0,   Just 71.0)
             , layer 71.0   214.65            3.9046834E-5 (Just -2.0) (Just 71.0,   Just 84.852)
             , layer 84.852 186.946           3.68501E-6   Nothing     (Just 84.852, Nothing)
             ]

layerTemperature :: (Floating a) => Layer a -> Length a -> ThermodynamicTemperature a
layerTemperature l dh = baseTemperature l + dh * temperatureGradient l

layerDelta :: (Floating a) => Layer a -> Length a -> Dimensionless a
layerDelta l h | Just tgrad <- temperatureGradient' l = pbase*(tbase/tlocal)**(_GMR/tgrad)
               | otherwise = pbase * exp (negate _GMR * dh / tbase)
  where
    dh = h - baseHeight l
    pbase = baseSpecificPressure l
    tbase = baseTemperature l
    tlocal = abs $ layerTemperature l dh -- abs is irrelevant, as it will always be positive unless that layer is also masked out by its indicator function

windowedLayerDelta :: (Floating a) => Layer a -> Length a -> Dimensionless a
windowedLayerDelta l = \h -> indicator l h * layerDelta l h

windowedLayerTemperature :: (Floating a) => Layer a -> Length a -> ThermodynamicTemperature a
windowedLayerTemperature l = (\h -> let dh = h - baseHeight l
                                     in indicator l h * layerTemperature l dh)

indicator :: (Fractional a) => Layer a -> Length a -> Dimensionless a
indicator = uncurry window . heightLimits

-- geopotential altitude input
theta :: (Floating a) => Length a -> Dimensionless a
theta h = (/ _TZERO) $ sum $ fmap (($ h) . windowedLayerTemperature) layerTable

-- geopotential altitude input
delta :: (Floating a) => Length a -> Dimensionless a
delta h = sum $ fmap (($ h) . windowedLayerDelta) layerTable

atmosphere :: (Floating a) => Length a -> (Dimensionless a, Dimensionless a, Dimensionless a)
atmosphere alt = (s,d,t)
  where
    t = theta h
    d = delta h
    s = d/t
    h = alt*_REARTH/(alt+_REARTH) -- geometric to geopotential altitude

signum' :: (Fractional a, KnownDimension d) => Quantity d a -> Dimensionless a
signum' = (*~ one) . signum . (/~ siUnit)

-- | @window min max x@ is `_1` if @x@ is between @min@ and @max@, `_0` if it isn't, and @_1 / _2@ on
-- the exact boundaries.
window :: (Fractional a, KnownDimension d) => Maybe (Quantity d a) -> Maybe (Quantity d a) -> Quantity d a -> Dimensionless a
window Nothing    Nothing    _ = _1
window (Just bot) Nothing    x = (signum' (x - bot) + _1) / _2
window Nothing    (Just top) x = (_1 - signum' (x - top)) / _2
window (Just bot) (Just top) x = (signum' (x - bot) - signum' (x - top)) / _2
