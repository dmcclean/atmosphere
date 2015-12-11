{-# OPTIONS_GHC -Wall #-}

module Atmosphere
       ( Atmos(..)
       , siAtmosphere
       , siAtmosphereSmooth
       , usAtmosphere
       , atmosphere
       , atmosphereSmooth
       , siAltitudeFromPressure
       ) where

import Atmosphere.Constants
import Math.Polynomial

data Atmos a = Atmos { atmosTemperature :: a
                     , atmosPressure :: a
                     , atmosDensity :: a
                     , atmosSpeedOfSound :: a
                     , atmosViscosity :: a
                     , atmosKinematicViscosity :: a
                     }

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
siAtmosphere :: (Floating a, Ord a) => a -> Atmos a
siAtmosphere = siAtmosphere' atmosphere

siAtmosphereSmooth :: (Floating a, Eq a) => a -> Atmos a
siAtmosphereSmooth = siAtmosphere' atmosphereSmooth

siAtmosphere' :: (Floating a) => (a -> (a,a,a)) -> a -> Atmos a
siAtmosphere' f alt_m =
  Atmos { atmosTemperature = temp
        , atmosPressure = pressure
        , atmosDensity = density
        , atmosSpeedOfSound = asound
        , atmosViscosity = viscosity
        , atmosKinematicViscosity = kinematicViscosity
        }
  where
    alt_km = 0.001*alt_m
    (sigma, delta, theta) = f alt_km
    temp = _TZERO * theta
    pressure = _PZERO * delta
    density = _RHOZERO * sigma
    asound = _AZERO * sqrt theta
    viscosity = metricViscosity theta
    kinematicViscosity = viscosity/density

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
usAtmosphere :: (Floating a, Ord a) => a -> Atmos a
usAtmosphere alt_ft =
  Atmos { atmosTemperature = temp
        , atmosPressure = pressure
        , atmosDensity = density
        , atmosSpeedOfSound = asound
        , atmosViscosity = viscosity
        , atmosKinematicViscosity = kinematicViscosity
        }
  where
    alt_km = 0.001*_FT2METERS*alt_ft
    (sigma, delta, theta) = atmosphere alt_km
    temp = _KELVIN2RANKINE*_TZERO*theta
    pressure = _PZERO*delta/47.88
    density = _RHOZERO*sigma/515.379
    asound = (_AZERO/_FT2METERS)*sqrt theta
    viscosity=(1.0/_PSF2NSM)*metricViscosity theta
    kinematicViscosity = viscosity/density

{- |
   Compute altitude at which the standard atmosphere has a certain pressure.

   Input: Pressure, N/m^2

   Output: Altitude in meters
-}
siAltitudeFromPressure :: (Floating a, Ord a) => a -> a
siAltitudeFromPressure pressureIn = 1000*alt
  where
    alt = _REARTH / (_REARTH/h - 1)
    deltaIn = pressureIn / _PZERO

    (htabI, tbase, ptabI, tgradI) = getI htpgTable
      where
        getI [htab'] = htab'
        getI (htab0:htab1@(_,_,delta',_):htabs)
          | deltaIn < delta'    = getI (htab1:htabs)
          | otherwise = htab0
        getI [] = error "something went wrong"

    h
      | 0.0 == tgradI = htabI - tbase / _GMR * (log (deltaIn / ptabI))
      | otherwise     = htabI + tbase/tgradI*((deltaIn/ptabI)**(-tgradI/_GMR) - 1)

metricViscosity :: (Floating a) => a -> a
metricViscosity theta = _BETAVISC*sqrt(t*t*t)/(t+_SUTH)
  where
    t = theta * _TZERO

{- |
   Compute temperature, density, and pressure in standard atmosphere.

   Correct to 86 km.  Only approximate thereafter.

   Input: alt geometric altitude, km.

   Output: (sigma, delta, theta)

   > sigma - density/sea-level standard density
   > delta - pressure/sea-level standard pressure
   > theta - temperature/sea-level std. temperature
-}
atmosphere :: (Floating a, Ord a) => a -> (a,a,a)
atmosphere alt = (sigma, delta, theta)
  where
    h = alt*_REARTH/(alt+_REARTH) -- geometric to geopotential altitude

    (htabI, tbase, ptabI, tgradI) = getI htpgTable
      where
        getI [htab'] = htab'
        getI (htab0:htab1@(h',_,_,_):htabs)
          | h > h'    = getI (htab1:htabs)
          | otherwise = htab0
        getI [] = error "something went wrong"

    deltah = h - htabI              -- height above local base
    tlocal = tbase + tgradI*deltah  -- local temperature

    theta  =  tlocal/_TZERO    -- temperature ratio

    delta
      | 0.0 == tgradI = ptabI*exp(-_GMR*deltah/tbase)
      | otherwise     = ptabI*(tbase/tlocal)**(_GMR/tgradI)
    sigma = delta/theta

{-

A smooth approximation of the model.
delta is a 9th order polynomial approximation derived using R, max error 0.002
theta is a 9th order polynomial approximation derived using R, max error 0.02
-}
atmosphereSmooth :: (Floating a, Eq a) => a -> (a,a,a)
atmosphereSmooth h = (sigma, delta, theta)
  where
    sigma = delta / theta
    delta = deltaH h
    theta = thetaH h
    deltaH = evalPoly $ poly LE [  1.000089e+00
                                , -1.171864e-01
                                ,  5.081180e-03
                                , -5.244425e-05
                                , -3.613845e-06
                                ,  1.669225e-07
                                , -3.317844e-09
                                ,  3.595433e-11
                                , -2.068676e-13
                                ,  4.959128e-16
                                ]
    thetaH = evalPoly $ poly LE [  1.001541e+00
                                , -1.131267e-02
                                , -5.555655e-03
                                ,  8.056253e-04
                                , -4.953411e-05
                                ,  1.683472e-06
                                , -3.350431e-08
                                ,  3.870836e-10
                                , -2.401129e-12
                                ,  6.182736e-15
                                ]
