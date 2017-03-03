#!/usr/bin/env python

import model
import numpy as np

a = model.runZonedModel(mapFilename="arenas/arena2roomsLargeSmall.png", nbSteps=1000, dt=1./3., nbFishes=5, nbZones=9, seed=42, kappa0=np.random.random(9))
print a

