# the Gaussian case_study1 model set is numerically unchanged

    {
      "type": "list",
      "attributes": {
        "names": {
          "type": "character",
          "attributes": {},
          "value": ["formulas", "modname", "delta.AICc", "delta.BIC", "wi.AICc", "wi.BIC", "vi.aic", "vi.bic"]
        }
      },
      "value": [
        {
          "type": "character",
          "attributes": {
            "names": {
              "type": "character",
              "attributes": {},
              "value": ["ZONE", "ZONE+complexity", "ZONE+complexity.by.ZONE", "ZONE+depth", "ZONE+depth.by.ZONE", "complexity", "depth", "null"]
            }
          },
          "value": ["~ZONE + s(site, bs = \"re\")", "~s(complexity, k = 3, bs = \"cr\") + ZONE + s(site, bs = \"re\")", "~s(complexity, by = ZONE, k = 3, bs = \"cr\") + ZONE + s(site, bs = \"re\")", "~s(depth, k = 3, bs = \"cr\") + ZONE + s(site, bs = \"re\")", "~s(depth, by = ZONE, k = 3, bs = \"cr\") + ZONE + s(site, bs = \"re\")", "~s(complexity, k = 3, bs = \"cr\") + s(site, bs = \"re\")", "~s(depth, k = 3, bs = \"cr\") + s(site, bs = \"re\")", "~s(site, bs = \"re\")"]
        },
        {
          "type": "character",
          "attributes": {},
          "value": ["ZONE", "ZONE+complexity", "ZONE+complexity.by.ZONE", "ZONE+depth", "ZONE+depth.by.ZONE", "complexity", "depth", "null"]
        },
        {
          "type": "double",
          "attributes": {},
          "value": [56.006, 0.506, 2.46, 54.094, 54.307, 0, 53.126, 55.121]
        },
        {
          "type": "double",
          "attributes": {},
          "value": [52.716, 0.883, 4.255, 57.244, 58.833, 0, 54.766, 49.929]
        },
        {
          "type": "double",
          "attributes": {},
          "value": [0, 0.375, 0.141, 0, 0, 0.483, 0, 0]
        },
        {
          "type": "double",
          "attributes": {},
          "value": [0, 0.365, 0.068, 0, 0, 0.567, 0, 0]
        },
        {
          "type": "double",
          "attributes": {
            "names": {
              "type": "character",
              "attributes": {},
              "value": ["complexity", "depth", "ZONE"]
            }
          },
          "value": [0.999, 0, 0.516]
        },
        {
          "type": "double",
          "attributes": {
            "names": {
              "type": "character",
              "attributes": {},
              "value": ["complexity", "depth", "ZONE"]
            }
          },
          "value": [1, 0, 0.433]
        }
      ]
    }

---

    {
      "type": "list",
      "attributes": {
        "names": {
          "type": "character",
          "attributes": {},
          "value": ["AICc", "BIC", "r2.vals", "edf"]
        }
      },
      "value": [
        {
          "type": "double",
          "attributes": {},
          "value": [205.46740938, 149.96699883, 151.92066796, 203.55462235, 203.76762372, 149.46106127, 202.58724059, 204.58255929]
        },
        {
          "type": "double",
          "attributes": {},
          "value": [211.7509325, 159.91795797, 163.29002684, 216.27924471, 217.86849559, 159.0352027, 213.80148126, 208.96419211]
        },
        {
          "type": "double",
          "attributes": {},
          "value": [0.01891, 0.59327, 0.59267, 0.16155, 0.18273, 0.59361, 0.15327, 0.07198]
        },
        {
          "type": "double",
          "attributes": {},
          "value": [3, 4.67, 5.68, 5.45, 6.26, 3.7, 4.59, 2]
        }
      ]
    }

# variable importance under VI.mods = 'all' is numerically unchanged

    {
      "type": "list",
      "attributes": {
        "names": {
          "type": "character",
          "attributes": {},
          "value": ["aic", "bic"]
        }
      },
      "value": [
        {
          "type": "double",
          "attributes": {
            "names": {
              "type": "character",
              "attributes": {},
              "value": ["complexity", "depth", "ZONE"]
            }
          },
          "value": [0.999, 0, 0.516]
        },
        {
          "type": "double",
          "attributes": {
            "names": {
              "type": "character",
              "attributes": {},
              "value": ["complexity", "depth", "ZONE"]
            }
          },
          "value": [1, 0, 0.433]
        }
      ]
    }

# the negative binomial case_study1 model set is numerically unchanged

    {
      "type": "list",
      "attributes": {
        "names": {
          "type": "character",
          "attributes": {},
          "value": ["formulas", "modname", "delta.AICc", "delta.BIC", "wi.AICc", "wi.BIC", "vi.aic", "vi.bic"]
        }
      },
      "value": [
        {
          "type": "character",
          "attributes": {
            "names": {
              "type": "character",
              "attributes": {},
              "value": ["ZONE", "ZONE+complexity", "ZONE+complexity.by.ZONE", "ZONE+depth", "ZONE+depth.by.ZONE", "complexity", "depth", "null"]
            }
          },
          "value": ["~ZONE + s(site, bs = \"re\")", "~s(complexity, k = 3, bs = \"cr\") + ZONE + s(site, bs = \"re\")", "~s(complexity, by = ZONE, k = 3, bs = \"cr\") + ZONE + s(site, bs = \"re\")", "~s(depth, k = 3, bs = \"cr\") + ZONE + s(site, bs = \"re\")", "~s(depth, by = ZONE, k = 3, bs = \"cr\") + ZONE + s(site, bs = \"re\")", "~s(complexity, k = 3, bs = \"cr\") + s(site, bs = \"re\")", "~s(depth, k = 3, bs = \"cr\") + s(site, bs = \"re\")", "~s(site, bs = \"re\")"]
        },
        {
          "type": "character",
          "attributes": {},
          "value": ["ZONE", "ZONE+complexity", "ZONE+complexity.by.ZONE", "ZONE+depth", "ZONE+depth.by.ZONE", "complexity", "depth", "null"]
        },
        {
          "type": "double",
          "attributes": {},
          "value": [64.445, 0.585, 0, 67.505, 63.856, 0.115, 65.091, 61.561]
        },
        {
          "type": "double",
          "attributes": {},
          "value": [61.218, 1.379, 0, 68.794, 68.605, 0.625, 62.761, 54.691]
        },
        {
          "type": "double",
          "attributes": {},
          "value": [0, 0.277, 0.372, 0, 0, 0.351, 0, 0]
        },
        {
          "type": "double",
          "attributes": {},
          "value": [0, 0.225, 0.448, 0, 0, 0.328, 0, 0]
        },
        {
          "type": "double",
          "attributes": {
            "names": {
              "type": "character",
              "attributes": {},
              "value": ["complexity", "depth", "ZONE"]
            }
          },
          "value": [1, 0, 0.649]
        },
        {
          "type": "double",
          "attributes": {
            "names": {
              "type": "character",
              "attributes": {},
              "value": ["complexity", "depth", "ZONE"]
            }
          },
          "value": [1.001, 0, 0.673]
        }
      ]
    }

---

    {
      "type": "list",
      "attributes": {
        "names": {
          "type": "character",
          "attributes": {},
          "value": ["AICc", "BIC", "r2.vals", "edf"]
        }
      },
      "value": [
        {
          "type": "double",
          "attributes": {},
          "value": [607.61919372, 543.7588198, 543.17405447, 610.67871272, 607.03038197, 543.28883907, 608.26535981, 604.73491505]
        },
        {
          "type": "double",
          "attributes": {},
          "value": [617.18401635, 557.34482215, 555.96614679, 624.76013634, 624.57130594, 556.59146314, 618.72708705, 610.65682158]
        },
        {
          "type": "double",
          "attributes": {},
          "value": [0.12484, 0.53131, 0.5259, 0.15369, 0.21775, 0.53454, 0.14571, 0.13446]
        },
        {
          "type": "double",
          "attributes": {},
          "value": [3, 4.95, 5.79, 5.07, 7.41, 4.51, 3.24, 2]
        }
      ]
    }

# the binomial gamm4/uGamm model set is numerically unchanged

    {
      "type": "list",
      "attributes": {
        "names": {
          "type": "character",
          "attributes": {},
          "value": ["formulas", "modname", "delta.AICc", "delta.BIC", "wi.AICc", "wi.BIC", "vi.aic", "vi.bic"]
        }
      },
      "value": [
        {
          "type": "character",
          "attributes": {
            "names": {
              "type": "character",
              "attributes": {},
              "value": ["Depth", "Depth+bleach.pres", "Depth.by.bleach.pres+bleach.pres", "av.wave", "av.wave+bleach.pres", "av.wave.by.bleach.pres+bleach.pres", "bleach.pres", "null"]
            }
          },
          "value": ["~s(Depth, k = 4, bs = \"cr\")", "~s(Depth, k = 4, bs = \"cr\") + bleach.pres", "~s(Depth, by = bleach.pres, k = 4, bs = \"cr\") + bleach.pres", "~s(av.wave, k = 4, bs = \"cr\")", "~s(av.wave, k = 4, bs = \"cr\") + bleach.pres", "~s(av.wave, by = bleach.pres, k = 4, bs = \"cr\") + bleach.pres", "~bleach.pres", "~1"]
        },
        {
          "type": "character",
          "attributes": {},
          "value": ["Depth", "Depth+bleach.pres", "Depth.by.bleach.pres+bleach.pres", "av.wave", "av.wave+bleach.pres", "av.wave.by.bleach.pres+bleach.pres", "bleach.pres", "null"]
        },
        {
          "type": "double",
          "attributes": {},
          "value": [84.731, 24.849, 27.963, 66.68, 17.656, 0, 22.36, 82.341]
        },
        {
          "type": "double",
          "attributes": {},
          "value": [73.368, 17.287, 27.963, 55.318, 10.093, 0, 7.186, 63.344]
        },
        {
          "type": "double",
          "attributes": {},
          "value": [0, 0, 0, 0, 0, 1, 0, 0]
        },
        {
          "type": "double",
          "attributes": {},
          "value": [0, 0, 0, 0, 0.006, 0.967, 0.027, 0]
        },
        {
          "type": "double",
          "attributes": {
            "names": {
              "type": "character",
              "attributes": {},
              "value": ["av.wave", "Depth", "bleach.pres"]
            }
          },
          "value": [1, 0, 1]
        },
        {
          "type": "double",
          "attributes": {
            "names": {
              "type": "character",
              "attributes": {},
              "value": ["av.wave", "Depth", "bleach.pres"]
            }
          },
          "value": [0.973, 0, 1]
        }
      ]
    }

---

    {
      "type": "list",
      "attributes": {
        "names": {
          "type": "character",
          "attributes": {},
          "value": ["AICc", "BIC", "r2.vals"]
        }
      },
      "value": [
        {
          "type": "double",
          "attributes": {},
          "value": [6042.03795297, 5982.15670182, 5985.26983452, 6023.98732938, 5974.96354128, 5957.30730046, 5979.66730171, 6039.6485183]
        },
        {
          "type": "double",
          "attributes": {},
          "value": [6057.35374356, 6001.27194899, 6011.94788123, 6039.30311997, 5994.07878845, 5983.98534716, 5991.17173701, 6047.32980248]
        },
        {
          "type": "double",
          "attributes": {},
          "value": [0.08436, 0.09093, 0.09335, 0.00781, 0.00171, 0.00663, 0.01085, 0]
        }
      ]
    }

# the cyclic case_study3 model set is numerically unchanged

    {
      "type": "list",
      "attributes": {
        "names": {
          "type": "character",
          "attributes": {},
          "value": ["formulas", "modname", "delta.AICc", "delta.BIC", "wi.AICc", "wi.BIC", "vi.aic", "vi.bic"]
        }
      },
      "value": [
        {
          "type": "character",
          "attributes": {
            "names": {
              "type": "character",
              "attributes": {},
              "value": ["Sex", "Sex+lunar.date", "Sex+lunar.date.by.Sex", "Sex+month", "Sex+month.by.Sex", "lunar.date", "lunar.date+month", "month", "null"]
            }
          },
          "value": ["~Sex", "~s(lunar.date, k = 5, bs = \"cc\") + Sex", "~s(lunar.date, by = Sex, k = 5, bs = \"cc\") + Sex", "~s(month, k = 5, bs = \"cc\") + Sex", "~s(month, by = Sex, k = 5, bs = \"cc\") + Sex", "~s(lunar.date, k = 5, bs = \"cc\")", "~s(lunar.date, k = 5, bs = \"cc\") + s(month, k = 5, bs = \"cc\")", "~s(month, k = 5, bs = \"cc\")", "~1"]
        },
        {
          "type": "character",
          "attributes": {},
          "value": ["Sex", "Sex+lunar.date", "Sex+lunar.date.by.Sex", "Sex+month", "Sex+month.by.Sex", "lunar.date", "lunar.date+month", "month", "null"]
        },
        {
          "type": "double",
          "attributes": {},
          "value": [32.435, 0, 3.527, 27.502, 30.875, 72.64, 61.283, 98.56, 104.945]
        },
        {
          "type": "double",
          "attributes": {},
          "value": [17.954, 0, 16.45, 26.806, 35.744, 67.644, 70.418, 92.599, 85.463]
        },
        {
          "type": "double",
          "attributes": {},
          "value": [0, 0.854, 0.146, 0, 0, 0, 0, 0, 0]
        },
        {
          "type": "double",
          "attributes": {},
          "value": [0, 1, 0, 0, 0, 0, 0, 0, 0]
        },
        {
          "type": "double",
          "attributes": {
            "names": {
              "type": "character",
              "attributes": {},
              "value": ["lunar.date", "month", "Sex"]
            }
          },
          "value": [1, 0, 1]
        },
        {
          "type": "double",
          "attributes": {
            "names": {
              "type": "character",
              "attributes": {},
              "value": ["lunar.date", "month", "Sex"]
            }
          },
          "value": [1, 0, 1]
        }
      ]
    }

---

    {
      "type": "list",
      "attributes": {
        "names": {
          "type": "character",
          "attributes": {},
          "value": ["AICc", "BIC", "r2.vals", "edf"]
        }
      },
      "value": [
        {
          "type": "double",
          "attributes": {},
          "value": [7527.90798135, 7495.47297697, 7499.00030583, 7522.97517297, 7526.34820251, 7568.11321644, 7556.75602209, 7594.03328521, 7600.41835192]
        },
        {
          "type": "double",
          "attributes": {},
          "value": [7542.92262741, 7524.9688947, 7541.41914216, 7551.77444232, 7560.71318409, 7592.61284172, 7595.38718878, 7617.56805332, 7610.4317424]
        },
        {
          "type": "double",
          "attributes": {},
          "value": [0.06823, 0.10237, 0.10061, 0.07891, 0.07596, 0.03517, 0.05183, 0.01312, 0]
        },
        {
          "type": "double",
          "attributes": {},
          "value": [2, 4.9, 7.49, 4.76, 5.88, 3.9, 6.73, 3.71, 1]
        }
      ]
    }

